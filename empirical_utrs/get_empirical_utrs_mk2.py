#!/usr/bin/env python

"""Usage:
    get_empirical_utrs [--log-level=<log-level>] [--max-utr-length=<max-utr-length>] <transcript-gtf-file> <cage-bam-file> <output-file>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals})
                                    [default: info].
--max-utr-length=<max-utr-length>   Maximum allowed length for empirical UTR
                                    [default: 1000].
<transcript-gtf-file>               File containing transcript definitions
                                    in GTF format.
<cage-bam-file>                     BAM file containing CAGE data.
<output-file>                       Output file container gene results

For each transcript of each gene, determine the region to be scanned for the
empirical transcription start site, namely from <max-utr-length> bases upstream
of the coding start of the transcript (if the transcript has no coding start
location defined, use the location of the first exon). Within these scan bounds
find the location with the greatest pile-up of reads mapped to the same strand
as the gene. In the case of ties, choose the rightmost location for genes on
the plus srand, and the leftmost location for genes on the minus strand.
"""

import docopt
import pysam
import schema
import sys

from collections import defaultdict

from . import gtf
from . import log
from . import options as opt
from .__init__ import __version__

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LOG_LEVELS.keys())

TRANSCRIPT_GTF_FILE = "<transcript-gtf-file>"
CAGE_BAM_FILE = "<cage-bam-file>"
OUTPUT_FILE = "<output-file>"
MAXIMUM_UTR_LENGTH = "--max-utr-length"

def validate_command_line_options(options):
    try:
        opt.validate_dict_option(
            options[LOG_LEVEL], log.LOG_LEVELS, "Invalid log level")
        opt.validate_file_option(
            options[TRANSCRIPT_GTF_FILE], "Transcript GTF file must exist")
        opt.validate_file_option(
            options[CAGE_BAM_FILE], "CAGE BAM file must exist")
        opt.validate_file_already_exist_option(
            options[OUTPUT_FILE], "Output file must not already exist")
        options[MAXIMUM_UTR_LENGTH] = opt.validate_int_option(
            options[MAXIMUM_UTR_LENGTH], "Maximum UTR length must be " +
            "greater than 0", 1)
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_tss_scan_bounds(gene, transcript, maximum_utr_length, logger):
    if transcript.coding_start is None:
        transcript_start = sys.maxint
        transcript_end = -sys.maxint - 1

        for exon in transcript.exons:
            if exon.start < transcript_start:
                transcript_start = exon.start
            if exon.end > transcript_end:
                transcript_end = exon.end

        transcript.coding_start = transcript_end if gene.strand == "-" \
            else transcript_start

        logger.debug("...no coding start, use the start position of the " +
                     "first exon {start} as coding start".format(
                         start=transcript.coding_start))

    if gene.strand == "-":
        return (transcript.coding_start, transcript.coding_start + maximum_utr_length + 500)

    return (transcript.coding_start - maximum_utr_length - 500, transcript.coding_start)


def _get_maximum_pileup_location(samfile, gene, scan_bounds, logger):
    logger.debug(("Get maximum pile up location, scan start {start}, " +
                  "end {end}, strand {strand}").
            format(start=scan_bounds[0], end=scan_bounds[1], strand=gene.strand))

    gene_is_reverse = gene.strand == "-"
    max_pileup = 0
    pileup_location = None

    for pileupcolumn in samfile.pileup(
            "chr" + str(gene.seqname), scan_bounds[0], scan_bounds[1],
            truncate=True):

        pileup_count = 0
        for pileupread in pileupcolumn.pileups:
            if pileupread.alignment.is_reverse == gene_is_reverse:
                pileup_count = pileup_count + 1

        if pileup_count > max_pileup or \
            (pileup_count == max_pileup and not gene_is_reverse):
            pileup_location = pileupcolumn.pos
            max_pileup = pileup_count

    return (max_pileup, pileup_location)


def _calculate_empirical_utrs(transcript_info, cage_bam, maximum_utr_length, logger):
    logger.info("Calculating empirical UTRs...")
    empirical_utrs = {}
    count = 0
    no_alignment_data = 0
    no_pileup_location = 0
    no_shortest_utr = 0
    utr_too_long = 0

    samfile = pysam.AlignmentFile(cage_bam)

    for gene_name, gene in transcript_info.iteritems():
        logger.debug("Calculating empirical UTR for {gene}.".
                format(gene=gene.name))

        for transcript in gene.transcripts.values():
            logger.debug("Calculating empirical UTR of {transcript}.".
                format(transcript=transcript.name))

            count = count + 1
            if count % 1000 == 0:
                logger.info("...processed {t} transcripts, {e} empirical UTRs found.".
                        format(t=count, e=len(empirical_utrs)))
                logger.info(("Found empirical UTRs for {num_transcripts} transcripts; " +
                             "no alignment data {no_align}, " +
                             "no pileup location {no_pileup}, " +
                             "no shortest UTR transcript {no_utr}, " +
                             "UTR too long {too_long}.").
                            format(num_transcripts=len(empirical_utrs),
                                   no_align=no_alignment_data,
                                   no_pileup=no_pileup_location,
                                   no_utr=no_shortest_utr,
                                   too_long=utr_too_long))

            chr_seqname = 'chr' + str(gene.seqname)
            if chr_seqname not in samfile.references:
                no_alignment_data += 1
                logger.debug("Found no alignment data for {transcript}".
                            format(transcript=transcript.name))
                continue

            scan_bounds = _get_tss_scan_bounds(
                    gene, transcript, maximum_utr_length, logger)

            max_pileup, pileup_location = _get_maximum_pileup_location(
                    samfile, gene, scan_bounds, logger)
            logger.debug("Pileup location {location}, max pileup {max_pileup}".
                    format(location=pileup_location, max_pileup=max_pileup))

            if pileup_location is None:
                no_pileup_location += 1
                logger.debug("Found no pileup location for {transcript}".
                            format(transcript=transcript.name))
                continue

            utr = pileup_location - transcript.coding_start
            logger.debug(("Transcript {transcript}, coding start {start}, " +
                          "UTR {utr}, strand {strand}").format(
                transcript=transcript.name, start=transcript.coding_start,
                utr=utr, strand=gene.strand))

            if utr < 0:
                utr = -utr

            if utr == 0:
                no_pileup_location += 1
                logger.debug("...pileup location overlaps with coding start. " +
                             "Skip this transcript.")
                continue

            if utr > maximum_utr_length:
                utr_too_long +=  1
                logger.debug("UTR too long {utr} for {transcript}".
                        format(utr=utr, transcript=transcript.name))
                continue

            empirical_utrs[transcript.name] = (gene_name, gene.seqname, gene.strand,
                                               pileup_location, max_pileup, utr)

    logger.info(("Found empirical UTRs for {num_transcripts} transcripts; " +
                 "no alignment data {no_align}, " +
                 "no pileup location {no_pileup}, " +
                 "no shortest UTR transcript {no_utr}, " +
                 "UTR too long {too_long}.").
                format(num_transcripts=len(empirical_utrs),
                       no_align=no_alignment_data,
                       no_pileup=no_pileup_location,
                       no_utr=no_shortest_utr,
                       too_long=utr_too_long))

    return empirical_utrs


def _print_empirical_utrs(empirical_utrs, logger, output_file):
    logger.info("Printing empirical UTRs...")

    f= open(output_file,"w+")
    f.write("gene,transcript,chr,strand,tss,pileup_depth,utr_length\n")

    for transcript in sorted(empirical_utrs.keys()):
        empirical_utr = empirical_utrs[transcript]

        f.write("{g},{t},{chr},{strand},{tss},{depth},{utr}\n".format(
            g=empirical_utr[0], t=transcript,
            chr=empirical_utr[1], strand=empirical_utr[2],
            tss=empirical_utr[3], depth=empirical_utr[4],
            utr=empirical_utr[5]))


def get_empirical_utrs(args):
    docstring = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(docstring,
                            version="get_empirical_utrs v" + __version__)
    validate_command_line_options(options)

    logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

    gtf_info = gtf.GtfInfo(options[TRANSCRIPT_GTF_FILE], logger)
    transcript_info = gtf_info.get_transcript_info()

    empirical_utrs = _calculate_empirical_utrs(
        transcript_info, options[CAGE_BAM_FILE],
        options[MAXIMUM_UTR_LENGTH], logger)

    _print_empirical_utrs(empirical_utrs, logger, options[OUTPUT_FILE])
