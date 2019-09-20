#!/usr/bin/env python

"""Usage:
    get_empirical_utrs [--log-level=<log-level>] [--max-utr-length=<max-utr-length>] <transcript-gtf-file> <cage-bam-file> <output-file>

-h --help                           Show this message.
-v --version                        Show version.
--log-level=<log-level>             Set logging level (one of {log_level_vals})
                                    [default: info].
--max-utr-length=<max-utr-length>   Maximum allowed length for empirical UTR
                                    [default: 500].
<transcript-gtf-file>               File containing transcript definitions
                                    in GTF format.
<cage-bam-file>                     BAM file containing CAGE data.
<output-file>                       Output file container gene results

Calculate empirical UTRs via the following procedure:

    i) For each gene in the supplied GTF file, determine the region to be
    scanned for the empirical transcription start site, namely from 1000 bases
    upstream of the earliest location of the start of any transcript, to the
    most downstream location of the end of any transcript.

    ii) Within these scan bounds find location with greatest pile up of reads
    mapped to the same strand as the gene. In the case of ties, choose the
    rightmost location for genes on the plus strand, and leftmost location for
    genes on the minus strand. This location is defined to be the empirical TSS
    for the gene.

    iii) Choose the transcript for which the empirical TSS is the smallest
    number of bases upstream (> 0) of the defined coding start location. If the
    TSS is downstream of the coding start location, ignore this transcript. If
    the transcript has no coding start location defined, use the location of
    the start of the first exon. The empirical 5' UTR is defined as the region
    between the empirical TSS and the coding start location of the chosen
    transcript.

    iv) Reject this empirical 5' UTR if it is longer than <max-utr-length> bases.
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
            "greater than 1", 1)
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_gene_bounds(gene, logger):
    gene_start = sys.maxint
    gene_end = -sys.maxint - 1

    for transcript in gene.transcripts.values():
        for exon in transcript.exons:
            if exon.start < gene_start:
                gene_start = exon.start
            if exon.end > gene_end:
                gene_end = exon.end

    logger.debug(("Gene {gene}, chromosome {chr}, start {start}, " +
                  "end {end}, strand {strand}").
            format(gene=gene.name, chr=gene.seqname,
                   start=gene_start, end=gene_end,
                   strand=gene.strand))

    return(gene_start, gene_end)


def _get_tss_scan_bounds(gene, gene_bounds, logger):
    return (gene_bounds[0], gene_bounds[1] + 1000) if gene.strand == "-" \
        else (gene_bounds[0] - 1000, gene_bounds[1])


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


def _get_shortest_utr(gene, gene_bounds, pileup_location, logger):
    logger.debug("Getting shortest UTR for {gene}, pileup {pileup_location}".
            format(gene=gene.name, pileup_location=pileup_location))

    shortest_utr = sys.maxint
    shortest_utr_transcript = None

    for transcript in gene.transcripts.values():
        logger.debug("Transcript {transcript}, coding start {start}".format(
            transcript=transcript.name, start=transcript.coding_start))

        if transcript.coding_start is None:
            transcript.coding_start = gene_bounds[1] if gene.strand == "-" \
                else gene_bounds[0]

            logger.debug("...no coding start, use the start position of the " +
                         "first exon {start} as coding start".format(
                             start=transcript.coding_start))

        utr = pileup_location - transcript.coding_start
        logger.debug(("Transcript {transcript}, coding start {start}, " +
                      "UTR {utr}, strand {strand}").format(
            transcript=transcript.name, start=transcript.coding_start,
            utr=utr, strand=gene.strand))

        if (gene.strand == "+" and utr > 0) or (gene.strand == "-" and utr < 0):
            logger.debug("...pileup location inconsistent with coding start")
            continue

        if utr < 0:
            utr = -utr

        if utr == 0:
            logger.debug("...pileup location overlaps with coding start. " +
                         "Skip this transcript.")
            continue

        if utr < shortest_utr:
            logger.debug("...new shortest UTR {utr}".format(utr=utr))
            shortest_utr = utr
            shortest_utr_transcript = transcript

    logger.debug("Found shortest UTR {utr} for {transcript}".
                format(utr=shortest_utr,
                       transcript="no transcript" \
                               if shortest_utr_transcript is None else \
                               shortest_utr_transcript.name))

    return (shortest_utr, shortest_utr_transcript)


def _calculate_empirical_utrs(transcript_info, cage_bam, maximum_utr_length, logger):
    logger.info("Calculating empirical UTRs...")
    empirical_utrs = {}
    count = 0
    no_alignment_data = 0
    no_pileup_location = 0
    no_shortest_utr = 0
    shortest_utr_too_long = 0

    samfile = pysam.AlignmentFile(cage_bam)

    for gene_name, gene in transcript_info.iteritems():
        logger.debug("Calculating empirical UTR for {gene}.".
                format(gene=gene.name))

        count = count + 1
        if count % 1000 == 0:
            logger.info("...processed {g} genes, {e} empirical UTRs found.".
                    format(g=count, e=len(empirical_utrs)))
            logger.info(("Found empirical UTRs for {num_genes} genes; " +
                         "no alignment data {no_align}, " +
                         "no pileup location {no_pileup}, " +
                         "no shortest UTR transcript {no_utr}, " +
                         "shortest UTR too long {too_long}.").
                        format(num_genes=len(empirical_utrs),
                               no_align=no_alignment_data,
                               no_pileup=no_pileup_location,
                               no_utr=no_shortest_utr,
                               too_long=shortest_utr_too_long))

        chr_seqname = 'chr' + str(gene.seqname)
        if chr_seqname not in samfile.references:
            no_alignment_data += 1
            logger.debug("Found no alignment data for {gene}".
                        format(gene=gene_name))
            continue

        gene_bounds = _get_gene_bounds(gene, logger)
        scan_bounds = _get_tss_scan_bounds(gene, gene_bounds, logger)

        max_pileup, pileup_location = _get_maximum_pileup_location(
                samfile, gene, scan_bounds, logger)
        logger.debug("Pileup location {location}, max pileup {max_pileup}".
                format(location=pileup_location, max_pileup=max_pileup))

        if pileup_location is None:
            no_pileup_location += 1
            logger.debug("Found no pileup location for {gene}".
                        format(gene=gene_name))
            continue

        shortest_utr, shortest_utr_transcript = _get_shortest_utr(
                gene, gene_bounds, pileup_location, logger)

        if shortest_utr_transcript is None:
            no_shortest_utr += 1
            logger.debug("Found no shortest UTR transcript for {gene}".
                        format(gene=gene_name))
            continue

        if shortest_utr > maximum_utr_length:
            shortest_utr_too_long +=  1
            logger.debug("Shorted UTR too long {utr} for {gene}".
                    format(utr=shortest_utr, gene=gene_name))
            continue

        empirical_utrs[gene_name] = (pileup_location, max_pileup,
                                     shortest_utr, shortest_utr_transcript.name,
                                     gene.strand, gene.seqname)

    logger.info(("Found empirical UTRs for {num_genes} genes; " +
                 "no alignment data {no_align}, " +
                 "no pileup location {no_pileup}, " +
                 "no shortest UTR transcript {no_utr}, " +
                 "shortest UTR too long {too_long}.").
                format(num_genes=len(empirical_utrs),
                       no_align=no_alignment_data,
                       no_pileup=no_pileup_location,
                       no_utr=no_shortest_utr,
                       too_long=shortest_utr_too_long))

    return empirical_utrs


def _print_empirical_utrs(empirical_utrs, logger, output_file):
    logger.info("Printing empirical UTRs...")

    f= open(output_file,"w+")
    f.write("gene,chr,strand,tss,pileup,utr,transcript\n")

    for gene in sorted(empirical_utrs.keys()):
        empirical_utr = empirical_utrs[gene]

        f.write("{g},{chr},{strand},{tss},{pileup},{utr},{transcript}\n".format(
            g=gene,
            tss=empirical_utr[0], pileup=empirical_utr[1],
            utr=empirical_utr[2], transcript=empirical_utr[3],
            strand=empirical_utr[4], chr=empirical_utr[5]))


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
