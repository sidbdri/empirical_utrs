#!/usr/bin/env python

"""Usage:
    get_empirical_utrs [--log-level=<log-level>] <transcript-gtf-file> <cage-bam-file>

-h --help                    Show this message.
-v --version                 Show version.
--log-level=<log-level>      Set logging level (one of {log_level_vals})
                             [default: info].
<transcript-gtf-file>        File containing transcript definitions in GTF
                             format.
<cage-bam-file>              BAM file containing CAGE data.

TODO: Calculate empirical UTRs.
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

def validate_command_line_options(options):
    try:
        opt.validate_dict_option(
            options[LOG_LEVEL], log.LOG_LEVELS, "Invalid log level")
        opt.validate_file_option(
            options[TRANSCRIPT_GTF_FILE], "Transcript GTF file must exist")
        opt.validate_file_option(
            options[CAGE_BAM_FILE], "CAGE BAM file must exist")
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_tss_scan_bounds(gene, logger):
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

    # THIS IS WRONG
    return(gene_start - 1000, gene_end)


def _get_maximum_pileup_location(samfile, gene, scan_start, scan_end, logger):
    gene_is_reverse = gene.strand == "-"
    max_pileup = 0
    pileup_location = None

    for pileupcolumn in samfile.pileup(
            "chr" + str(gene.seqname), scan_start, scan_end, truncate=True):

        pileup_count = 0
        for pileupread in pileupcolumn.pileups:
            if pileupread.alignment.is_reverse == gene_is_reverse:
                pileup_count = pileup_count + 1

        if pileup_count > max_pileup or \
            (pileup_count == max_pileup and not gene_is_reverse):
            pileup_location = pileupcolumn.pos
            max_pileup = pileup_count

    return (max_pileup, pileup_location)


def _get_shortest_utr(gene, pileup_location, logger):
    logger.debug("Getting shortest UTR for {gene}".format(gene=gene.name))

    shortest_utr = sys.maxint
    shortest_utr_transcript = None

    for transcript in gene.transcripts.values():
        logger.debug("Transcript {transcript}, coding start {start}".format(
            transcript=transcript.name, start=transcript.coding_start))

        if transcript.coding_start is None:
            logger.debug("...no coding start")
            continue

        utr = pileup_location - transcript.coding_start
        logger.debug("Transcript {transcript}, coding start {start}, UTR {utr}, strand {strand}".format(
            transcript=transcript.name, start=transcript.coding_start,
            utr=utr, strand=gene.strand))

        if (gene.strand == "+" and utr > 0) or (gene.strand == "-" and utr < 0):
            logger.debug("...pileup location inconsistent with coding start")
            continue

        if utr < 0:
            utr = -utr

        if utr < shortest_utr:
            logger.debug("...new shortest UTR {utr}".format(utr=utr))
            shortest_utr = utr
            shortest_utr_transcript = transcript

    if shortest_utr_transcript is None or shortest_utr > 500:
        logger.debug("No shortest UTR found.")
        return (None, None)

    logger.debug("Found shortest UTR {utr} for {transcript}".
                format(utr=shortest_utr, transcript=shortest_utr_transcript.name))
    return (shortest_utr, shortest_utr_transcript)


def _calculate_empirical_utrs(transcript_info, cage_bam, logger):
    logger.info("Calculating empirical UTRs...")
    empirical_utrs = {}
    count = 0
    no_alignment_data = 0
    no_pileup_location = 0
    no_shortest_utr = 0

    samfile = pysam.AlignmentFile(cage_bam)

    for gene_name, gene in transcript_info.iteritems():
        count = count + 1
        if count % 1000 == 0:
            logger.info("...processed {g} genes.".format(g=len(empirical_utrs)))

        chr_seqname = 'chr' + str(gene.seqname)
        if chr_seqname not in samfile.references:
            no_alignment_data += 1
            logger.debug("Found no alignment data for {gene}".
                        format(gene=gene_name))
            continue

        scan_start, scan_end = _get_tss_scan_bounds(gene, logger)

        max_pileup, pileup_location = _get_maximum_pileup_location(
                samfile, gene, scan_start, scan_end, logger)

        if pileup_location is None:
            no_pileup_location += 1
            logger.debug("Found no pileup location for {gene}".
                        format(gene=gene_name))
            continue

        shortest_utr, shortest_utr_transcript = _get_shortest_utr(
                gene, pileup_location, logger)

        if shortest_utr_transcript is None:
            no_shortest_utr += 1
            logger.debug("Found no shortest UTR for {gene}".
                        format(gene=gene_name))
            continue

        empirical_utrs[gene_name] = (pileup_location, max_pileup,
                                     shortest_utr, shortest_utr_transcript.name)


    logger.info("Found empirical UTRs for {num_genes} genes; no alignment data {no_align}, no pileup location {no_pileup}, no shortest UTR {no_utr}.".
                format(num_genes=len(empirical_utrs),
                       no_align=no_alignment_data,
                       no_pileup=no_pileup_location,
                       no_utr=no_shortest_utr))

    return empirical_utrs


def _print_empirical_utrs(empirical_utrs, logger):
    logger.info("Printing empirical UTRs...")

    print("gene,tss,pileup,utr,transcript")

    for gene in sorted(empirical_utrs.keys()):
        empirical_utr = empirical_utrs[gene]

        print("{g},{tss},{pileup},{utr},{transcript}".format(
            g=gene,
            tss=empirical_utr[0], pileup=empirical_utr[1],
            utr=empirical_utr[2], transcript=empirical_utr[3]))


def get_empirical_utrs(args):
    docstring = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(docstring,
                            version="get_empirical_utrs v" + __version__)
    validate_command_line_options(options)

    logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

    gtf_info = gtf.GtfInfo(options[TRANSCRIPT_GTF_FILE], logger)
    transcript_info = gtf_info.get_transcript_info()

    empirical_utrs = _calculate_empirical_utrs(
        transcript_info, options[CAGE_BAM_FILE], logger)

    _print_empirical_utrs(empirical_utrs, logger)
