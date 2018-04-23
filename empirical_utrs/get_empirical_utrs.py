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


def _calculate_empirical_utrs(transcript_info, cage_bam, logger):
    logger.info("Calculating empirical UTRs...")
    empirical_utrs = {}

    for gene_name, gene in transcript_info.iteritems():
        gene_start = sys.maxint
        gene_end = -sys.maxint - 1
        exon_starts = defaultdict(list)
        exon_ends = defaultdict(list)

        for transcript in gene.transcripts.values():
            for exon in transcript.exons:
                if exon.start < gene_start:
                    gene_start = exon.start
                if exon.end > gene_end:
                    gene_end = exon.end

                exon_starts[exon.start].append(exon)
                exon_ends[exon.end].append(exon)

        logger.info("Gene {gene}, chromosome {chr}, start {start}, end {end}".
                format(gene=gene, chr=gene.seqname,
                       start=gene_start, end=gene_end))

        empirical_utrs[gene_name] = 0

        if len(empirical_utrs) % 1000 == 0:
            logger.info("...processed {g} genes.".format(g=len(gene_lengths)))

    return empirical_utrs


def _print_empirical_utrs(empirical_utrs, logger):
    logger.info("Printing empirical UTRs...")

    print("gene,...")

    for gene in sorted(empirical_utrs.keys()):
        empirical_utr = empirical_utrs[gene]

        print("{g},{eu}".format(
            g=gene, eu=empirical_utr))


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
