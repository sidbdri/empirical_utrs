#!/usr/bin/env python

"""Usage:
    get_utr_sequences [--log-level=<log-level>] <input-file> <output-file>

-h --help                    Show this message.
-v --version                 Show version.
--log-level=<log-level>      Set logging level (one of {log_level_vals})
                             [default: info].
<input-file>                 Result from get_empirical_utrs.py
<output-file>                Output file containing the sequence of the utrs

Output 5' UTR sequences corresponding to the output of get_empirical_utrs.py.
"""

import sys
import xml.etree.ElementTree as etree
import urllib2
import docopt
import schema
import math

from . import log
from . import options as opt
from .__init__ import __version__

LOG_LEVEL = "--log-level"
LOG_LEVEL_VALS = str(log.LOG_LEVELS.keys())

INPUT_FILE = "<input-file>"
OUTPUT_FILE = "<output-file>"

FASTA_SEQ_LINE_LIMIT = 60
UCSC_DAS_BASE_PATH = "http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment="


def validate_command_line_options(options):
    try:
        opt.validate_dict_option(
            options[LOG_LEVEL], log.LOG_LEVELS, "Invalid log level")
        opt.validate_file_option(
            options[INPUT_FILE], "Input file must exist")
        opt.validate_file_already_exist_option(
            options[OUTPUT_FILE], "Output file must not already exist")
    except schema.SchemaError as exc:
        exit(exc.code)


def _get_sequence(chromosome, start, end, logger):
    url = UCSC_DAS_BASE_PATH + chromosome + ':' + str(start) + ',' + str(end)
    logger.debug("querying: {url}".format(url=url))

    xml = urllib2.urlopen(url).read()
    if xml != '':
        dom = etree.fromstring(xml)
        sequence = dom[0][0].text.replace('\n','')
    else:
        sequence = 'THE SEQUENCE DOES NOT EXIST FOR GIVEN COORDINATES'

    logger.debug("Query result from ucsc: {result}\n".format(result=sequence))
    return sequence


def _get_reverse_sequence(seq):
    complement_table = {"A":"T", "T":"A", "C":"G", "G":"C"}
    complement = "".join([complement_table.get(base, "N") for base in seq.upper()])
    return complement[::-1]


def _write_sequences(input_file, output_file, logger):
    out = open(output_file, "w+")

    with open(input_file, "r") as ip:
        #skip the first line as header
        next(ip)

        count = 0
        for line in ip:
            count += 1
            if count % 1000 == 0:
                logger.info("Processed {utr} UTRs.".format(utr=count))

            gene, transcript, chrom, strand, tss, pileup, utr = \
                line.strip('\n').split(',')
            tss = int(tss)
            utr = int(utr)

            logger.debug("Finding UTR sequence for gene {gene}: {line}".
                         format(gene=gene,line=line))

            if strand == '-':
                sequence = _get_sequence('chr' + chrom, tss - utr + 1, tss, logger)
                sequence = _get_reverse_sequence(sequence)
            else:
                sequence = _get_sequence('chr' + chrom, tss, tss + utr - 1, logger)

            chunks = math.ceil(float(len(sequence))/FASTA_SEQ_LINE_LIMIT)
            sequence_str = '\n'.join([sequence[i:i + FASTA_SEQ_LINE_LIMIT] \
                    for i in range(0, len(sequence), FASTA_SEQ_LINE_LIMIT)])

            out.write((">{gene_name},{transcript},{chrom},{strand},{tss},{pileup}," +
                       "{utr}\n{seq}\n").format(
                gene_name=gene,
                chrom=chrom, strand=strand, tss=tss, pileup=pileup,
                utr=utr, transcript=transcript, seq=sequence_str.upper()))


def get_utr_sequences(args):
    docstring = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(docstring,
                            version="get_utr_sequences v" + __version__)
    validate_command_line_options(options)

    logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

    _write_sequences(options[INPUT_FILE],options[OUTPUT_FILE],logger)
