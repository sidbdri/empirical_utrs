#!/usr/bin/env python
# http://biodas.open-bio.org/documents/spec-1.53.html


"""Usage:
    get_sequence [--log-level=<log-level>] <input-file> <output-file>

-h --help                    Show this message.
-v --version                 Show version.
--log-level=<log-level>      Set logging level (one of {log_level_vals})
                             [default: info].
<input-file>                 Result from get_empirical_utrs.py
<output-file>                Output file containing the sequence of the utrs

TODO: Calculate empirical UTRs.
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

FASTA_SEQ_LINE_LIMIT=80
UCSC_DAS_BASE_PATH="http://genome.ucsc.edu/cgi-bin/das/mm10/dna?segment="

def validate_command_line_options(options):
    try:
        opt.validate_dict_option(
            options[LOG_LEVEL], log.LOG_LEVELS, "Invalid log level")
        opt.validate_file_option(
            options[INPUT_FILE], "Input file must exist")
        opt.validate_file_already_exist_option(
            options[OUTPUT_FILE], "Output file must not exist")
    except schema.SchemaError as exc:
        exit(exc.code)


def _getSequence(chromosome, start, end,logger):
    base = UCSC_DAS_BASE_PATH
    url = base + chromosome + ':' + str(start) + ',' + str(end)
    logger.debug("querying: {url}".format(url=url))

    xml = urllib2.urlopen(url).read()
    if xml != '':
        dom =  etree.fromstring(xml)
        sequence = dom[0][0].text.replace('\n','')
    else:
        sequence = 'THE SEQUENCE DOES NOT EXIST FOR GIVEN COORDINATES'

    logger.debug("Query result from ucsc: {result}\n".format(result=sequence))
    return sequence


def _getReverseSeq(seq):
    complement_table = {"A":"T", "T":"A", "C":"G", "G":"C"}
    complement = "".join([complement_table.get(base, "N") for base in seq.upper()])
    return complement[::-1]



def _write_result(input,output,logger):
    out=open(output,"w+")

    with open(input, "r") as ip:
        #skip the first line as header
        next(ip)
        for line in ip:
            line=line.strip('\n')
            gene,chr,strand,tss,pileup,utr,transcript = line.split(',')
            logger.info("Finding UTR sequence for gene {gene}: {line}".format(gene=gene,line=line))
            if strand=='-':
                sequence = _getSequence('chr' + str(chr),int(tss)-int(utr),int(tss),logger)
                sequence = _getReverseSeq(sequence)
            else:
                sequence = _getSequence('chr' + str(chr),int(tss),int(tss)+int(utr),logger)

            chunks, chunk_size = math.ceil(float(len(sequence))/FASTA_SEQ_LINE_LIMIT), FASTA_SEQ_LINE_LIMIT
            sequence_str = '\n'.join([ sequence[i:i+chunk_size] for i in range(0, int(chunks), chunk_size) ])

            out.write(">{gene_name} [organism={org}] {chr},{strand},{tss},{pileup},{utr},{transcript}\n{seq}\n".format(
                gene_name=gene,org="Mus musculus",
                chr=chr,strand=strand,tss=tss,pileup=pileup,
                utr=utr,transcript=transcript,
                seq=sequence_str.upper()))



def get_sequence(args):
    docstring = __doc__.format(log_level_vals=LOG_LEVEL_VALS)
    options = docopt.docopt(docstring,
                            version="get_empirical_utrs v" + __version__)
    validate_command_line_options(options)

    logger = log.get_logger(sys.stderr, options[LOG_LEVEL])

    _write_result(options[INPUT_FILE],options[OUTPUT_FILE],logger)