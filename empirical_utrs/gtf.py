import pandas

from . import feature

SOURCE = "transcript_utils"
MISSING_VALUE = "."


class GtfRow(object):
    SEQNAME_COL = 0
    FEATURE_COL = 2
    START_COL = 3
    END_COL = 4
    STRAND_COL = 6
    ATTRIBUTES_COL = 8

    EXON_FEATURE = "exon"

    GENE_ID_ATTRIBUTE = "gene_id"
    TRANSCRIPT_ID_ATTRIBUTE = "transcript_id"

    @classmethod
    def from_file(cls, row_data):
        strip_quotes = lambda x: x.replace('"', '')
        attr_str = row_data[GtfRow.ATTRIBUTES_COL]
        attr_dict = {attr: strip_quotes(val) for attr, val in
                     [av.split(" ", 1) for av in attr_str.split("; ")]}

        return GtfRow(row_data, attr_dict)

    @classmethod
    def from_values(cls, seqname, feature_type, start, end,
                    strand, gene, transcript):

        row_data = [seqname, SOURCE, feature_type, start, end,
                    MISSING_VALUE, strand, MISSING_VALUE]
        attr_dict = {GtfRow.GENE_ID_ATTRIBUTE: gene,
                     GtfRow.TRANSCRIPT_ID_ATTRIBUTE: transcript}

        return GtfRow(row_data, attr_dict)

    def __init__(self, row_data, attr_dict):
        self.row_data = row_data
        self.attr_dict = attr_dict

    def get_seqname(self):
        return self.row_data[GtfRow.SEQNAME_COL]

    def get_feature(self):
        return self.row_data[GtfRow.FEATURE_COL]

    def get_start(self):
        return self.row_data[GtfRow.START_COL]

    def get_end(self):
        return self.row_data[GtfRow.END_COL]

    def get_strand(self):
        return self.row_data[GtfRow.STRAND_COL]

    def get_gene(self):
        return self.attr_dict[GtfRow.GENE_ID_ATTRIBUTE]

    def get_transcript(self):
        return self.attr_dict[GtfRow.TRANSCRIPT_ID_ATTRIBUTE]

    def is_exon(self):
        return self.get_feature() == GtfRow.EXON_FEATURE

    def __str__(self):
        fields = list(self.row_data)

        attr_str = "; ".join(["{k} \"{v}\"".format(k=k, v=v)
                            for k, v in self.attr_dict.iteritems()])
        fields.append(attr_str)

        return "\t".join([str(field) for field in fields])


class GtfInfo(object):
    def __init__(self, gtf_file, logger):
        self.gtf_file = gtf_file
        self.data = pandas.read_csv(
            gtf_file, sep="\t", header=None, comment="#")
        self.logger = logger

    def rows(self):
        for index, row in self.data.iterrows():
            yield GtfRow.from_file(row)

    def get_transcript_info(self):
        self.logger.info("Reading transcript info...")

        transcript_info = {}
        lines_processed = 0

        for row in self.rows():
            lines_processed += 1
            if lines_processed % 10000 == 0:
                self.logger.debug("Processed {l} GTF lines.".format(l=lines_processed))

            if not row.is_exon():
                continue

            gene_name = row.get_gene()

            gene = None
            if gene_name in transcript_info:
                gene = transcript_info[gene_name]
            else:
                gene = feature.Gene(row)
                transcript_info[gene_name] = gene

            transcript = gene.add_transcript(row)
            transcript.add_exon(row)

        self.logger.info("...read transcript information for {g} genes".format(
            g=len(transcript_info)))

        return transcript_info
