class Gene(object):
    def __init__(self, gtf_row):
        self.name = gtf_row.get_gene()
        self.seqname = gtf_row.get_seqname()
        self.transcripts = {}

    def add_transcript(self, gtf_row):
        transcript_name = gtf_row.get_transcript()

        transcript = None
        if transcript_name in self.transcripts:
            transcript = self.transcripts[transcript_name]
        else:
            transcript = Transcript(gtf_row)
            self.transcripts[transcript_name] = transcript

        return transcript


class Transcript(object):
    def __init__(self, gtf_row):
        self.name = gtf_row.get_transcript()
        self.strand = gtf_row.get_strand()
        self.exons = []

    def add_exon(self, gtf_row):
        self.exons.append(Exon(gtf_row))


class Exon(object):
    def __init__(self, gtf_row):
        self.start = gtf_row.get_start()
        self.end = gtf_row.get_end()

    def length(self):
        return self.end - self.start + 1
