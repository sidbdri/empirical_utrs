The script ``get_empirical_utrs`` calculates empirical 5' UTRs from CAGE data, using a procedure adapted from Khajuria et al., _Ribosome Levels Selectively Regulate Translation and Lineage Commitment in Human Hematopoiesis_, Cell (2018):

i) For each gene in the supplied GTF file, determine the region to be scanned for the empirical transcription start site, namely from 1000 bases upstream of the earliest location of the start of any transcript, to the most downstream location of the end of any transcript.

ii) Within these scan bounds find the location with greatest pile up of reads mapped to the same strand as the gene. In the case of ties, choose the rightmost location for genes on the plus strand, and leftmost location for genes on the minus strand. This location is defined to be the empirical TSS for the gene.

iii) Choose the transcript for which the empirical TSS is the smallest number of bases upstream (> 0) of the defined coding start location. If the TSS is downstream of the coding start location, ignore this transcript. If the transcript has no coding start location defined, use the location of the start of the first exon. The empirical 5' UTR is defined as the region between the empirical TSS and the coding start location of the chosen transcript.

iv) Reject this empirical 5' UTR if it is longer than 500 bases.

The script ``get_utr_sequences`` can be used to obtain a FASTQ file of UTR sequences from the output of ``get_empirical_utrs``.
