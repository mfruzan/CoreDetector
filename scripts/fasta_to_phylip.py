#!/usr/bin/env python3
from Bio import SeqIO

records = SeqIO.parse("concatinated_msa.fa", "fasta")
count = SeqIO.write(records, "concatinated_msa.phylip", "phylip")
print("Converted %i records" % count)
