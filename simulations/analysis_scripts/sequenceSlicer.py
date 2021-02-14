from Bio import SeqIO

for record in SeIO.parse("yeast_transporters_protein.fasta"):
	print(record.seq)