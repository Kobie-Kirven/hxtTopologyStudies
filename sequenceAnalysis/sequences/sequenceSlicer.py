from Bio import SeqIO

sequences = []
for record in SeqIO.parse("yeast_transporters_protein.fasta", "fasta"):
	sequences.append(list(record.seq))

#Nterm = ''.join(sequences[0][:59])


indList = [524
,503
,518
,518
,519
,519
,520]

for i in indList:
        print(sequences[0][i+1])
