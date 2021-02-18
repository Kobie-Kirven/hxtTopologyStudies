from Bio import SeqIO

sequences = []
for rec in SeqIO.parse("SNF3_protein.fsa", "fasta"):
	sequences.append(str(rec.seq))

for rec in SeqIO.parse("RGT2_protein.fsa", "fasta"):
	sequences.append(str(rec.seq))

motifs = []
for i in range(len(sequences[0])):
	for z in range(len(sequences[1])):
		if sequences[0][i] == sequences[1][z]:
			flag = True
			index = i
			mini = [i]
			while flag == True and index < (len(sequences[1])-1):
				if sequences[0][index + 1]==sequences[1][index + 1]:
					index +=1
				else:
					mini.append(index)
					if (mini[1] - mini[0]):
						motifs.append(mini)
					flag = False


print(motifs)
