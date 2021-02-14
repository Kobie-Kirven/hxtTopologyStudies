#Imports 
from Bio import PDB, SeqIO
from PeptideBuilder import Geometry
import PeptideBuilder

def rev(x):
	reverse = ''
	for t in range(len(x)):
		reverse += x[-(t+1)]
	return reverse

#N-terminus end C-terminus start
nEnd = 62
cStart = 512
cEnd = 570

#Transporter
Transporter = "hxt1"

sequences = []

for rec in SeqIO.parse("/Users/kobiekirven/Documents/GitHub/yeastGlucoseTransporters/sequenceAnalysis/sequences/yeast_transporters_protein.fasta", 'fasta'):
	sequences.append(str(rec.seq))


def buildTerminus(sequence, start, stop, transporter, NC):
	if NC=='C':
		term = PeptideBuilder.make_extended_structure(rev(sequence[start:stop]))
		out = PDB.PDBIO()
		out.set_structure(term)
	else:
		term = PeptideBuilder.make_extended_structure(sequence[start:stop])
		out = PDB.PDBIO()
		out.set_structure(term)
	if NC == 'N':
		out.save("out/" + transporter + "_N_terminus.pdb")
	else:
		out.save("out/" +transporter + "_C_terminus.pdb")


#Build the N-terminus
buildTerminus(sequences[0], 0, (nEnd-1), "hxt1", NC="N")

#Build the C-terminus
buildTerminus(sequences[0], (cStart-1), (cEnd), "hxt1", NC="C")

pdb_io = PDB.PDBIO()
pdb_parser = PDB.PDBParser()
pdbfile = 'out/hxt1_C_terminus.pdb'
structure = pdb_parser.get_structure(" ", pdbfile)

new_resnums = [i + (cStart) for i in range((cEnd - cStart)+1)]

for model in structure:
    for chain in model:
        for i, residue in enumerate(chain.get_residues()):
            res_id = list(residue.id)
            res_id[1] = new_resnums[-(i+1)]
            residue.id = tuple(res_id)

pdb_io.set_structure(structure)
pdb_io.save(pdbfile)


# fn = open(pdbfile)
# lines = fn.readlines()
# fn.close()

