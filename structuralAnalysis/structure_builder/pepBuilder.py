#Imports 
from Bio import PDB, SeqIO
from PeptideBuilder import Geometry
import PeptideBuilder


#N-terminus end C-terminus start
nEnd = 60
cStart = 512
cEnd = 570

#Transporter
Transporter = "hxt1"

sequences = []

for rec in SeqIO.parse("../sequenceAnalysis/sequences/yeast_transporters_protein.fasta", 'fasta'):
	sequences.append(str(rec.seq))


def buildTerminus(sequence, start, stop, transporter, NC):
	term = PeptideBuilder.make_extended_structure(sequence[start:stop])
	out = PDB.PDBIO()
	out.set_structure(term)
	if NC == 'N':
		out.save(transporter + "_N_terminus.pdb")
	else:
		out.save(transporter + "_C_terminus.pdb")


#Build the N-terminus
buildTerminus(sequences[0], 0, 60, "hxt1", NC="N")

#Build the C-terminus
buildTerminus(sequences[0], 512, 569, "hxt1", NC="C")

pdb_io = PDB.PDBIO()
pdb_parser = PDB.PDBParser()
pdbfile = 'hxt1_C_terminus.pdb'
structure = pdb_parser.get_structure(" ", pdbfile)

new_resnums = [i + cStart for i in range(cEnd - cStart)]

for model in structure:
    for chain in model:
        for i, residue in enumerate(chain.get_residues()):
            res_id = list(residue.id)
            res_id[1] = new_resnums[i]
            residue.id = tuple(res_id)

pdb_io.set_structure(structure)
pdb_io.save(pdbfile)