from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio import SeqIO
from Bio.SeqUtils import ProtParamData
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import string

sequences = []

matplotlib.rc('font',family='Arial')

fig, axs = plt.subplots(4, sharex=True, sharey=True)


titles = ['Hxt1p','Hxt2p','Hxt3p','Hxt4p']

axs[0].set_xticks(np.arange(0, 570, 50))


fn = open('hydrophobicOutput.txt', 'w')
for rec in SeqIO.parse("../sequenceAnalysis/sequences/yeast_transporters_protein.fasta", 'fasta'):
	sequences.append(str(rec.seq))

for i in range(len(sequences)):
	X = ProteinAnalysis(sequences[i])
	kd = X.protein_scale(ProtParamData.kd, 19, edge=0.1)
	num_1 = [g + 1 for g in range(len(kd))]
	axs[i].plot(num_1, kd)
	axs[i].plot(num_1, [1.6 for g in range(len(kd))], color='red', linewidth=0.4)
	axs[i].grid(True)
	axs[i].margins(tight=True, x=-0, y=None)
	
	for k in kd:
		if k >= 1.6:
			fn.write(titles[i] + '\t'  + sequences[i][kd.index(k)] 
			 + '\t'+ str(kd.index(k)+1) + '\t'+  str(k) + '\n')

b = 0
for ax in axs:
    ax.text(8, 2, string.ascii_uppercase[b], 
                            size=10, weight='bold', fontname="Arial") 
    b += 1 
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', top=False, bottom=False, left=False, right=False)
plt.ylabel("Hydropathy Index", labelpad=10, fontweight='bold')
plt.xlabel("Residue Position", labelpad=10, fontweight='bold')
plt.subplots_adjust(hspace=0.1)
plt.title("Hydrophobicity of Hxt1-4p", fontweight='bold',fontsize=14)

fn.close()

plt.tight_layout()

plt.savefig('Hydrophobicity.png', dpi=800)