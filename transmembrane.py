seqs = list("MNSTPDLISPQKSNSSNSYELESGRSKAMNTPEGKNESFHDNLSESQVQPAVAPPNTGKGVYVTVSICCVMVAFGGFIFGWDTGTISGFVAQTDFLRRFGMKHHDGSHYLSKVRTGLIVSIFNIGCAIGGIVLAKLGDMYGRRIGLIVVVVIYTIGIIIQIASINKWYQYFIGRIISGLGVGGITVLSPMLISEVAPSEMRGTLVSCYQVMITLGIFLGYCTNFGTKNYSNSVQWRVPLGLCFAWALFMIGGMMFVPESPRYLVEAGRIDEARASLAKVNKCPPDHPYIQYELETIEASVEEMRAAGTASWGELFTGKPAMFQRTMMGIMIQSLQQLTGDNYFFYYGTIVFQAVGLSDSFETSIVFGVVNFFSTCCSLYTVDRFGRRNCLMWGAVGMVCCYVVYASVGVTRLWPNGQDQPSSKGAGNCMIVFACFYIFCFATTWAPIAYVVISECFPLRVKSKCMSIASAANWIWGFLISFFTPFITGAINFYYGYVFMGCMVFAYFYVFFFVPETKGLSLEEVNDMYAEGVLPWKSASWVPVSKRGADYNADDLMHDDQPFYKSLFSRK")

seqs=seqs[:230]
sosui = [[59,81],[114,136],[144,166],[174,196],[203,225]]
tmpred = [[37,63],[74,99],[115,140],[153,175],[203,221],[252,275],[285,308]]
tmap= [[43,69],[73,97],[112,140],[148,176],[201,229],[254,274],[281,301]]
tmhmm = [[39, 61], [74,96],[111,133],[153,175],[202,224],[253,275],[285,307]]
pred_tmr2 = [[37,54],[84,105],[114,133],[153,175],[203,221],[253,275]]
split = [[36,63],[71,99],[113,137],[153,173],[202,227],[250,275],[284,307]]
phobius = [[40,63],[75,95],[115,133],[153,173],[203,228],[249,275],[287,308]]
phdhtm = [[38,61],[76,96],[113,132],[154,171],[204,225],[254,271],[286,306]]
from fpdf import FPDF
list_of_tests = []
# list_of_tests.append(uniprot)
list_of_tests.append(sosui)
list_of_tests.append(tmpred)
list_of_tests.append(tmap)
list_of_tests.append(tmhmm)
list_of_tests.append(pred_tmr2)
list_of_tests.append(split)
list_of_tests.append(phobius)
list_of_tests.append(phdhtm)

tools = ["SOSUI","TMpred", "Tmap", "TMhmm","Pred-TMR2", "SPLIT", "Phobius", "PHDhtm"]

list_of_colors = [[205,92,92], [100,149,237], [50,205,50],[244, 208, 63], [218, 247, 166],[187, 143, 206],[131, 249, 84], [84, 249, 244]]
pdf = FPDF('L')
pdf.add_page()


for x in range(len(seqs)):
	# pdf.set_font("Arial", style = 'B', size = 0.9)
	# pdf.text(x/1.5, 3, str(x+1))
	pdf.set_font("Arial", style = 'B', size = 2)
	pdf.text(x/1.5, 4, seqs[x])


for i in range(len(list_of_tests)):
	for ele in list_of_tests[i]:
		pdf.set_fill_color(list_of_colors[i][0],list_of_colors[i][1],list_of_colors[i][2])
		pdf.rect((ele[0]-1)/1.5, (i*2) +5,(ele[1]-ele[0])/1.5, 1.5, style='F')

pdf.set_font("Arial", style = 'B', size = 14)
for g in range(len(list_of_colors)):
	pdf.set_fill_color(list_of_colors[g][0],list_of_colors[g][1],list_of_colors[g][2])
	pdf.rect(200, (g*11) +100,10, 10, style='F')
	pdf.text(212, (g*11)+106.3, tools[g])
pdf.output('transm.pdf')
