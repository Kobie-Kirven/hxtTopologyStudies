# hxtTopologyStudies
<h3>Predictions of transmembrane regions for hxt1-4 in Yeast,Saccharomyces cerevisiae</h3>
<p> Several different tools were used for the prediction of transmembrane regions of hxt1-4. The results from each prediction tool were aggregated for each
protein into a '.tmsr' file, which stands for transmembrane spanning regions. This file type resembles the FASTA filetype, with a few exceptions.
The '@' symbol is placed before each prediction tool used. The transmembrane regions are then placed on the line(s) following. The transmembrane
regions are indicated by the residue where the region begins, followed by a '-', followed by the residue where the region ends. Transmembrane 
regions are seperated by commas. An example line is show below</p>

@SOSUI \
59-81,114-136,144-166,174-196,203-225,236-258 \
360-382,390-412,429-451,470-492,496-518

After the .tmsr files were genrated, they were then used as input in tmComparisonPlot.py. This script generates a plot of the predicted 
transmembrane regions for each tool for easy comparision. 

