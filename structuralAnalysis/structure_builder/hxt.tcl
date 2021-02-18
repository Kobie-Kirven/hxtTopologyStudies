#Preparation of PDB and PSF files for NAMD molecular
#dynamics simulations of yeast glucose transporters

#The input pdb file name and output name
# set filename full_model_54606_1.pdb
set outPsf out/hxt1.psf
set outPdb out/hxt1.pdb

#Indicate the anchoring residues
set nEnd 61
set cStart 512

#set the distance between the alpha carbons of the 
# anchoring residues (Angstroms) 
set Distance 10

#Load the psfgen package and specify the residue /
# atom conversions
package require psfgen
topology top_all36_prot.rtf
topology 
pdbalias residue HIS HSE
pdbalias atom ILE CD1 CD

#Create a psf file and pdb file with H's for the N termius
segment N {pdb out/hxt1_N_terminus.pdb}
coordpdb out/hxt1_N_terminus.pdb N
guesscoord
writepdb out/nterm_h.pdb
writepsf out/nterm.psf

#Reset the psf file
resetpsf 

#Create a psf file and a pdb file with H's for the C terminus
segment C {pdb out/hxt1_C_terminus.pdb}
coordpdb out/hxt1_C_terminus.pdb C
guesscoord
writepdb out/cterm_h.pdb
writepsf out/cterm.psf

#Combine the individual pdb and psf files into one
package require topotools
set midlist {}
set mol [mol new out/nterm.psf waitfor all]
mol addfile out/nterm_h.pdb
lappend midlist $mol
set mol [mol new out/cterm.psf waitfor all]
mol addfile out/cterm_h.pdb $mol
lappend midlist $mol

set mol [::TopoTools::mergemols $midlist]
animate write psf $outPsf $mol
animate write pdb $outPdb $mol

#Delete the molecule
mol delete top

#Open the file with both termani
mol new $outPdb

##################################################
# Fix Specified Atoms
##################################################

mol new out/hxt1.pdb waitfor all
#Set the beta value for all of the atoms to zero
set sel [atomselect top all]
$sel set beta 0



#Set the beta value to 1 for the alpha carbons of each
# of the anchoring residues
set nFix [atomselect top "resid $nEnd and name CA"]
$nFix set beta 1
set cFix [atomselect top "resid $cStart and name CA"]
$cFix set beta 1


set c_term [atomselect top "resid $cStart to 570"]
#Get the x coordinate for the alpha carbon in the fixed
# residue of the N and C-termani
set nX [$nFix get {x}]
set cX [$cFix get {x}]

# #Calculate the distance the x value needs to be moved to
# # satisfy the value set by the user
set x_dist [expr $nX - $cX]

# #Get the Y coordinate for the alpha carbon in the fixed
# # residue of the N and C-termani
set nY [$nFix get {y}]
set cY [$cFix get {y}]

# #Calculate the distance the Y value needs to be moved so
# # that the Y values are the same between the alpha carbons
# # of the fixed residues
set y_dist [expr $nY - $cY]

# #Get the Z coordinate for the alpha carbon in the fixed
# # residue of the N and C-termani
set nZ [$nFix get {z}]
set cZ [$cFix get {z}]

# #Calculate the distance the Z value needs to be moved so
# # that the Z values are the same between the alpha carbons
# # of the fixed residues
set z_dist [expr ($nZ - $cZ) - $Distance]


# #Add each of the coordinates to the coordinates vector
lappend vect $x_dist
lappend vect $y_dist
lappend vect $z_dist

# #Move the C-terminus by the specified values
$c_term moveby $vect

$sel writepdb $outPdb

quit

