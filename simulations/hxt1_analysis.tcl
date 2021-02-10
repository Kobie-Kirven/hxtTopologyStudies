#Script for analyzing the RMSD from namd output
#Kobie Kirven - 1/30/2021

set psfFile hxt1.psf
set trajFile hxt1_simulations.dcd

#Load the protein structure file
mol new $psfFile type psf

#Load the trajectory file
mol addfile $trajFile type dcd first 0 last -1 step 1 waitfor all

#Set a reference frame for the rmsd calculation
set reference [atomselect top "protein and backbone" frame 0]

#Set the frame to compare to
set compare [atomselect top "protein and backbone"]

#Get the number of frames in the trajectory file
set num_steps [molinfo top get numframes]

#Set output file
set outfile [open hxt1_rmsd.dat w]

for {set frame 0} {$frame < $num_steps} {incr frame} {
			#get the correct frame
			$compare frame $frame

			set trans_mat [measure fit $compare $reference]

			$compare move $trans_mat

			set rmsd [measure rmsd $compare $reference]
	puts $outfile "$frame 	$rmsd"
}

close $outfile


#RMSF calculation
set num [expr {$num_steps -1}]

set outfile [open rmsf.dat w]
set sel [atomselect top "protein and name CA"]
set rmsf [measure rmsf $sel first 0 last $num step 1]
for {set i 0} {$i < [$sel num]} {incr i} {
	puts $outfile "[expr {$i+1}] [lindex $rmsf $i]"
}
close $outfile


#hydrogen_bond analysis
package require hbonds
hbonds -sel1 [atomselect top protein] -writefile yes -plot no

quit
