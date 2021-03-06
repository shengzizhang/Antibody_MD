# Antibody_MD
These scripts are for running molecular dynamics simulation of antibody. We suggest to run in Linux.

###MD.pl###

To run MD simulation using MD.pl, please install the following software:
CUDA (compatable with Amber installation)

Amber (must be callable in commandline without absolute path, pmemd.cuda is required)

FoldX (for generating mutations by MD.pl, optional)

xmgrace (for plotting RMSD during MD simulation, optional)

Torque (for submitting jobs to queue, required for running mutliple MD runs at the same time)

MIPCH (for running pmemd.cuda.MPI, optional)


run perl MD.pl to see parameters required.

If you prefer to generate Amber topology files using other methods, please name the required topology and coordinate files to Temp.prmtop and Temp.inpcrd respectively.

####Traj.R####

The Traj.R script used R to analyze antibody features (VH-VL angle, elbow angle, buried accessible surface area, PCA, etc.) from MD trajectory.

The script requires the following programs installed and can be called from command line.

R packages: bio3d,ncdf4,igraph,ggplot2,gplots,reshape2

PISA from ccp4 package: https://www.ccp4.ac.uk/ 

PyMOL with elbow_angle script:https://pymolwiki.org/index.php/Elbow_angle

ANARCI: http://opig.stats.ox.ac.uk/webapps/newsabdab/sabpred/anarci/#:~:text=ANARCI%20is%20a%20tool%20for,T%2Dcell%20receptor%20variable%20domains.&text=TCR%20sequences%20can%20only%20be,antibody%20and%20TCR%20domain%20types

TMalign: https://zhanggroup.org/TM-align/

To run the script, please first add the installation paths of the above programs to the Traj.R script.

  chmod +x Traj.R
  
  cd 'MD_trajectory_folder generated by MD.pl'
  
  Traj.R
  
If topology files and trajectory files were not generated by Amber, please provide a pdb file for MD simulation and a trajectory file in ncdf format, and name them Temp_MD.pdb (with H for heavy chain and L for light chain) and production.mdcrd respectively.
