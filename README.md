# Antibody_MD
These scripts are for running molecular dynamics simulation of antibody as well as antibody-antigen complex. We suggest to run in Linux.

To run MD simulation using MD.pl, please install the following software:
CUDA (compatable with Amber installation)
Amber (must be callable in commandline without absolute path, pmemd.cuda is required)
FoldX (for generating mutations by MD.pl, optional)
xmgrace (for plotting RMSD during MD simulation, optional)
Torque (for submitting jobs to queue, required for running mutliple MD runs at the same time)
MIPCH (for running pmemd.cuda.MPI, optional)

run perl MD.pl to see parameters required.

If you prefer to generate Amber topology files using other methods, please name the required topology and coordinate files to Temp.prmtop and Temp.inpcrd respectively.
