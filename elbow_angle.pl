#!/usr/bin/perl
use strict;

# This script is used to calculate sequence identity between germline V, antibody gene and reads
if(@ARGV%2>0||@ARGV==0){die "Usage: identity.pl 
	-f pdbfile 
	-h heavy chain 
	-l light chain
	-hend last residue of heavy chain variable domain 
	-lend last residue of light chain variable domain  
	\n";}
my %para=@ARGV;
if(!$para{'-f'}){die'no input file\n';}
if(!$para{'-h'}){$para{'-h'}='H';}
if(!$para{'-l'}){$para{'-l'}='L';}
if(!$para{'-hend'}){$para{'-hend'}=113;}
if(!$para{'-lend'}){$para{'-lend'}=107;}
if(!$para{'-pymol'}){$para{'-pymol'}='pymol';}
my $pdbname=$para{'-f'};
$pdbname=~s/.pdb$//;
if($pdbname=~/\//){my @l=split/\//,$pdbname;$pdbname=$l[$#l];}
###############Processing###################
print &pymol(),"\n";

###########################################
sub pymol{
    open HH,">pymol.pml";	
print HH "import elbow_angle
load $para{'-f'}
elbow_angle $pdbname, heavy=\'$para{'-h'}\',light=\'$para{'-l'}\',limit_h=$para{'-hend'},limit_l=$para{'-lend'} 
quit
";
my @l=`$para{'-pymol'} -c pymol.pml`;
my $angle=0;
foreach(@l){
	if($_=~/Elbow angle: (\d+) degrees/){
	   $angle=$1;	
	}	
}
system("rm pymol.pml");
if($angle==0){die "something wrong with calculaiton: $para{'-f'}\n";}
else{return $angle;}
}
