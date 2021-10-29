#!/usr/bin/perl
use strict;
use Cwd qw(abs_path);
if(!@ARGV||@ARGV%2>0){
die"
Usage:
-pdb  pdb file
-mut mutations, e.g. QH39L,YH95F
-t number of threads
-aaid amino acid positions in the structure:1-427
-GaMD whether to perform accelerated MD, 1 or 0 
-nsteps default 200ns 
-c reorder pdb positions and generate mutations, 1 for yes, 2 for no, default :1
-MPI using multiple GPU
-tleap whether to use tleap to generate topology files, 1 for yes, 2 for no. Default 1
-vamber, version of Amber software, default: 18
-foldx path to the FoldX installation folder 
";
}

my %para=@ARGV;
if(!$para{'-dev'}){$para{'-dev'}='0';}
if(!$para{'-t'}){$para{'-t'}=1;}
if(!$para{'-tleap'}){$para{'-tleap'}=1;}
if(!$para{'-c'}){$para{'-c'}=1;}
if(!$para{'-MPI'}){$para{'-MPI'}=0;}
if(!$para{'-nsteps'}){$para{'-nsteps'}=100000000;}
if(!$para{'-vamber'}){$para{'-vamber'}=18;}
if($para{'-mut'} && !$para{'-foldx'}){die "FoldX is needed for mutating structure\n";}
if(!-d $para{'-foldx'}){die "FoldX folder not exist\n";}
my $nsteps=$para{'-nsteps'};
my $store=$nsteps/1000;
my $aan=0;
my $path=abs_path($0);


if($para{'-c'}==1){
  &clean($para{'-pdb'},$para{'-mut'});#clean pdb file and generate mutations
  if(-e 'Temp.pdb'){
     $aan=&aanum('Temp.pdb');
     print "total amino acids $aan\n";
   }
  else{die "Temp.pdb not found\n";}

}
else{
  $aan=&aanum($para{'-pdb'});
}

if($para{'-tleap'}==1){
	&tleap('Temp.pdb');
}
&min1();

&min2();

&heating();
&printrmsd('heating.mdcrd','heating');

&equil(2000000);
&printrmsd('equil.mdcrd','equil');
&production();
&printrmsd('production.mdcrd','production');

&cppj('production');

&hbond('production');
&convert_image('production.mdcrd','production');
system("mv production_reimage.nc production.mdcrd");


############sub routines################

sub cppj{
    my ($name)=@_;
    my $eq='equil';
    open YY,">rmsd.cpptraj";
    print YY "trajin heating.mdcrd
    trajin $eq.mdcrd
    trajin $name.mdcrd
    reference Temp_MD.pdb
    autoimage
    rms reference :1-$aan\@CA mass out heating_production.rms

";system("cpptraj -p Temp.prmtop -i rmsd.cpptraj &> rmsd.log ");
}

sub hbond{#calculate hydrogen bonds in snapshots
    my ($name)=@_;
    open YY,">hbond.cpptraj";
    print YY "
    trajin $name.mdcrd
    reference min2.rst
    autoimage
    hbond out hbond_total_per_frame.txt avgout hbond_each_per_frame.txt series uuseries hbond_solute_timeseries.txt uvseries hbond_water.txt

";system("cpptraj -p Temp.prmtop -i hbond.cpptraj &> hbond.log ");
&Hbond_pos();
}
sub production{
    open HH,">production.in";
    print HH "
   Production
 \&cntrl
  imin=0,
  ntx=5, !Read coordinates and velocities from unformatted inpcrd file
  irest=1,!Restart previous MD run [This means velocities are expected in the inpcrd file and will be used to provide initial atom velocities]
  nstlim=$nsteps, !steps
  dt=0.002,!time step
  ntf=2,
  ntc=2,
  temp0=300.0,
  ntpr=$store,
  ntwx=$store,
  ntwr=$store,
  ntb=2, !use periodic boundary conditions with constant pressure
  ntp=1, !Use the Berendsen barostat for constant pressure simulation
  ntt=3,
  iwrap=1,
  gamma_ln=1.0,
  cut=8,
  taup=5,
  rgbmax=15,
 /
"; 
 #system("export CUDA_VISIBLE_DEVICES=$para{'-dev'};mpirun -np $para{'-t'} pmemd.cuda.MPI -O -i production.in -o production.out -c equil.rst -p Temp.prmtop -r production.rst -x production.mdcrd -inf production.mdinfo -suffix production ");
 if($para{'-MPI'}){system("mpirun -np $para{'-t'} pmemd.cuda.MPI -O -i production.in -o production.out -c equil.rst -p Temp.prmtop -r production.rst -x production.mdcrd -inf production.mdinfo -suffix production ");}
else{system("pmemd.cuda -O -i production.in -o production.out -c equil.rst -p Temp.prmtop -r production.rst -x production.mdcrd -inf production.mdinfo -suffix production ");}
print "production run finished\n";
}


sub heating{
open HH,">heating.in";
print HH "
Heating up the system equilibration stage 1
 \&cntrl
  nstlim=2500000, dt=0.002, ntx=1, irest=0, ntpr=20000, ntwr=20000, ntwx=20000,
  tempi =0.0, temp0=300.0, ntt=3, gamma_ln=1.0,taup=5,
 
  ntb=1, ntp=0,
  cut=8,
  ntc=2, ntf=2,
  nmropt=1, 
  iwrap=1,
  nrespa=1,!frequency of steps to evaluate slowly-varying terms
/

&wt type='TEMP0', istep1=0, istep2=900000, value1=0.0, value2=300.0 /
&wt type='TEMP0', istep1=900001, istep2=2500000, value1=300.0, value2=300.0 /
&wt type='END' /

";
close HH;
    if($para{'-MPI'}){system("mpirun -np $para{'-t'} pmemd.cuda.MPI -O -i heating.in -o heating.out -c min2.rst -x heating.mdcrd -p Temp.prmtop -r heating.rst -inf heating.mdinfo -suffix heating ");}
else{system("pmemd.cuda -O -i heating.in -o heating.out -c min2.rst -x heating.mdcrd -p Temp.prmtop -r heating.rst -inf heating.mdinfo -suffix heating ");}
print "heating finished\n";
}

sub equil{
  my ($steps)=@_;
open HH,">equil.in";
print HH "
Equilibration stage 2
 \&cntrl
  imin=0,
  ntx=5, !Read coordinates and velocities from unformatted inpcrd file
  irest=1,!Restart previous MD run [This means velocities are expected in the inpcrd file and will be used to provide initial atom velocities]
  nstlim=$steps, !steps
  dt=0.002,!time step
  ntf=2,
  ntc=2,
  temp0=300.0,
  ntpr=20000,
  ntwx=20000,
  ntwr=20000,
  ntb=2, !use periodic boundary conditions with constant pressure
  ntp=1, !Use the Berendsen barostat for constant pressure simulation
  ntt=3,
  cut=8,
  !fswitch=9,
  gamma_ln=1.0,
  iwrap=1,
  !nsnb= 20, !frequency of nonbonded list updates
  !nbflag =1, !pairlist of atoms update whenever any atome moves 1/2 skinnb angstrom
  !skinnb =2, 
  taup=5,
  rgbmax=15,
 /
";
close HH;
    if($para{'-MPI'}){system("mpirun -np $para{'-t'} pmemd.cuda.MPI -O -i equil.in -o equil.out -c heating.rst -x equil.mdcrd -p Temp.prmtop -r equil.rst -suffix equil -inf equil.mdinfo");}
else{system("pmemd.cuda -O -i equil.in -o equil.out -c heating.rst -x equil.mdcrd -p Temp.prmtop -r equil.rst -suffix equil -inf equil.mdinfo");}
print "Equilibration finished\n";
}

sub min1{
    open HH,">min1.in";
    print HH "
    Minimize hydrogens
 \$cntrl
  imin   = 1,! (0 for simulation without minimization)
  maxcyc = 10000,
  ncyc   = 5000, !(after steps of steepest descent, switch to conjugate gradient)
  ntmin  = 1, !(0, full conjugate;1, steepest descent and conjugate;2,steepest descent;3, XMIN;4,LMOD)
  ntpr=10, !(print energy frequency)
  ntr=1, !(turn on Cartesian restraints)
  restraint_wt=9999.0, !(force constant for restraint)
  restraintmask=\'\:1-$aan\', !(atoms in residues 1-58 restrained)
  ntwx   = 100, !write crd file every ntwx steps
  ntwe   = 100, !write energy
  ntwf   = 100, !write force
  ntwv   = 100, !write velocity
  ntwprt = 0, ! how many atoms to write to output crd and vel file,0 is for all atoms
  idecomp = 0, !energy decomposition
/
    ";
    close HH;
    system("pmemd.cuda -O -i min1.in -o min1.out -c Temp.inpcrd -p Temp.prmtop -x min1.mdcrd -r min1.rst -inf min1.mdinfo -ref Temp.inpcrd -suffix min1 ");
print "minimization solvent finished\n";
}

sub min2{
    open HH,">min2.in";
    print HH "
    Minimize fully
 \$cntrl
  imin   = 1,! (0 for simulation without minimization)
  maxcyc = 10000,
  ntpr   = 100,
  ncyc   = 5000, !(after steps of steepest descent, switch to conjugate gradient)
  ntmin  = 1, !(0, full conjugate;1, steepest descent and conjugate;2,steepest descent;3, XMIN;4,LMOD)
  ntwx   = 100, !write crd file every ntwx steps
  ntwe   = 100, !write energy
  ntwf   = 100, !write force
  ntwv   = 100, !write velocity
  ntwprt = 0, ! how many atoms to write to output crd and vel file,0 is for all atoms
  idecomp = 0, !energy decomposition
/
    ";
    close HH;
    system("pmemd.cuda -O -i min2.in -o min2.out -c min1.rst -p Temp.prmtop -x min2.mdcrd -inf min2.mdinfo -r min2.rst -suffix min2 ");
print "minimization whole system finished\n";
}

sub clean{
    my($pdb,$mut)=@_;
	if($mut){
		open M,">individual_list.txt";
		print M "$mut\;\n";
		system("cp $para{'-foldx'}/rotabase.txt .");
		system("$para{'-foldx'}/foldx --command=BuildModel --pdb=$pdb --mutant-file=individual_list.txt");
		my $newpdb=$pdb;
		$newpdb=~s/.pdb/\_1.pdb/;
		if(-e $newpdb){
			open NP,"$newpdb";
		open PO,">mutant_reformat.pdb";
			<NP>;
			<NP>;
			<NP>;
		while(<NP>){
 		  print PO "$_";
		 }
		 print PO "\nTER\n";
		close PO;close NP;
   		    system("pdb4amber -i mutant_reformat.pdb --reduce --no-reduce-db -d -y -l pdb4amber.log -o Temp.pdb");
		}
		else{die "foldx run wrong\n";}
		unlink "rotabase.txt";
	}
    	else{
		`pdb4amber -i $pdb --reduce --no-reduce-db -l pdb4amber.log -d -y -o Temp.pdb`;
		
	}
}

sub tleap{
    my ($pdb)=@_;
open HH,">tleap.in";
print HH '
source leaprc.protein.ff14SB #Source leaprc file for ff14SB protein force field
source leaprc.GLYCAM_06j-1
source leaprc.DNA.OL15

source leaprc.lipid17
source leaprc.water.tip3p
source leaprc.gaff2
 
mol = loadpdb ',$pdb,'  #Load PDB file for protein-ligand complex

';
close HH;
if($para{'-b'})
 {
       system("cat $para{'-b'} >> tleap.in" );
 }
open HH,">>tleap.in";
print HH '
solvatebox mol TIP3PBOX 10 #Solvate the complex with a cubic water box
addions2 mol Cl- 0 #Add Cl- ions to neutralize the system
addions2 mol Na+ 0 #Add Na+ ions to neutralize the system

saveamberparm mol Temp.prmtop Temp.inpcrd #Save AMBER topology and coordinate files
quit #Quit tleap program
';
close HH;
system("tleap -f tleap.in ");
open YY,">parmed.in";
print YY "
parm Temp.prmtop
addPDB Temp.pdb 
outparm Temp_chainIDadded.prmtop
";
close YY;
system("parmed -i parmed.in");
#system("ambpdb -p Temp.prmtop -c Temp.inpcrd -ext -bres -conect >Temp_MD.pdb");
#system("ambpdb -p Temp_chainIDadded.prmtop -c Temp.inpcrd -ext -bres -conect >Temp_MD.pdb");
system("ambpdb -p Temp_chainIDadded.prmtop -c Temp.inpcrd -ext -conect >Temp_MD.pdb");
open ZZ,"Temp_MD.pdb";
open ZZZ,">Temp_MD1.pdb";
	while(<ZZ>){
	~s/CYX/CYS/;
 	 print ZZZ "$_";
	}
close ZZ;
close ZZZ;
system("mv Temp_MD1.pdb Temp_MD.pdb");
}

sub aanum{
    my ($pdb)=@_;
   open HH,"$pdb" or die "$pdb not found\n";
   open YY,">Temp1.pdb";
	my $aanum=0;
	my $Hatom='';
   while(<HH>){
   	if($_=~/ATOM/){
	   $aanum=substr($_,length('ATOM   3326  OXT GLU L'),4);
	   $aanum=~s/ +//g;
	   if(substr($_,length('ATOM   3309  NE  ARG '),1) eq 'H'){
		$Hatom=substr($_,length('ATOM  '),5);
		$Hatom=~s/ +//g;
	   }	
	  print YY $_;
	}
	elsif($_=~/CONECT/ && $para{'-vamber'}==16){
		chomp;
		my $atom1=substr($_,6,5);
		my $atom2=substr($_,11,5);
		$atom1=~s/ +//g;
		$atom2=~s/ +//g;
		if($atom1>$Hatom){
		   $atom1++;
		   $atom2++;
		}
		while(length($atom1)<5){
		  $atom1=" $atom1";
		}
		while(length($atom2)<5){
		  $atom2=" $atom2";
		}
         #$atom2=~s/3358/3359/;
	print YY "CONECT$atom1$atom2\n";
	}
	else{print YY $_;}

   }
close YY;
system("mv Temp1.pdb Temp.pdb");
   return $aanum;
}

sub printrmsd{
    my ($crdfile,$name)=@_;
    open HH,">$name\_cpptrajrmsd.in";
print HH "
trajin $crdfile
reference Temp_MD.pdb [min]
autoimage
rms min out $name.rms \@CA
run
";
close HH;
system("cpptraj -p Temp.prmtop -i $name\_cpptrajrmsd.in >> rmsd.log ");
system("xmgrace $name.rms -hardcopy -printfile $name\_rmsd.eps");
}
sub convert_image{
    my ($crdfile,$name)=@_;
    open HH,">$name\_cpptrajrmsd.in";
print HH "
trajin $crdfile
reference Temp_MD.pdb [min]
autoimage
trajout production_reimage.nc
run
";
system("cpptraj -p Temp.prmtop -i $name\_cpptrajrmsd.in >> rmsd.log ");
close HH;
}
sub Hbond_pos{
    open HH,"hbond_solute_timeseries.txt" or die "Hbond file not found\n";
    my $title=<HH>;
    my %Hbond=();
	my %pos_total=();
	chomp $title;
    my @names=split/[\t ]+/,$title;
    for(my $i=1;$i<@names;$i++){
       my @res=split/\-/,$names[$i];
        $res[0]=~/([A-Z]+)\_(\d+)/;
        $res[0]="$2$1";
	#my $pos1=$2; 
        $res[1]=~/([A-Z]+)\_(\d+)/;
        $res[1]="$2$1";	
	#my $pos2= 
      $names[$i]=join '_',sort {$a<=>$b} ($res[0],$res[1]); 
	$pos_total{$res[0]}=1;
	$pos_total{$res[1]}=1; 
     }
    while(<HH>){
	chomp;
      my @l=split/[\t ]+/,$_;
     	for(my $i=2;$i<@l;$i++){
  	    my @pos=split/\_/,$names[$i];
	    $Hbond{$l[1]}{$pos[0]}+=$l[$i];
	    $Hbond{$l[1]}{$pos[1]}+=$l[$i];
	}
    }
  open YY,">Hbond_per_pos_total.txt";
	foreach(sort {$a<=>$b} keys %Hbond){
		my $f=$_;
		foreach(sort {$a<=>$b} keys %pos_total){
			if($Hbond{$f}{$_}){
				print YY "$f\t$_\t$Hbond{$f}{$_}\n";
			}
			else{print YY "$f\t$_\t0\n";}
		}
	}

	close HH;
	close YY;

}
