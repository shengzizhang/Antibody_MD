#!/usr/bin/perl
use strict;

if(@ARGV<1||@ARGV%2>0){die "
	Usage: 
	-f interface ASA file 
	-o output file 	
 
	 ";}

my %para=@ARGV;
	 
my %bonds=();
if(!$para{'-o'}){$para{'-o'}="interface_Hbond_trajectory.txt";}
open HH,"$para{'-f'}" or die "No input file\n";
my $i=0;

while(<HH>){
  chomp;
  $i++;
  	my @line=split/\t/,$_;
  	if($line[6] ne 'NA'){
  	  my @l=split/\;/,$line[6];
  	  #pop @l;
  	  
  	  foreach(@l){
  	  	$bonds{'H'}{$_}{$line[0]}=1;
  	  }	
  	}
  	if($line[7] ne 'NA'){
  	  my @l=split/\;/,$line[7];
  	  foreach(@l){
  	  	$bonds{'SB'}{$_}{$line[0]}=1;
  	  }	
  	}
}

open YU,">$para{'-o'}";
print YU "Bondtype\tresidue_pair\ttrajectory\tfound\n";
foreach('H','SB'){
  my $b=$_;
  foreach(sort keys %{$bonds{$b}}){
  	my $HB=$_;
  	foreach(1..$i){
  	  if($bonds{$b}{$HB}{$_}){
  	  	print YU "$b\t$HB\t$_\t1\n";
  	  }	
  	  else{
  	  	print YU "$b\t$HB\t$_\t0\n";
  	  }
  	}
  }	
}