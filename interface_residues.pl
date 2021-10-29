#!/usr/bin/perl
use strict;

if(@ARGV<1||@ARGV%2>0){die "
	Usage: 
	-f pdb file
	-c1 first chain of interaction
	-c2 second chain of interaction 
	-o output file 
	-ASA output buried surface area file
	-frame output id for each record in the ASA file
	 ";}

my %change=("ALA","A", "CYS","C", "CYX","C","ASP","D", "GLU","E", "PHE","F", "GLY","G", "HIS","H", "ILE","I","LYS","K", "LEU","L", "MET","M", "ASN","N", "PRO","P", "GLN","Q", "ARG","R", "SER","S", "THR","T","VAL", "V", "TRP","W", "TYR","Y","MSE","M","HSE","H","ASX","N"); 
my %para=@ARGV;
my $pdb_name=$para{'-f'};
$pdb_name=~s/\.pdb//;	 
my %bonds=();
my $randnm=int(rand(99999999));
my $pisaprojname=$pdb_name.$randnm;
if(!$para{'-o'}){$para{'-o'}=$pdb_name."_interface_residue.txt";}
if(!$para{'-ASA'}){$para{'-ASA'}="Buried_interface_area.txt";}
#system("cp $para{'-f'} change.pdb");
#copy($para{'-f'}, "change.pdb");
system("pisa $pisaprojname -analyse $para{'-f'} /home/sheng/work/software/ccp4/share/pisa/pisa.cfg");
system("pisa $pisaprojname -list interfaces /home/sheng/work/software/ccp4/share/pisa/pisa.cfg >temp_interface.txt");
open YY,">$para{'-o'}";
open ASA, ">>$para{'-ASA'}";
#print ASA "frame\tinterface\tH_num_res\tL_num_res\tH_bas\tL_bsa\tHbonds\tSaltb\n";
my ($interface_ids,$chains)=&interface('temp_interface.txt',0);
my $i=0;
foreach(@$interface_ids){
	my $interface_id=$_;
	
	my $chain_paire=$chains->[$i];$i++;
print "Interface ID is $interface_id\n";
system("pisa $pisaprojname -detail interface $interface_id /home/sheng/work/software/ccp4/share/pisa/pisa.cfg >temp_interface.txt");
my $residues=&residues('temp_interface.txt');
my %seq_pos=&pdb_pos($para{'-f'});
&interface('temp_interface.txt',1);


foreach(sort keys %{$residues}){
	 my $chain=$_;
	 
	 foreach(sort {$a<=>$b} keys %{$residues->{$chain}}){
	     my $pos=$_;
	     if($residues->{$chain}{$pos}=~/[0-9]/){
	     	$residues->{$chain}{$pos}='HDR';
	     }	
	     print YY "$pdb_name\t$chain\t$chain_paire\t$i\t$pos\t$residues->{$chain}{$pos}\n";#\t$seq_pos{$chain}{$pos}
	 }
}

}

unlink 'temp_interface.txt';
system("rm -r /home/sheng/CCP4/$pisaprojname");
########sub##################

sub interface{
    my ($file,$label)=@_;
    open HH, "$file" or die "No interface files\n";	

    my $mark=0;
    my @id=();
    my @chain=();
    my @sasa=();
    my %bond=();
    my $bondtype='';
    my @num_interfaceresidues=();
    my $chainid='';
	  while(<HH>){
	  	if($_=~/Selection range/){
	  		my @matches=split/[\| \t]+/,$_;
	  		$chainid=$matches[3].$matches[4];
	  	}
	  	elsif($_=~/Residues in the interface/){
	  		  my @matches=$_=~/([0-9\.]+)/g;
	  		  @num_interfaceresidues=($matches[0],$matches[2]);
	  		}
	  	elsif($_=~/Buried ASA/){
	  		  my @matches=$_=~/([\d\.]+)/g;
	  		  @sasa=($matches[1],$matches[3]);
	  		}
	  	 if($_=~/\-\+/){$mark=1;next;}
	  	 if($_=~/\-\'/){$mark=0;$bondtype='';next;}
	  	 if($_=~/Hydrogen/){$bondtype='H';}
	  	 if($_=~/Salt/){$bondtype='SB';}
	  	 if($bondtype){
	  	   	if($_=~/\:/){
	  	   		chomp;
	  	   	  	my @line=split/\|/,$_;
	  	   	    $line[1]=~/([^ \t]+)\:([\w]+)[ \t]+(\d+)/;
	  	   	    my $chain1=$1;
	  	   	    if(!$chain1){next;}
	  	   	    my $aa1=$change{$2};
	  	   	    my $pos1=$3;
	  	   	    $line[3]=~/([^ \t]+)\:([\w]+)[ \t]+(\d+)/;
	  	   	    my $chain2=$1;
	  	   	    my $aa2=$change{$2};
	  	   	    my $pos2=$3;
	  	   	    $bond{$bondtype}.=$chain1.$aa1.$pos1.'_'.$chain2.$aa2.$pos2.';';  	
	  	   	    $bonds{$bondtype}{$chain1.$aa1.$pos1.'_'.$chain2.$aa2.$pos2}{$para{'-frame'}}=1;   	    
	  	   	}
	  	 }
	  	 if($mark==1){
	  	 	chomp;
	  	     my @line=split/[ \t\|]+/,$_;	
	  	     my $a=$line[3];
	  	     my $b=$line[4];
	  	     if($line[3]=~/([0-9A-Za-z])\:/){$a=$1;}
	  	     if($line[4]=~/([0-9A-Za-z])\:/){$b=$1;}
	  	     
	  	     if(($a =~/[$para{'-c1'}]/ && $b =~/[$para{'-c2'}]/)||($a =~/[$para{'-c2'}]/ && $b =~/[$para{'-c1'}]/)){
	  	       #$id=$line[1];	
	  	       push @chain,"$line[3]"."$line[4]";
	  	       push @id,$line[1];
	  	    }
	  	}
	  }
	  if($label){
	 	 	if(!$bond{'SB'}){$bond{'SB'}='NA'}
	  	if(!$bond{'H'}){$bond{'H'}='NA'}
	  	print ASA "$para{'-frame'}\t$chainid\t$num_interfaceresidues[0]\t$num_interfaceresidues[1]\t$sasa[0]\t$sasa[1]\t$bond{'H'}\t$bond{'SB'}\n";
	  }
	  return \@id,\@chain;
}

sub residues{
    my ($file)=@_;
    open HH, "$file" or die "No interface files\n";		  
  	my $mark=0;
  	my $i=0;
  	my %residues=();
  	while(<HH>){
  	   if($_=~/Interfacing Residues/)	{
  	     	 $i=4;next;
  	   }
  	   if($_=~/\'/){$i=0;next;}
  	   if($i>1){$i--;next;}
  	   if($i==1){
  	     my @line=split/[ \t\|]+/,$_;	
  	     my @chain=split/\:/,$line[3];
  	     my $p=substr($_,length('    1 |I| K:NAG'),4);
  	     $p=~s/[ ]+//;
  	     my $type=substr($_,length('    1 |I| K:NAG1041  |'),2);
  	     $type=~s/[ ]+//;
  	     if(!$type){$type='HPH';}
  	      if($line[2] eq 'I'){
  	      	$residues{$chain[0]}{$p}=$type;
  	      }
  	  }
  	}
  	return \%residues;
}


sub pdb_pos{
    my ($pdb)=@_;
  open HH,"$pdb" or die "PDB file not found\n";
  my $seq='';
  my %pos=();
  my $mark=0;
  my $chain='';
  my $start=0;
  my $pos=0;
  foreach(<HH>){
    
    if($_=~/^ATOM/&&substr($_,length('ATOM   1243  '),3) eq "C  "){
    	$pos=substr($_,22,5);
    	$pos=~s/ +//g;
    	if($chain ne substr($_,21,1)){
    	   $chain=substr($_,21,1); 
    	   $start=$pos;   
    	}
       $pos{$chain}{$pos}=$start;
       $start++;
    }
  }
    return %pos;

}
