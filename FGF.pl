#!/usr/bin/perl
use strict;
use lib "/home/users/bio1411/perl5/lib/perl5";
use Algorithm::NeedlemanWunsch;

###### Program description
#
# Authors:
#	 Robert Anton Hafthorsson
#	
# Description:
#	This program is used for predicting the Fibrinogen gene in mammals. Using a conserved exon method we create a
#	template to search for the exons of other mammals using a pairwise search algorithm. The program then predicts the protein structure     
#	of the gene.
# 
# Modules:
#	nmw - The pairwise algorithm.
#		
# Procedure:
#	1. Read in the DNA sequence.
#	2. Start searching for exons:
#			1) Use a regular expression to search for exons based on patterns observed from already sequenced mammalian exons.
#			   The regular expression is used first, do to the fact its very fast but it is highly selective.
#			2) Use pairwise algorithm to search the exons if the regular expression fails. This process is slower but more rigorous. 
#			3) Translate the exons.
#			4) Use pairwise algorithm to check the translated protein against a conserved protein sequence from different mammals.
#	3. Create the protein from the highest scoring exons.
#	5. Output the characteristics of protein sequence, DNA sequence and final exons, in a "program.gtf" file.
#
# usage: program.pl (input file)


#################### Input file into program.

my $in = $ARGV[0]; #Read in the file from argv.
chomp $in;
my $out = $in;
$out =~ s/\.\w{3,5}/\.gtf/;  #File format Convertion

my $fullseq;
my $all="";
my $protein="";
my $codon;
my $exon1=""; my $exon2="";my $exon3="";my $exon4=""; my $exon5="";my $exon6="";my $exon7="";my $exon8=""; my$exon9="";my $exon10="";

open my $FH, '<', $in or die "Cannot open $in : $!";
open my $OUT,'>', $out or die "Cannot open $out : $!";

while ( my $read = <$FH> ){
	next if $read =~ />/;
	$read =~ s/\n//g;
	$read =~ s/\t//g;
	$read =~ s/\r//g;
	$fullseq .= $read;
}
	#Exon 1
	if ($fullseq =~ /(ATG[AC][AGC]TT[GT]GTC.+CA[AG]C[GA][TG]G[TC][GC]T[AG][GA][GC][AT])/g) {
		#expected length 72
		$exon1 = $1;
	}elsif( $exon1 == "" ){
		my $exon = "ATGAGTTGGTCCTTGCACCCCCGGAGTTTAATTCTCTNCTTCTGTGCTCTTTTANTGCTCTCTTCAACATGCCTGGCA";
		print STDERR "Searching for exon 1 \n";
		($exon1) = pairwise ($exon, $fullseq);		
	}
	#Exon 2
	if ($fullseq =~ /((CAG)?TA[TC][GA][TC][TG]GC[TC][AC].+AGA[TC][GC][AG]{5}TT[TC])/g) {
		#exp length 45
		$exon2 = $1;
	}elsif ($exon2 == "" ){
		my $exon = "TATGTTGCTACCAGAGACAACTGCTGCATCTTAGATGAAAGATTC";
		print STDERR "Searching for exon 2 \n";
		($exon2) = pairwise ($exon, $fullseq);		
	}
	#Exon 3
	if ($fullseq =~ /(G?G[CT]AGTT[AGT][TC]TG[TC][CT][CT][AG][AG][CT][TCG]ACC.+?[GC][AG][AG][TC][TC][ACG][TC][CT][GA][TGA][TA][AG][ACT][CT][ACG][AG])/g){
		#exp length 184
		$exon3 = $1;	
	}elsif ($exon3 == "" ){
		my $exon = "GGTAGTTANTGCCCAACNACCTGTGGCATNGCAGATTTCCTGNNTNCTTACCAAACCNACGTNGACAATGATCTNCNNACTNTGGAAGACATCTTANATCNNGNTGAAAACANAACNNCAGAAGCCAAAGANCTGATCAAAGCAATCCANGTNNACTANAACCCNGANCAACCNNCAAAGCCAN";
		print STDERR "Searching for exon 3 \n";
		($exon3) = pairwise ($exon, $fullseq);
	}
	#Exon 4
	if ($fullseq =~ /([AG]TA[TGA][GC]ATA[CG]A[AGC][GA][GC]TGC[TC]AC[TC].+GA[GAC][TA]C[AG]A[GC][TACG]AT[TC]CG)/g) {
		#exp length 94
		$exon4 = $1;	
	}elsif ($exon4 == "" ){
		my $exon = "NTATGATAGANNGTGCNACTCAGAAGTCNAAGAAGATGGTAGAAGAAATTATGAAATANGAAGCACTGATATTAACNCATGAGTCAAGTATTCG";
		print STDERR "Searching for exon 4 \n";
		($exon4) = pairwise ($exon, $fullseq);
	}
	#Exon 5
	if ($fullseq =~ /([GA]T[AT]TTT[AG]CA[GA]GA[CA]AT[CA]TA[CT][GA].+[AC]AC[TC]GG[AG]A[GA]AG)/g) {
		#exp length 131
		$exon5 = $1;	
	}elsif ($exon5 == "" ){
		my $exon = "NTATTTNCAGGAAATNTATAATTCAAATAATCANAAGATCNNTAACCTNAAANAGAAGGTNGCCCAGCTTGAAGCACAGTGCCAGGAGCCTTGCAANGACNCTGTGCAAATCCATGANANAACTGGNAAAG";
		print STDERR "Searching for exon 5 \n";
		($exon5) = pairwise ($exon, $fullseq);
	}
	#Exon 6
	if ($fullseq =~ /(ATTG[TC]CA[AG]GA.+C[CT]GT[AG][TC]T[TG]CAGAA[AG])/g) {
		#exp length 134
		$exon6 = $1;
	}elsif ($exon6 == "" ){
		my $exon = 'ATTGTCANGATATTGCCAANAAGGGNGCCAAAGANAGTGGACTTTACTTNATTNGACCTNTGAAAGCTAAGCAGCAGTTCTTAGTNTACTGTGAAATCGANGGNTCTGGAAATGGATGGACNGTNTTNCAGAAG';
		print STDERR "Searching for exon 6 \n";
		($exon6) = pairwise ($exon, $fullseq);
	}
	#Exon 7
	if ($fullseq =~ /(AG[GA][GCA]T[TG]GA[CT]GG.+[AG][GA]TGGCA[AG][GA]A[GC]CAG)/g) {
		#exp length 185
		$exon7 = $1;
	}elsif ($exon7 == "" ){
		my $exon = "AGGCTTGATGGCAGTGTGGATTTCAAGAANAACTGGATTCANTATAAAGAAGGATTTGGACANCTGTCTCCTACTGGCACCACAGANTTTTGGCTGGGAAATGAGAAGATTCATTTGATAAGCANNCAGTCNACCATCCCATATGCANTNAGAATACAGCTNNAAGACTGGAATGGCAGAACCAG";
		print STDERR "Searching for exon 7 \n";
		($exon7) = pairwise ($exon, $fullseq);
	}
	#Exon 8
	if ($fullseq =~ /([TC]AC[TC]GC[ACG]GA[TC]TA[TC][GT]C.+T[GT][GT][TA]T[TC]A[TC]C[AC]AG)/g) {
		#exp length 278
		$exon8 = $1;
	}elsif ($exon8 == "" ){
		my $exon = "CACTGCNGACTATGCCATGTTCANGGTGGGNCCTGAATCTGACAAATACCGCCTGACNTATGCCTACTTCATTGGTGGAGATGCNGGNGATGCCTTNGANGGCTACGATTTTGGCGATGATCCNAGTGACAAGTTTTTCACATCCCANAANGGCATGCAGTTCAGTACCTGGGACAATGACAANGATAAGTTTGAAGGCAACTGTGCTGAACAGGATGGATCTGGNTGGTGGATGAACAANTGTCACGCTGGCCACCTCAATGGAGTTTATTACCAAG";
		print STDERR "Searching for exon 8 \n";
		($exon8) = pairwise ($exon, $fullseq);
	}
	#Exon 9
	if ($fullseq =~ /(G[TG]GG[CTG]ACTTACTC.+[ATCG][ACG][ACG][TG][TGA]{3,6}A[GT]CCAAA[CA]AG)/g) {
		#exp length 170
		$exon9 = $1;
	}elsif ($exon9 == "" ){
		my $exon = "GTGGCACTTACTCAAAATCATCTACTCCTAATGGTTATGANAATGGCATTATTTGGGCCACNTGGAAAACCCGNTGGTATTCCATGAAGNAAACCACCATGAAGATAATNCCNTTCAACAGACTCTCCATTGGAGANGGACAGCANCANCACNTGGGGGGANCCAAACAG";
		print STDERR "Searching for exon 9 \n";
		($exon9) = pairwise ($exon, $fullseq);
		
	}
	#Exon 10
	if ($fullseq =~ /(G[CT]TGGA[AG]A[TC][GA][GTC][TGCA]T[TAG][GAT])/g) {
		##exp length 15
		$exon10 = $1;
	}elsif ($exon10 == "" ){
		my $exon = "GCTGGAGACGTTTAA";
		print STDERR "Searching for exon 10 \n";
		($exon10) = pairwise ($exon, $fullseq);	
	}

################### Search Subroutine
sub pairwise {
my $exonIN = $_[0];  #Take in the exon to check 
my $seqIN = $_[1];
my %scores;
my %scores2;
my %sectors;
my %seq;
my $score;
my @unarr;
my $max;
my $highsector;
my $hsp;
my $highsector;
my @keys;
my $temp2;
my $seqence;	
#First see if i can find it in a big section. 
my $lengd = (length $exonIN);

if ($lengd>30){
	for(my $i=0; $i<(length $seqIN);$i+=30){
		my $temp= substr($seqIN,$i,(length $exonIN));
		($score) = nmw($temp,$exonIN);
			if($score>$max){
				$sectors{$score}=$i;
				$scores{$score}=$temp;
			}
			push @unarr,$score;
			my @sortarr = sort {$a <=> $b} @unarr;
			my $max = int ($sortarr[-1]/2);
			if ($max<0){$max = ($sortarr[-1])-10;}
	}
	my @keys = sort {$a <=> $b}(keys(%scores));
	$hsp =  $keys[-1];
	$highsector = $sectors{$hsp};
	my $temp2 = substr($seqIN,($highsector - 30),((length $exonIN)*2)); 
	my %scores2;
	for(my $j=0; $j<(length $temp2);$j++){
		my $temp3= substr($temp2,$j,(length $exonIN));
		($score) = nmw($temp3,$exonIN);
		if($score>$max){
			my $temp4 = substr($temp3,0,6);
			my $temp5 = substr($exonIN,0,6);
			my ($score2) = nmw ($temp3,$exonIN);
			if (($score2>30) and ($score>$max)){
				$scores2{$score2}=$temp3;
			}
		}
		if ($score<-10){
			last;
		}	
		} 
	my @keys2 = sort {$a <=> $b}(keys(%scores2));
	my $hsp2 =  $keys2[-1];
	my $seqence = $scores2{$hsp2};
	return $seqence;	
}
else{
	my $max;
	for(my $i=0; $i<(length $seqIN);$i++){
		my $temp= substr($seqIN,$i,(length $exonIN));
		($score) = nmw($temp,$exonIN);
		$scores{$score}=$temp;
		push @unarr,$score;
		my @sortarr = sort {$a <=> $b} @unarr;
		$max = int ($sortarr[-1]/2);
		if ($max<0){
			$max = ($sortarr[-1])-10;
		}
	}
	my @keys2 = sort {$a <=> $b}(keys(%scores));
	my $hsp2 =  $keys2[-1];
	my $seqence = $scores{$hsp2};
	return $seqence;
}

}
################### Needleman Wunsch
sub nmw {

#Here I create some stuff to align. 
my $a = @_[0];
my @a = (split '', $a);
my $b = @_[1];
my @b = (split '', $b);

sub score_sub {
        if (!@_) {
            return -1; # gap penalty
        }
				##mismatch scores -1, match +1
        return ($_[0] eq $_[1]) ? 1 : -1;
}
## callbacks that print something useful
## prints an 'alignment string' in the order of the  
## recursion of the dynamic programming algorithm 
## print "-" only on match
sub on_align  {  "align", " " , $a[$_[0]], ($a[$_[0]] eq $b[$_[1]] ) ? "-" : " ", $b[$_[1]], "\n" }; 
sub on_shift_a { "gap  ", "" , $a[$_[0]], "\n" };
sub on_shift_b {  "gap  ", "   " , $b[$_[0]], "\n"};

### Dumb select, need to return one of the keys for alternative 
### alignments with equal score. Here, we always take the first option, but don't print it.

sub on_select_align {"print (select_align)\n"; return (keys (%{$_[0]})) [0]};
## one gets the same behaviour with not assigning on_select_align, I am too lazy to implement this callback correctly ...
	
	
    my $matcher = Algorithm::NeedlemanWunsch->new(\&score_sub);
    my $score = $matcher->align(
               \@a,
               \@b,
               {   align     => \&on_align,
                   shift_a => \&on_shift_a,
                   shift_b => \&on_shift_b,
                   #select_align => \&on_select_align #This is not needed for us i thinks.
               });
			   
	return $score;
}


################### translation tables.
sub codon2aa{
my($codon)=@_;
$codon=uc $codon;
my %aacode = (
  TTT => "F", TTC => "F", TTA => "L", TTG => "L",
  TCT => "S", TCC => "S", TCA => "S", TCG => "S",
  TAT => "Y", TAC => "Y", TAA => "", TAG => "",
  TGT => "C", TGC => "C", TGA => "", TGG => "W",
  CTT => "L", CTC => "L", CTA => "L", CTG => "L",
  CCT => "P", CCC => "P", CCA => "P", CCG => "P",
  CAT => "H", CAC => "H", CAA => "Q", CAG => "Q",
  CGT => "R", CGC => "R", CGA => "R", CGG => "R",
  ATT => "I", ATC => "I", ATA => "I", ATG => "M",
  ACT => "T", ACC => "T", ACA => "T", ACG => "T",
  AAT => "N", AAC => "N", AAA => "K", AAG => "K",
  AGT => "S", AGC => "S", AGA => "R", AGG => "R",
  GTT => "V", GTC => "V", GTA => "V", GTG => "V",
  GCT => "A", GCC => "A", GCA => "A", GCG => "A",
  GAT => "D", GAC => "D", GAA => "E", GAG => "E",
  GGT => "G", GGC => "G", GGA => "G", GGG => "G",
); # this is the hash table for the amino acids

if(exists $aacode{$codon})
{
return $aacode{$codon};
}
}
my @prot = (
	"MSWSLQPRSFILCXALLLLSPTGLA", #Exon 1
	"QYVATRDNCCILDERF", #Exon 2
	"GSYCPTTCGIADFLXSYQTDVDKDLXTLEDILXXXENXTXEAKELIKAIQVYYNPDQPXKP", #Exon 3
	"MIXSATQKSKKMVEEIMKYEALXLTHESSI", #Exon 4
	"YLQEIYNSNNQKIXNLKQKVAQLEAQCQEPCKDXVQIHDTTGK", #Exon 5
	"CQDIANKGAKESGLYFIRPLKAKQQFLVYCEIDGSGNGWTVXQK", #Exon 6
	"RLDGSXDFKKNWIQYKEGFGHLSPTGTTEFWLGNEKIHLISXQSTIPYALRIQLXDWNGRT", #Exon 7
	"TADYAMFXVGPESDKYRLTYAYFIGGDAGDAFDGYDFGDDPSDKFFTSHNGMQFSTWDNDNDKFEGNCAEQDGSGWWMNKCHAGHLNGVYYQ", #Exon 8
	"GTYSKSSTPNGYDNGIIWATWKTRWYSMKETTMKIIPFNRLSIGEGQQHHMGGSKQ", #Exon 9
	"AGDV" #Exon 10
);
################### For printing 

$all = $exon1.$exon2.$exon3.$exon4.$exon5.$exon6.$exon7.$exon8.$exon9.$exon10; # Here i concatonate all the exons.
my @exons =($exon1,$exon2,$exon3,$exon4,$exon5,$exon6,$exon7,$exon8,$exon9,$exon10); # I create a referance to the exons in an array.
my $exon_ref = \@exons;
my @phase;
my @Exonsphase;
################### Findig the phase and translating
my $score;
print "Exon nr\t\tScore\tORF\tAmino Acid\n";
print $OUT "Prediction information table\nExon nr\t\tScore\tORF\tAmino Acid\n";

for (my $j=0; $j<(scalar @exons); $j++){
	my $temp='';
	### ORF +1
	for(my $i=0;$i<(length($$exon_ref[$j])-2);$i+=3){
		$codon=substr($$exon_ref[$j],$i,3);
		$temp.= &codon2aa($codon);	
	}
	($score) = nmw($temp, $prot[$j]);
	if ($score >=0){
		$phase[$j] = 0;
		my $phaseF='';
		my $phaseF1='';
		my $phaseF2='';
		my $phaseprotF='';
		if ($phase[$j-1] == 2 and $phase[$j] == 0){
			$phaseF1 = substr($$exon_ref[$j-1],((length $$exon_ref[$j-1])-2),2);
			$phaseF2 = substr($$exon_ref[$j+1],0,1);
			$phaseF = $phaseF1.$phaseF2;
			if(length $phaseF == 1 or length $phaseF == 2){
				$phaseF='';
			}else{
				$phaseprotF.=&codon2aa($phaseF);
			}
			$protein.=$temp;
			$protein.=$phaseprotF;
			my $protein2=$temp.$phaseprotF;
			push @Exonsphase, $protein2;
		}else {
			$protein.=$temp;
			my $protein2=$temp;
			push @Exonsphase, $protein2;
		}
		print "Exon [",$j+1,"]\t$score\t+1\t$phaseF   $phaseprotF\n";
		print $OUT "Exon [",$j+1,"]\t$score\t\t+1\t$phaseF   $phaseprotF\n";
	}
	elsif($score<0){
		### ORF +2
		$temp='';
		for(my $i=1;$i<(length($$exon_ref[$j])-1);$i+=3){
		$codon=substr($$exon_ref[$j],$i,3);
		$temp.=&codon2aa($codon);
	}
	($score) = nmw($temp, $prot[$j]);
		if ($score>=0){
			$phase[$j] = 1;
			my $phaseF;
			my $phaseF1;
			my $phaseF2;
			my $phaseprotF;
			if($phase[$j]==1 and $phase[$j-1]==2){
				$phaseF1 = substr($$exon_ref[$j-1],((length $$exon_ref[$j-1])-2),2);
				$phaseF2 = substr($$exon_ref[$j],0,1);
				$phaseF = $phaseF1.$phaseF2;
				$phaseprotF.=&codon2aa($phaseF);
				$protein.=$phaseprotF;
				$protein.=$temp;
				my $protein2=$temp.$phaseprotF;
				push @Exonsphase, $protein2;
			}else{
			$protein.=$temp;
			my $protein2=$temp;
			push @Exonsphase, $protein2;
			}
			print "Exon [",$j+1,"]\t$score\t+2\t$phaseF   $phaseprotF\n";
			print $OUT "Exon [",$j+1,"]\t$score\t\t+2\t$phaseF   $phaseprotF\n";
		}elsif($score<0){
			### ORF +3
			$temp='';
			for(my $i=2;$i<(length($exons[$j]));$i+=3){
				$codon=substr($$exon_ref[$j],$i,3);
				$temp.=&codon2aa($codon);
			}
			($score) = nmw($temp, $prot[$j]);
			if ($score>=0){
				$phase[$j] = 2;
				my $phaseF;
				my $phaseF1;
				my $phaseF2;
				my $phaseprotF;
				if ($phase[$j]==2 and $phase[$j-1]==0){
					$phaseF1 = substr($$exon_ref[$j-1],((length $$exon_ref[$j-1])-1),1);
					$phaseF2 = substr($$exon_ref[$j],0,2);
					$phaseF = $phaseF1.$phaseF2;
					$phaseprotF.=&codon2aa($phaseF);
					$protein.=$phaseprotF;
					$protein.=$temp;
					my $protein2=$temp.$phaseprotF;
					push @Exonsphase, $protein2;
				}
				elsif($phase[$j]==2 and $phase[$j-1]==1){
					$phaseF1 = substr($$exon_ref[$j-1],((length $$exon_ref[$j-1])-1),1);
					$phaseF2 = substr($$exon_ref[$j],0,2);
					$phaseF = $phaseF1.$phaseF2;
					$phaseprotF.=&codon2aa($phaseF);
					$protein.=$phaseprotF;
					$protein.=$temp;
					my $protein2=$temp.$phaseprotF;
					push @Exonsphase, $protein2;
				}else {
				$protein.=$temp;
				my $protein2=$temp;
				push @Exonsphase, $protein2;
				}
				print "Exon [",$j+1,"]\t$score\t+3\t$phaseF   $phaseprotF\n";		
				print $OUT "Exon [",$j+1,"]\t$score\t\t+3\t$phaseF   $phaseprotF\n";		
			}else{
				$protein.="";
				my $protein2="";
				push @Exonsphase, $protein2;
			}
		}
	 }
 }
my $length = length $protein;
$protein =~ s/(.{80})/$1\n/g;
$all =~ s/(.{80})/$1\n/g;

if ($length<400){
print "No prediction was made. Check if dna is for a Fibrinogen Gamma chain and try again\n";
}else{
print $OUT "\nThe complete protein and length: ",$length,"";
print $OUT "\n$protein\n\nThe Translated Exons\n"; 
for(my $i= 0; $i < scalar@Exonsphase;$i++){print $OUT "Exon nr[",$i+1,"]\t:\t",$Exonsphase[$i],"\n"};
print $OUT "\nThe DNA sequence and length :",length $all,"\n$all\n";
print $OUT "\nThe Exons Shown Separately:\nExon 1\t:\t$exon1\nExon 2\t:\t$exon2\nExon 3\t:\t$exon3\nExon 4\t:\t$exon4\nExon 5\t:\t$exon5\nExon 6\t:\t$exon6\nExon 7\t:\t$exon7\nExon 8\t:\t$exon8\nExon 9\t:\t$exon9\nExon 10\t:\t$exon10\n";
}

close $FH;
close $OUT; 
