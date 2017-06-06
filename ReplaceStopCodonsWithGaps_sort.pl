#!/usr/bin/perl -w
##
## use this to remove stop codons from an alignment
## typically, this would be done to calculate dN/dS in HYPHY
## Usage: perl ../Scripts/ReplaceStopWithGaps.pl -nuc 104D5.fasta
## use this to replace stop codons from the nucleotide alignment 

 
use strict;
use Getopt::Long;
use Bio::SeqIO;
 
my ($innuc,$output, $i, %stop);
&GetOptions(
'nuc:s' => \$innuc,
);
  
$output = $innuc;
$output.="_nostop.fasta";
#print "output = $output\n";
my $nuc  = Bio::SeqIO->new(-file => "$innuc" , '-format' => 'fasta');
my $out = Bio::SeqIO->new(-file => ">$output" , '-format' => 'fasta');

my %idhash = ();

while (my $nucseq = $nuc->next_seq()){
	my $nuc_id=$nucseq->id();
	my $sequence_string=uc($nucseq->seq);
	my $n = 3; 
	my @codons = unpack "a$n" x (length( $sequence_string ) /$n ), $sequence_string;
	for my $i (0 ... $#codons){
		$codons[$i] = "---" if $codons[$i] eq "TAA" || $codons[$i] eq "TAG" || $codons[$i] eq "TGA";
	}
	$nucseq->seq(join("",@codons));
	$idhash{$nuc_id} = $nucseq;
#	$out->write_seq($nucseq);
}

my @ids = sort { $a cmp  $b } (keys %idhash);

foreach my $nuc_id (@ids){
	$out->write_seq($idhash{$nuc_id});	
}



