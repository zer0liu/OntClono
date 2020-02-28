#!/usr/bin/perl -w

#use Data::Dumper;
use 5.010;
use strict;
use warnings;

use List::Util qw(shuffle);
use Bio::SeqIO;
use Net::Random;
use File::Basename;
use Smart::Comments;

# This script can be used to Mutate a DNA Sequence,only mismatch.
my $usage = << "EOS";
Generate mutated DNA sequence which contains a number of mismatches,

Usage:
  mutation.pl <fin> <fout> <num_seq> <prefix><num_base>
Args:
  fin	    A FASTA format file. Only the FIRST sequence will be operated.
  fout	    A multi-FASTA sequence file.
  num_seq   Number of sequences to be generated.
            Default 1.
  prefix    Prefix of generated sequence ID. 
            Default 'RND_'.
  num_base	Number of bases to be mutated.			
EOS

my $fin		= shift or die $usage;
my $fout	= shift or die $usage;
my $num_seq = shift // 1;
my $prefix  = shift // 'RND_';
my $num_mut_sites   = shift or die $usage;

# Generate output file basenames for gaps and mismatches
my ($basename, $dir, $suffix)   = fileparse($fout, qr/\..*/);


my $fmis    = $basename . '_mis.txt';

# Get input sequence
my $o_seqi  = Bio::SeqIO->new(
    -file   => $fin,
    -format => 'fasta',
);

my $o_seq   = $o_seqi->next_seq;

# Sequence string
my $seq_str = $o_seq->seq;

# Sequence length
my $seq_len = $o_seq->length;


# Generate output multi-FASTA sequences
# Convert sequence into a hash
my @bases = split //, $seq_str;

my $pos = 0;    # Base location start at '0'
my %base_pos;

for my $base (@bases) {
	$base_pos{$pos} = $base;
	$pos++;
}

# Generate output file
open my $fh_out, ">", $fout
	or die "[ERROR] Create output file '$fout' failed!\n$!\n";



open my $fh_mis, ">", $fmis
    or die "[ERROR] Create output gap file '$fmis' failed!\n$!\n";

## %base_pos
my $rnd_seq_num	= 0;

while ($rnd_seq_num < $num_seq) {
	my %work_base_pos	= %base_pos;
	
	#my @rnd_pos = shuffle(0..($seq_len - 1));
	 my @rnd_pos = qw(30 96 105 106 108 116 193 214 262 265 310 377 395 456 546 569 595 639 648 684 698 700 771 841 848 866 867 877 895 934 1002 1076 1123 1129 1192 1206 1207 1242 1244 1273 1277 1293 1304 1334 1354 1424 1497 1503 1504 1552 1594 1663 1734 1757 1773 1867 1906 1929 1938 1972 2039 2054 2082 2106 2108 2193 2264);
	 my @rnd_mis_pos	= @rnd_pos[0..$num_mut_sites-1];
	
	
	# Modify mismatches
	for my $mis_pos	(@rnd_mis_pos) {
		my $cur_base	= $work_base_pos{$mis_pos};
		
		my $mis_base	= site_mutation($cur_base);
		
		## $mis_pos
		## $cur_base
		## $mis_base
		
		$work_base_pos{$mis_pos}	= $mis_base;
		say $fh_mis $cur_base, "\t", $mis_pos, "\t", $mis_base;
	}
	
	
	my $out_seq;
	
	for my $pos (sort {$a<=>$b} keys %work_base_pos) {
		$out_seq .= $work_base_pos{$pos};
	}
	
	# Output
	print $fh_out ">",$prefix, $rnd_seq_num, "\n";
	print $fh_out $out_seq, "\n";
	
	$rnd_seq_num++;
	
    }

   

close $fh_mis;
close $fh_out;

# Name:	site_mutation
# Desc:	Mutate a nucleotide 
# Args:	A character, one of "ACGT", the current base
# Ret:	A mutated base character

sub site_mutation {
	my ($cur_base)	= @_;
	
	my $bases	= 'AGCT';
	
	# Remove current base
	$bases	=~ s/$cur_base//;
	
	@bases	= split //, $bases;
	
	my @rnd_pos	= shuffle(0..2);
	
	return $bases[ $rnd_pos[0] ];
}
