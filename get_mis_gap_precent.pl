#!/usr/bin/perl


use 5.010;
use strict;
use warnings;

use Bio::AlignIO;
use Bio::SeqIO;
use Smart::Comments;

my $usage = << "EOS";
根据给出的位置信息，计算序列的插入、缺失和错配；
Usage:
  get_mis_gap_precent.pl <fblast> <fin_query> <fsite> <fout>
Args:
  <fblast>  NCBI BLAST+ report file. "-outfmt 2"
  <fsite>   Site locations
  <fout>	Output 
Note:
  1. This script assumes these is ONLY ONE BLAST result in the report.
EOS

my $fblast  = shift or die $usage;
my $fsite   = shift or die $usage;
my $fout	= shift or die $usage;

# Generate output file basenames for gaps and mismatches
my $fdelet    = $fout . '_delet.txt';
my $finsert  = $fout . '_insert.txt';
my $fmis    = $fout . '_mis.txt';

open(my $fh_blt, "<", $fblast)
    or die "[ERROR] Open BLAST report failed!\n$!\n";

# Get all sequence IDs
my $ra_seqids	= get_all_seqids($fh_blt);

my $seq_num		= scalar( @{ $ra_seqids } );

say "[NOTE] Total sequence number:\t$seq_num\n";


# Initialize control hash, to record hit sequence status
my %seq_status	= ();

$seq_status{ $_ } = 0 for ( @{ $ra_seqids } );

# for my $seqid ( @{ $ra_seqids } ) {
# 	$seq_status{ $seqid }	= 0;
# }

# Initialize alignment hash
my %aln	= ();

$aln{ $_ } = '' for ( @{ $ra_seqids } );

# for my $seqid ( @{ $ra_seqids } ) {
# 	$aln{ $seqid }	= '';
# }

# Query status ALWAYS "1"
$seq_status{"Query_1"}	= 1;

# my %blast_info  = ();   # Blast information
# my %aln = ();   # Alignment

# Reset file handle
# seek($fh_blt, 0, SEEK_SET);
seek($fh_blt, 0, 0);


# Generate output file
open my $fh_delet, ">", $fdelet
    or die "[ERROR] Create output gap file '$fdelet' failed!\n$!\n";

open my $fh_insert, ">", $finsert
    or die "[ERROR] Create output gap file '$finsert' failed!\n$!\n";
    
open my $fh_mis, ">", $fmis
    or die "[ERROR] Create output gap file '$fmis' failed!\n$!\n";

my %ref_gap = ();

# Parse file and get alignment
while (<$fh_blt>) {
	if (/^(?:BLASTN|BLASTP|BLASTX|TBLASTN|TBLASTX)/) {
		chomp;
		say "[NOTE] ", $_;
	}
    elsif (/^Database:\s(.+)$/) {  # BLAST database
        # $blast_info{'db'}   = $1;
		say "[NOTE] BLAST database: ", $1;
    }
    elsif (/^Query= (.+)$/) {   # Query 
        # $blast_info{'qry'}  = $1;
		say "[NOTE] Query: ", $1;
    }
    elsif (/^Length=(\d+)$/) {  # Query length
        # $blast_info{'qry_len'}  = $1;
		say "[NOTE] Length: ", $1;
    }
	# Begin to parse alignment block
    elsif (/^Query_1/) {
		# Init sequence status hash
		my %chk_seq	= %seq_status;
		
        # if (/^Query_1\s+[\d\s]{6}([ACGT\-\s]+?)\s{2}\d+$/) {
		if (/^(Query_1\s+)(\d+\s+)([ACGT\-]+\s*)(\s{2}\d+)$/i) {
			# Split Query line into *4* parts
            my $seqid  		= $1;
			my $aln_start	= $2;
            my $aln_seq 	= $3;
			my $aln_end		= $4;
			
			# Get indexes of each part
			# Index for sequence ID is ALWAYS 0
			my $idx_seqid		= 0;
			my $idx_aln_start	= length($seqid);
			my $idx_aln			= length($seqid) + length($aln_start);
			my $idx_aln_end		= length($seqid) + length($aln_start) + length($aln_seq);
			
			$seqid			=~ s/\s//g;
			
			
			$aln_seq	=~ s/\s//g;	# Remove ANY spaces
			my $seq		= $aln_seq;			
			
			my $aln_len	= length($seq);	# Length of current alignment
			

            $aln{ $seqid }   .= $seq;

            # Read whole alignment block
            while (<$fh_blt>) {
				##$_
				
                last if /^\s*$/;	# A blank line is the END of an alignment BLOCK
				#next if /^\s+$/;		# Skip InDel lines
				
				$ref_gap{$seqid}=0 unless ($ref_gap{$seqid});

				if (/^\s+/) {
					my $row = $_;
					if (my $m = ( $row =~ s/[ACGT]//gi ) ) {
					    $ref_gap{$seqid} = $ref_gap{$seqid} + $m;				
					}
					next;
				}

				chomp;
				
				$seqid	= substr($_, 0, $idx_aln_start);
				$seq		= substr($_, $idx_aln, $idx_aln_end - $idx_aln);
				
			
				# Remove any spaces in sequence ID
				$seqid		=~ s/\s//g;
				## $seqid
				
				# Replace possible spaces in alignment sequence with 'N'
				$seq		=~ s/\s/N/g;
				
				$aln{ $seqid }	.= $seq;
				
				$chk_seq{ $seqid }	= 1;
            }
			
			# Fill sequences not existed in current alignment with 'N' x aln_len
			my $fill_seq	= 'N' x $aln_len;
			
			for my $seqid (keys %chk_seq) {
				next if ($chk_seq{ $seqid } == 1);
				
				# Sequence status NOT SET, i.e., 0
				$aln{ $seqid }	.= $fill_seq;
			}
        }
        else {
            die "[ERROR] Wrong Query line!\n";
        }
    }
}

close $fh_blt;

my $ftemp    = $fout . '_temp.fasta';
# Output to a temporary multi-FASTA alignment file
open my $fh_tmp, ">", $ftemp
    or die "[ERROR] Create temporary alignment file failed!\b$!\n";

for my $seq_id (keys %aln) {
    say $fh_tmp ">", $seq_id;
    say $fh_tmp $aln{$seq_id};
}

close $fh_tmp;

# Read site location file
open my $fh_site, "<", $fsite
	or die "[ERROR] Read input site location file failed!\n$!\n";

my @sites	= ();
	
while (<$fh_site>) {
	next if /^#/;
	next if /^\s*$/;
	chomp;
	
	push @sites, $_;
}

close $fh_site;

my $o_alni	= Bio::AlignIO->new(
	-file	=> $ftemp,
	-format	=> 'fasta',
);

my $o_aln	= $o_alni->next_aln;	# There is only ONE alignment in the file

my %site_aln	= ();

#sequence length
my $seq_len = @sites;


for my $site ( @sites ) {
	# my $o_slice;
	
	my $o_slice	= $o_aln->slice($site, $site, 1);	# Get the whole column of '$site'

	for my $o_seq ( $o_slice->each_seq) {
		my $seqid	= $o_seq->id;
		my $seq		= $o_seq->seq;
		
		$site_aln{ $seqid }	.= $seq;
	}
}

my $fnon_haplo    = $fout . '_non_haplo.fasta';
# Output FASTA format to STDOUT
open my $fh_out, ">", $fnon_haplo
	or die "[ERROR] Create output haplotype file failed!\n$!\n";

for my $seqid (keys %site_aln) {
	say $fh_out ">", $seqid;
	say $fh_out $site_aln{$seqid};
}

close $fh_out;
=pod
# Get input query sequence
my $query_seqi  = Bio::SeqIO->new(
    -file   => $fin_query,
    -format => 'fasta',
);
my $query_seq   = $query_seqi->next_seq;
# Sequence string

=cut
# Get query sequence
my $data_seqi  = Bio::SeqIO->new(
    -file   => $fnon_haplo,
    -format => 'fasta',
);

my $query_seq;

while (my $data_seq   = $data_seqi->next_seq) {
	my $data_seq_str = $data_seq->seq;
	
	my $data_seq_id = $data_seq->id;
	
	if ($data_seq_id =~ /Query_1/){
        $query_seq = $data_seq_str;
        last;
	}else{
        next;
	}	
}

# Convert sequence into a 
my @query_bases = split //, $query_seq;	

# Get reads sequence
$data_seqi  = Bio::SeqIO->new(
    -file   => $fnon_haplo,
    -format => 'fasta',
);


while (my $data_seq   = $data_seqi->next_seq) {
    my %results = ('N' => 0, 'gap' => 0, 'mismatch' => 0);
    
	my $data_seq_str = $data_seq->seq;
	my $data_seq_id = $data_seq->id;
	
	my @data_bases = split //, $data_seq_str;
	
	for (my $idx=0; $idx < $seq_len; $idx++) {
		if ($query_bases[$idx] eq $data_bases[$idx]) {
			next;
		}
		elsif ($data_bases[$idx] eq '-') {
			$results{'gap'}++;
		}
		elsif ($data_bases[$idx] eq 'N') {
			$results{'N'}++;
		}
		else {
			$results{'mismatch'}++;
		}
	}

	## $seq_len
    ## %results
    ## %ref_gap
	say $fh_mis $data_seq_id, "\t", $results{'mismatch'}/($seq_len - $results{'N'} +$ref_gap{$data_seq_id})*100;
	say $fh_insert $data_seq_id, "\t", $results{'gap'}/($seq_len - $results{'N'}+$ref_gap{$data_seq_id})*100;
    say $fh_delet $data_seq_id, "\t", $ref_gap{$data_seq_id}/($seq_len - $results{'N'}+$ref_gap{$data_seq_id})*100;
}
close $fh_mis;
close $fh_delet;
close $fh_insert;


say "[DONE]";

exit 0;

sub get_all_seqids {
	my ($fh)	= @_;
	
	my %seqids;
	
	# Parse alignment BLOCK
	# This BLOCK starts with "Query_" and end with a blank line
	while (<$fh>) {
		if (/^Query_1/) {
			$seqids{'Query_1'}	= '';
			
			while (<$fh>) {
				last if /^\s*$/;
				
				if (/^(\w+)\s/) {
					$seqids{ $1 }	= '';
				}
			}
			
		}
	}
	
	my @srt_seqids	= sort(keys %seqids);
	
	return \@srt_seqids;
}
