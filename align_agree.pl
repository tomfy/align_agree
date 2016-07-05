#!/usr/bin/perl -w
use strict;
use Getopt::Long;
use List::Util qw (min max sum);

use File::Basename 'dirname';
use Cwd 'abs_path';
my ( $bindir, $libdir );

BEGIN {	    # this has to go in Begin block so happens at compile time
   $bindir =
     dirname( abs_path(__FILE__) ) ; # the directory containing this script (i.e. Orthologger/bin )
   $libdir = $bindir . '/../lib';
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
# use lib '/home/tomfy/Orthologger_2014_11_28/lib/';
use lib $libdir;
use lib '/home/tomfy/Orthologger/lib/';

use Algnmnt;
use Column;

use TomfyMisc qw 'timestring ';

# compare two alignments, A and B, of the same set of sequences (e.g. A done with mafft, B with muscle)
# for each pair of cols (one from alignment A, the other from B) a score is defined 
# (for now just add 1 for each sequence where the two cols agree (same AA or both gaps), 0 if disagree)
# try to find the optimal correspondence between the cols of the two alignments.
# If alignment A has N_A cols, and B is longer with N_A + M cols, then there are (N_A + M)!/(N_A! * M!)
# possible ways of lining the two alignments up by inserting M cols of gaps in A, or equivalently by 
# deleting M cols from B.
# The extreme cases are all M inserted cols could be at the beginning (so A's col i gets matched with B's
# col i+M), or all M inserted cols are at the end, (so A's col i gets matched with B's col i). So we need to
# consider matching A's col i with i, i+1, ... i+M, etc. (M+1 possibilities)
# So, read in two alignments, decide which is shorter, and

# read in a file with alignments for each of many families
# The format is for each family:
# 1) 1 line of info about family as whole
# 2) fasta for all aligned sequences in family (this is omitted in cases of families which do not
#    meet certain criteria; at present by default if they don't include sequences from >= 3 monocot species
# 3) one blank line

my $input_alignment_file_A = shift;
my $input_alignment_file_B = shift;
my $min_nongaps = 10;
my $min_frac_nongaps = 0.05;
my $min_alal_frac = 0.0;
# my $way = 'string';

GetOptions(
           'min_nongap=i' => \$min_nongaps,
           'min_frac_nongap=f' => \$min_frac_nongaps,
           'min_alal_frac=f' => \$min_alal_frac,
	  );

my $start_time = time();

open my $fhA, "<", "$input_alignment_file_A" or die "couldn't open $input_alignment_file_A for reading. \n";
open my $fhB, "<", "$input_alignment_file_B" or die "couldn't open $input_alignment_file_B for reading. \n";
while (1) {
   my ($idlineA, $fasta_A, $idline_famsize_A, $famsize_A) = next_align($fhA);
   my $align_obj_A = Algnmnt->new( $fasta_A, $min_frac_nongaps, $min_nongaps );

   my ($idlineB, $fasta_B, $idline_famsize_B, $famsize_B) = next_align($fhB);
   my $align_obj_B = Algnmnt->new( $fasta_B, $min_frac_nongaps, $min_nongaps );

   my $col_objects_A = $align_obj_A->array_of_col_objects();
   my $col_objects_B = $align_obj_B->array_of_col_objects();

   my ($alalA, $alalB, $scores, $aas) = NeedlemanWunsch($col_objects_A, $col_objects_B);
   print_2_arrays_of_cols($alalA, $alalB, $scores, $aas, $min_alal_frac);

   last;                    # for now just read one pair of aligments.
}
my $end_time = time();
print "# ended at: ", timestring($end_time) , ". Run time ", $end_time - $start_time, " seconds. \n";

###########################################################################

# read in next alignment from $fh (alfastas file format);
# expect Id ... line, then fasta lines (w gaps),
# then blank line (or EOF?)
# if there is no next alignment (i.e. hit EOF before Id ... line), then returns ('', '', undef, 0)
sub next_align{
   my $fh = shift;
   my $id_line = '';
   my $idline_famsize;
   my $count_famsize= 0;
   my $fasta = '';
   while (<$fh>) {
      next if(/^\s*#/);
      if (/^Id.*fam_size:\s*(\d+)/) {
         $id_line = $_;
         $idline_famsize = $1;
         last;
      } else {
         warn "Expected next line to start with 'Id ', but didn't.\n";
      }
   }
   while (<$fh>) {
      last if(/^\s*$/);         # blank line - done
      $count_famsize++ if(/^>/);
      $fasta .= $_;
   }
   return ($id_line, $fasta, $idline_famsize, $count_famsize);
}

sub hd{
   return (($_[0] ^ $_[1]) =~ tr/\001-\255//); #  ^ will be non-zero iff chars differ, tr counts non-zeroes.
}

sub print_array_of_cols{
   my $array_of_cols = shift;
   for (@$array_of_cols) {
      print "$_ \n";
   }
}

  sub print_2_arrays_of_cols{
     my $array_of_cols_A = shift;
     my $array_of_cols_B = shift;
     my $array_of_scores = shift;
     my $aas = shift;
     my $min_frac = shift || 0;
     my $i = 0;
     for (@$array_of_cols_A) {
        my $frac = ($aas->[$i] > 0)? $array_of_scores->[$i]/$aas->[$i] : 0;
        if($frac >= $min_frac){
           my $Astring = $array_of_cols_A->[$i];
           my $Bstring = $array_of_cols_B->[$i];
        printf("%4i    ", $i);
        print $Astring, "   ", $Bstring, "\n";
           my $BdiffA_string = BdiffA($Astring, $Bstring);
        print "        ", $BdiffA_string, "   ";
        printf("%4i  %5.3f \n", , $array_of_scores->[$i], $frac); 
     }
        $i++;
     }
  }

    sub BdiffA{
       my $sA = shift;
       my $sB = shift;
       my $BdiffA = '';

       for (my $i=0; $i<length $sA; $i++) {
          my $cA = substr($sA,$i,1);
          my $cB = substr($sB,$i,1);
          $BdiffA .= ($cA eq $cB)? '.' : $cB;
       }
       return $BdiffA;
    }

  sub NeedlemanWunsch{
     my $A1 = shift;            # array ref of col objects
     my $A2 = shift;

     my $L1 = scalar @$A1;
     my $L2 = scalar @$A2;
     my $Nseq = length $A1->[0]->{col_str}; # print "Nseq: $Nseq \n"; exit;

     my $Gap = 0;               # gap penalty

     # Initialize:
     my @matrix;
     $matrix[0][0]{score} = 0;
     $matrix[0][0]{pointer} = 'none'; # upper-left corner - nowhere to go.
     for (my $j=1; $j <= $L1; $j++) {
        $matrix[0][$j]{score} = $Gap * $j;
        $matrix[0][$j]{pointer} = 'left';
     }
     for (my $i=1; $i <= $L2; $i++) {
        $matrix[$i][0]{score} = $Gap * $i;
        $matrix[$i][0]{pointer} = 'up';
     }

     # Fill:
     for (my $i=1; $i <= $L2; $i++) {
           my $string2 = $A2->[$i-1]->{col_str};
           $string2 =~ s/-/~/g;
        for (my $j=1; $j <= $L1; $j++) {
           my $string1 = $A1->[$j-1]->{col_str};
           my ($diag_score, $left_score, $up_score);

           $diag_score = $matrix[$i-1][$j-1]{score} + $Nseq - hd($string1, $string2);

           $up_score = $matrix[$i-1][$j]{score} + $Gap;
           $left_score = $matrix[$i][$j-1]{score} + $Gap;

           if ($diag_score >= $up_score) {
              if ($diag_score >= $left_score) {
                 $matrix[$i][$j]{score} = $diag_score;
                 $matrix[$i][$j]{pointer} = 'diag';
              } else {
                 $matrix[$i][$j]{score} = $left_score;
                 $matrix[$i][$j]{pointer} = 'left';
              }
           } else {
              if ($up_score >= $left_score) {
                 $matrix[$i][$j]{score} = $up_score;
                 $matrix[$i][$j]{pointer} = 'up';
              } else {
                 $matrix[$i][$j]{score} = $left_score;
                 $matrix[$i][$j]{pointer} = 'left';
              }
           }
        }
     }

     # Traceback:
     my @align1 = ();
     my @align2 = ();
     my @scores = ();
     my @aas = ();
     my $j = $L1;
     my $i = $L2;

     my $gap_col_string = '****************';
     while (length $gap_col_string < $Nseq) {
        $gap_col_string .= $gap_col_string;
     }
     $gap_col_string = substr($gap_col_string, 0, $Nseq);

     while (1) {
        if ($matrix[$i][$j]{pointer} eq 'diag') {
           unshift @align1, $A1->[$j-1]->{col_str};
           unshift @align2, $A2->[$i-1]->{col_str};
           my $s2 = $A2->[$i-1]->{col_str};
           $s2 =~ s/-/~/g;
           unshift @scores, $Nseq - hd($A1->[$j-1]->{col_str}, $s2);
           unshift @aas, min($A1->[$j-1]->get_nongap_count('-'), $A2->[$i-1]->get_nongap_count('-'));
           $i--; $j--;
        } elsif ($matrix[$i][$j]{pointer} eq 'left') {
           unshift @align1, $A1->[$j-1]->{col_str};
           unshift @align2, $gap_col_string;
           unshift @scores, 0;
           unshift @aas, 1;
           $j--;
        } elsif ($matrix[$i][$j]{pointer} eq 'up') {
           unshift @align1, $gap_col_string;
           unshift @align2, $A2->[$i-1]->{col_str}; #substr($seq2, $i-1, 1);
           unshift @scores, 0;
           unshift @aas, 1;
           $i--;
        } elsif ($matrix[$i][$j]{pointer} eq 'none') {
           last;
        } else {
           die "Unknown NW pointer value: " . $matrix[$i][$j]{pointer}, "\n";
        }
     }
     map(s/~/-/g, @align2);     # restore - as gap character.
     return (\@align1, \@align2, \@scores, \@aas);
  }

#########################################################

