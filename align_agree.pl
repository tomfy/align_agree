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
   $libdir = $bindir;                # . '/../lib';
   $libdir = abs_path($libdir);	# collapses the bin/../lib to just lib
}
# use lib '/home/tomfy/Orthologger_2014_11_28/lib/';
use lib $libdir;
use lib '/home/tomfy/Orthologger/lib/';

use Algnmnt;
use Column;
use AlignAlign;
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
my $interleaved = 0;

GetOptions(
           'min_nongap=i' => \$min_nongaps,
           'min_frac_nongap=f' => \$min_frac_nongaps,
           'min_alal_frac=f' => \$min_alal_frac,
           'interleaved!' => \$interleaved,
	  );

my $start_time = time();

open my $fhA, "<", "$input_alignment_file_A" or die "couldn't open $input_alignment_file_A for reading. \n";
open my $fhB, "<", "$input_alignment_file_B" or die "couldn't open $input_alignment_file_B for reading. \n";

while (1) {
   print STDERR "Calling next_align(A)\n";
   my ($idlineA, $fasta_A, $idline_famsize_A, $famsize_A) = next_align($fhA);
   my $align_obj_A = Algnmnt->new( $fasta_A, {min_nongap_fraction => $min_frac_nongaps, min_nongap_chars => $min_nongaps} );

   print STDERR "Calling next_align(B)\n";
   my ($idlineB, $fasta_B, $idline_famsize_B, $famsize_B) = next_align($fhB);
   my $align_obj_B = Algnmnt->new( $fasta_B, {min_nongap_fraction => $min_frac_nongaps, min_nongap_chars => $min_nongaps} );

   my $alal_obj = AlignAlign->new($align_obj_A, $align_obj_B, {gg_score => 0.25});
   print  $alal_obj->alal_cols($min_alal_frac);

    print $alal_obj->alal_fasta($min_alal_frac, $interleaved);

#   print STDERR $alal_obj->alal_cols(), "\n";
   last;                   # for now just read one pair of alignments.
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
      if (/^Id/){  #.*fam_size:\s*(\d+)/) {
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

# sub print_array_of_cols{
#    my $array_of_cols = shift;
#    for (@$array_of_cols) {
#       print "$_ \n";
#    }
# }

#   sub print_2_arrays_of_cols{
#      my $array_of_cols_A = shift;
#      my $array_of_cols_B = shift;
#      my $array_of_scores = shift;
#      my $aas = shift;
#      my $min_frac = shift || 0;
#      my $i = 0;
#      for (@$array_of_cols_A) {
#         my $frac = ($aas->[$i] > 0)? $array_of_scores->[$i]/$aas->[$i] : 0;
#         if ($frac >= $min_frac) {
#            my $Astring = $array_of_cols_A->[$i];
#            my $Bstring = $array_of_cols_B->[$i];
#            printf("%4i    ", $i);
#            print $Astring, "   ",
#              #       $Bstring, 
#              "\n";
#            my $BdiffA_string = BdiffA($Astring, $Bstring);
#            print "        ", $BdiffA_string, "   ";
#            printf("%4i  %5.3f \n", , $array_of_scores->[$i], $frac); 
#         }
#         $i++;
#      }
#   }

#   sub BdiffA{
#      my $sA = shift;
#      my $sB = shift;
#      my $BdiffA = '';

#      for (my $i=0; $i<length $sA; $i++) {
#         my $cA = substr($sA,$i,1);
#         my $cB = substr($sB,$i,1);
#         $BdiffA .= ($cA eq $cB)? '.' : $cB;
#      }
#      return $BdiffA;
#   }

