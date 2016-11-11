package AlignAlign;
use strict;
use List::Util qw ( min max sum );

# The idea here is the object is constructed from two alignments of the same data,
# e.g. one done with mafft and another done with muscle.
# This class knows how to align the two objects with each other,
# using the Needleman-Wunsch algorithm
# and then output the alignment of the alignments in various
# formats.
# The scoring of the matches between pairs of alignment columns
# is controlled by several parameters

my $NO_ID = 'NOT_AN_ACTUAL_ID';

sub  new {
   my $class = shift;
   my $al_obj_a = shift;
   my $al_obj_b = shift;
   my $args = shift;
   my $default_args = {
                       # aa score is 1
                       gg_score => 0, # gap-gap score
                       ab_score => 0, # score for dissimilar amino acids
                       ag_score => 0, # score for amino-acid-gap
                      };
   my $self = bless $default_args, $class;

   $self->{align_obj_a} = $al_obj_a;
   $self->{align_obj_b} = $al_obj_b;
   $self->{ids} = $al_obj_a->{ids};
   #  print "Algnmnt object fields: \n";
   for my $option ( keys %$args ) {
      warn "Unknown option: $option in Algnmnt constructor.\n" if ( !exists $self->{$option} );
      if ( defined $args->{$option} ) { # if arg is undef, leaves default in effect
         $self->{$option} = $args->{$option};
         print STDERR "setting $option to ", $self->{$option}, "\n";
      }
   }
   $self->NeedlemanWunsch();
   return $self;
}

sub NeedlemanWunsch{
   my $self = shift;
   my $A1 = $self->{align_obj_a}->{cols}; # array ref of col objects
   my $A2 = $self->{align_obj_b}->{cols}; 

   my @ids1 = @{$self->{align_obj_a}->{ids}};
   my @ids2 = @{$self->{align_obj_b}->{ids}};
   my $ids12_agree = (scalar @ids1 == scalar @ids2);
   if ($ids12_agree) {
      for (0..scalar @ids1-1) {
         if ($ids1[$_] ne $ids2[$_]) {
            $ids12_agree = 0;
            last;
         }
      }
   }
   die "Id lists for the two alignments do not agree!\n" if(! $ids12_agree);

   my $L1 = scalar @$A1;                  # length of alignment 1
   my $L2 = scalar @$A2;                  # length of alignment 2
   my $Nseq = length $A1->[0]->{col_str}; # print "Nseq: $Nseq \n"; exit;

   my $Gap = 0;           # gap penalty - negative number is a penalty

   print STDERR "Begin NW initialization.\n";
   # Initialize:
   my @matrix;
   $matrix[0][0][0] = 0;
   $matrix[0][0][1] = 'n';      # upper-left corner - nowhere to go.
   for (my $j=1; $j <= $L1; $j++) {
      $matrix[0][$j][0] = $Gap * $j;
      $matrix[0][$j][1] = 'l';
   }
   for (my $i=1; $i <= $L2; $i++) {
      $matrix[$i][0][0] = $Gap * $i;
      $matrix[$i][0][1] = 'u';
   }

   print STDERR "Begin NW fill.\n";
   # Fill:
   for (my $i=1; $i <= $L2; $i++) {
      my $string2 = $A2->[$i-1]->{col_str};
      my $string2_tilde = $string2;
      $string2_tilde =~ s/-/~/g; # tildes for gaps
      for (my $j=1; $j <= $L1; $j++) {
         my $string1 = $A1->[$j-1]->{col_str};
         my ($diag_score, $left_score, $up_score);
         my ($ccs, $ccinfo) = $self->col_col_score($string1, $string2, $string2_tilde);
         $diag_score =  $matrix[$i-1][$j-1][0] + $ccs;

         $up_score = $matrix[$i-1][$j][0]; # + $Gap; # Since gap is zero, don't do the addition.
         $left_score = $matrix[$i][$j-1][0]; # + $Gap;

         if ($diag_score >= $up_score) {
            if ($diag_score >= $left_score) {
               $matrix[$i][$j][0] = $diag_score;
               $matrix[$i][$j][1] = 'd';
            } else {
               $matrix[$i][$j][0] = $left_score;
               $matrix[$i][$j][1] = 'l';
            }
         } else {
            if ($up_score >= $left_score) {
               $matrix[$i][$j][0] = $up_score;
               $matrix[$i][$j][1] = 'u';
            } else {
               $matrix[$i][$j][0] = $left_score;
               $matrix[$i][$j][1] = 'l';
            }
         }
      }
   }
   print STDERR "Begin NW traceback.\n";
   # Traceback:
   my @align1 = ();
   my @align2 = ();
   my @scores = (); # array of col-col scores for the optimal alignment-alignment alignment.
   my @colcolmatch_infostrings = ();
   #   my @aas = ();
   my $j = $L1;
   my $i = $L2;

   my $gap_col_string = '****************';
   while (length $gap_col_string < $Nseq) {
      $gap_col_string .= $gap_col_string;
   }
   $gap_col_string = substr($gap_col_string, 0, $Nseq);

   while (1) {
      if ($matrix[$i][$j][1] eq 'd') {
         unshift @align1, $A1->[$j-1]->{col_str};
         unshift @align2, $A2->[$i-1]->{col_str};
         my $s2 = $A2->[$i-1]->{col_str};
         my $s2_tg = $s2;
         $s2_tg =~ s/-/~/g;
         #        my $s2_sng = $s2;
         #        $s2_sng =~ s/[^-]/*/g;
         my ($ccs, $ccinfo) = $self->col_col_score($A1->[$j-1]->{col_str}, $s2, $s2_tg);
         unshift @scores,  $ccs; #        $Nseq - hd($A1->[$j-1]->{col_str}, $s2);
         unshift @colcolmatch_infostrings, $ccinfo;
         #       unshift @aas, min($A1->[$j-1]->get_nongap_count('-'), $A2->[$i-1]->get_nongap_count('-')); # min of nongap counts of the two aligned cols
         $i--; $j--;
      } elsif ($matrix[$i][$j][1] eq 'l') {
         unshift @align1, $A1->[$j-1]->{col_str};
         unshift @align2, $gap_col_string;
         unshift @scores, 0;
         unshift @colcolmatch_infostrings, "0 0 0 0";
         #         unshift @aas, 1;
         $j--;
      } elsif ($matrix[$i][$j][1] eq 'u') {
         unshift @align1, $gap_col_string;
         unshift @align2, $A2->[$i-1]->{col_str}; #substr($seq2, $i-1, 1);
         unshift @scores, 0;
         unshift @colcolmatch_infostrings, "0 0 0 0";
         #         unshift @aas, 1;
         $i--;
      } elsif ($matrix[$i][$j][1] eq 'n') {
         last;
      } else {
         die "Unknown NW pointer value: " . $matrix[$i][$j][1], "\n";
      }
   }
   print STDERR "NW traceback done. Returning. \n";
   map(s/~/-/g, @align2);       # restore - as gap character.
   $self->{Acols} = \@align1;
   $self->{Bcols} = \@align2;
   $self->{cc_scores} = \@scores;
   $self->{cc_info} = \@colcolmatch_infostrings;
   #  return (\@align1, \@align2, \@scores, \@colcolmatch_infostrings);
}

sub col_col_score{
   my $self = shift;
   my $str1 = shift;           # string containing an alignment column
   my $str2 = shift;
   my $str2_tildegap = shift;
   my $Nseq = length $str1;
   my $N0 = hd($str1, $str2);   # 
   my $N1 = hd($str1, $str2_tildegap);
   my $ng1 = ($str1 =~ tr/-//);
   my $ng2 = ($str2 =~ tr/-//);
   my $n_AA = $Nseq - $N1;
   my $n_gapgap = $N1 - $N0;
   my $n_Agap = $ng1 + $ng2 - 2*$n_gapgap;
   my $n_AB = $N0 - $n_Agap;
   my $score = $n_AA + $n_gapgap*$self->{gg_score} + $n_Agap*$self->{ag_score} + $n_AB*$self->{ab_score}; 
   my $info_str = "$n_AA $n_gapgap $n_Agap $n_AB";
   return ($score, $info_str);
}

sub alal_cols{
   my $self = shift;
   my $min_frac = shift || 0;
   my $array_of_cols_A = $self->{Acols};
   my $array_of_cols_B = $self->{Bcols};
   my $array_of_scores = $self->{cc_scores};
   my $array_of_ccinfo = $self->{cc_info};
   my $alal_col_string = '';
   my $i = 0;
   for (@$array_of_cols_A) {
      my ($n_AA, $n_gg, $n_Ag, $n_AB) = split(" ", $array_of_ccinfo->[$i]);
      my $N = length $array_of_cols_A->[$i];
      # print STDERR "$N $n_AA $n_gg $n_Ag $n_AB \n";
      my $sum_n = $n_AA + $n_gg + $n_Ag + $n_AB;
      die "!!!! $N $n_AA $n_gg $n_Ag $n_AB \n" if($sum_n > 0  and  $N != $sum_n);
      my $frac = $array_of_scores->[$i]/$N; # ($n_AA > 0)? $n_AA/($n_AA + $n_AB + 0.5*$n_Ag) : 0.0;
      if ($frac >= $min_frac) {
         my $Astring = $array_of_cols_A->[$i];
         my $Bstring = $array_of_cols_B->[$i];
         $alal_col_string .= sprintf("%4i   %s   %4i  %5.3f\n", $i, BagreeA($Astring, $Bstring), $array_of_scores->[$i], $frac);
      }
      $i++;
   }
   return $alal_col_string;
}

sub alal_fasta{
   my $self = shift;
   my $min_frac = shift || 0.0;
   my $interleaved = shift || 0;
   my $array_of_cols_A = $self->{Acols};
   my $array_of_cols_B = $self->{Bcols};
   my $array_of_scores = $self->{cc_scores};
   my $array_of_ccinfo = $self->{cc_info};
   my $Nseq = length $array_of_cols_A->[0];
   my $seq_length = scalar @$array_of_cols_A;
   my $min_score = $min_frac * $Nseq;
   my %id_string = ();
   for my $id (@{$self->{ids}}) {
      $id_string{$id} = '';
   }
   my $fasta_string = '';
   if ($interleaved) {
      for my $col (0 .. $seq_length-1) {
         if ($array_of_scores->[$col] >= $min_score) {
            my @Achars = split('',  $array_of_cols_A->[$col]);
            my @Bchars = split('',  $array_of_cols_B->[$col]);
            for my $row (0 .. $Nseq-1) {
               my $id = $self->{ids}->[$row];
               $id_string{$id} .= $Achars[$row];
               $id_string{$id} .= $Bchars[$row];
            }
         }
      }
   } else {
      for my $col (0 .. $seq_length-1) {
         if ($array_of_scores->[$col] >= $min_score) {
            my @Achars = split('',  $array_of_cols_A->[$col]);
            for my $row (0 .. $Nseq-1) {
               my $id = $self->{ids}->[$row];
               $id_string{$id} .= $Achars[$row];
            }
         }
      }
      for my $col (0 .. $seq_length-1) {
         if ($array_of_scores->[$col] >= $min_score) {
            my @Bchars = split('',  $array_of_cols_B->[$col]);
            for my $row (0 .. $Nseq-1) {
               my $id = $self->{ids}->[$row];
               $id_string{$id} .= $Bchars[$row];
            }
         }
      }
   }

   my @ids = sort keys %id_string;
   # while (my ($id, $str) = each %id_string) {
   for my $id (@ids) {
      my $str = $id_string{$id};
      if ($str ne '') {
         $fasta_string .= sprintf(">%s\n", $id);
         $fasta_string .= sprintf("%s\n", $str);
      }
   }
   return $fasta_string;
}

sub hd{
   return (($_[0] ^ $_[1]) =~ tr/\001-\255//); #  ^ will be non-zero iff chars differ, tr counts non-zeroes.
}

sub BagreeA{
   # construct a string from strings A and B which has same char as A and B where
   # they agree, and X for two distinct aa's, x for aa and gap.
   my $sA = shift;
   my $sB = shift;
   my $BagreeA = '';

   for (my $i=0; $i<length $sA; $i++) {
      my $cA = substr($sA,$i,1);
      my $cB = substr($sB,$i,1);
      my $agree_char;
      if ($cA eq $cB) {
         $agree_char = $cA;
      } elsif ($cA eq '-' or $cB eq '-') {
         $agree_char = 'x';
      } else {
         $agree_char = 'X';
      }
      $BagreeA .= $agree_char;
   }
   return $BagreeA;
}
