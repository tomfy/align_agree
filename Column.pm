package Column;
use strict;
use List::Util qw ( min max sum );

use constant F => 0.7;

sub  new {
   my $class = shift;
   my $id_seq = shift;      # hashref; keys: sequence ids, values: aligned sequences (with gaps)
   my $ids = shift;         # array ref of sorted ids
   my $position = shift;    # index of column
   #my $id_chars = shift;
   my $args= {};
   my $self = bless $args, $class;
   #  $self->{col_str} = '';

   my $str = '';
   my %aa_count = ('-' => 0, '~' => 0); # counts aa's and also any other chars such as - or ~
   my $qual_sum = 0;
   for my $id (@$ids) {
      my $length = length $id_seq->{$id};
      my $char = substr($id_seq->{$id}, $position, 1);
      $str .= $char;
      $aa_count{$char}++;
      my $min_dist_to_gap = 5;
      for my $d (0..4) {
         if ($position-$d >= 0) {
            $min_dist_to_gap = $d if(substr($id_seq->{$id}, $position-$d, 1) eq '-');
         }
         if ($position+$d < $length) {
            $min_dist_to_gap = $d if(substr($id_seq->{$id}, $position+$d, 1) eq '-');
         }
         last if($min_dist_to_gap == $d);
      }
      $qual_sum += (1.0 - F**$min_dist_to_gap)

   }
   $self->{position} = $position;
   $self->{qual_sum} = $qual_sum;
   $self->{col_str} = $str;
   $self->{aa_count} = \%aa_count;

   return $self;
}

sub get_col_string{
   my $self = shift;
   return $self->{col_str};
}

sub get_char_count{ # how many of a particular character are there in the column
   my $self = shift;
   my $char = shift;
   #  print "char: $char \n";
   return (exists $self->{aa_count}->{$char})? $self->{aa_count}->{$char} : 0;
}

sub get_col_length{
   my $self = shift;
   return length $self->{col_str};
}

sub get_gap_count{
   my $self = shift;
   my $gap_char = shift;
   #   print "gap char:  $gap_char \n";
   return  $self->get_char_count($gap_char);
}

sub get_nongap_count{
   my $self = shift;
   my $gap_char = shift;
   #  print "gap char:  $gap_char \n";
   return $self->get_col_length() - $self->get_gap_count($gap_char);
}

sub get_qual_sum{
   my $self = shift;
   return (exists $self->{qual_sum})? $self->{qual_sum} : 0;
}

sub get_position{
   my $self = shift;
   return $self->{position};
}

1;
