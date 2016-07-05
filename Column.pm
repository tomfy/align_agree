package Column;
use strict;
use List::Util qw ( min max sum );

sub  new {
   my $class = shift;
   my $id_seq = shift; # hashref; keys: sequence ids, values: aligned sequences (with gaps)
   my $ids = shift;    # array ref of sorted ids
   my $position = shift;        # index of column
   my $args= {};
   my $self = bless $args, $class;
   $self->{col_str} = '';
   $self->{aa_count} = {'-' => 0, '~' => 0}; 

   for my $id (@$ids) {
      my $char = substr($id_seq->{$id}, $position, 1);
      $self->{col_str} .= $char;
      $self->{aa_count}->{$char}++;
   }

   return $self;
}

sub get_col_string{
my $self = shift;
return $self->{col_str};
}

sub get_char_count{
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

1;
