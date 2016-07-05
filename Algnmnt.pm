package Algnmnt;
use strict;
use List::Util qw ( min max sum );

my $NO_ID = 'NOT_AN_ACTUAL_ID';

sub  new {
  my $class = shift;
  my $arg = shift; # either filename or string with contents of fasta file
  my $min_nongap_fraction = shift || 0;
  my $min_nongap_chars = shift || 0;
  my $seed = shift || undef;
  my $args= {};
  my $self = bless $args, $class;

  if (defined $seed) {
    srand($seed);
  } else {
    srand();
  }
  my @ids = ();
  my %id_overlapseq = ();
  my %id_sequence = ();
  my @lines = ();
  #  print "arg: $arg\n";
  #  print "arg is file? ", -f $arg, "\n";
  #  print "XXX: ",  $arg =~ /\n/, "\n";
  #  if ((! $arg =~ /\n/) and 
  if (!($arg =~ /\n/) and -f $arg) { # $arg is filename
    open my $fhin, "<$arg";
    @lines = <$fhin>;
    #    print "XXXXXXXlines: ", join("\n", @lines), "\n";
    close $fhin;
	
  } else {			# treat arg as string
    @lines = split("\n", $arg);
  }
  my $id = $NO_ID;
  my $sequence = '';
  while (@lines) {
    my $line = shift @lines;
    #   print $line;
    if ($line =~ /^>/) {
      if ($id ne $NO_ID) {
	$id_sequence{$id} = $sequence;
	push @ids, $id;
      }
      $id = $line;
      $id =~ s/^>\s*//;
      $id =~ s/\s+$//;
      $sequence = '';
    } else {
      $line =~ s/^\s+//;
      $line =~ s/\s+$//;
      $sequence .= $line;
    }
  }
  if (! exists $id_sequence{$id}  and  $sequence ne '') { # take care of the last id-sequence pair.
    $id_sequence{$id} = $sequence;
    push @ids, $id;
  }
# sort the ids
@ids = sort @ids;
  $self->{id_seq} = \%id_sequence;
  $self->{ids} = \@ids;
  $self->{cols} = []; # array of column objects.

  my $n_sequences = scalar @ids;
  my $n_required = ($min_nongap_fraction >= 1)? $n_sequences: int ($min_nongap_fraction * $n_sequences) + 1;
  $n_required = $min_nongap_chars if($min_nongap_chars > $n_required);
  $self->{n_required} = $n_required;
  my $seq_length = length $id_sequence{$ids[0]};
  $self->{align_length} = $seq_length;
  $self->{n_sequences} = $n_sequences;
  my @position_counts = ((0) x $seq_length);
  my @position_aas = (('') x $seq_length);
  my $max_position_counts = -1;
  foreach my $id (@ids) {
    my $sequence = $id_sequence{$id};
    my $seql = length $sequence;
    die "Non-equal sequence lengths in alignment: $seq_length, $seql. id: $id \nsequence: $sequence\n" if($seql ne $seq_length);

    for (my $i = 0; $i < $seq_length; $i++) {
      my $aa = substr($sequence, $i, 1);
      if ($aa ne '-') {
	$position_counts[$i]++;
	#print "ZZZZZZZZ: $id $i $aa ", $position_counts[$i], "\n";
	$max_position_counts = $position_counts[$i] if($position_counts[$i] >= $max_position_counts);
	if (!($position_aas[$i] =~ /$aa/)) {
	  $position_aas[$i] .= $aa; # not invariant
	}
      }
    }
  }


#print "$min_nongap_fraction $min_nongap_chars  $n_required  $n_sequences \n";# exit;
  for (my $i=0; $i < $seq_length; $i++){
     my $col_obj = Column->new($self->{id_seq}, $self->{ids}, $i);
#     print STDERR $col_obj->{col_str}, "\n";
#print "nongap count:  ", $col_obj->get_nongap_count('-'),"   nongaps required: ", $n_required, "\n";
     push $self->{cols}, $col_obj if($col_obj->get_nongap_count('-') >= $n_required);
  }

  my $n_invariant = 0;
  foreach (@position_aas) {
    # print "[$_]\n";
    $n_invariant++ if(length $_ == 1);
  }
  #  print "n_invariant: $n_invariant  align length: ", scalar @position_aas, "\n"; 
  #  print "pinv: ", $n_invariant/$seq_length, "\n";
  $self->{position_counts} = \@position_counts;
 
  my $overlap_length = 0;
  my %id_overlapnongapcount = ();
  my $overlap_n_invariant = 0;
  foreach my $id (@ids) {
    $id_overlapseq{$id} = ''; $id_overlapnongapcount{$id} = 0;
  }
  foreach my $position (0..@position_counts-1) {
    my $count = $position_counts[$position];
    #   print "XXXX: $count $n_required [", ($count >= $n_required)? '1': '0', "]\n";
    #    exit;
    if ($count >= $n_required) {
      $overlap_length++;
      foreach my $id (@ids) {
	my $char = substr($id_sequence{$id}, $position, 1);
	$id_overlapseq{$id} .= $char;
	$id_overlapnongapcount{$id}++ if($char ne '-');	
      }
      $overlap_n_invariant++ if(length $position_aas[$position] == 1);
    }
  }
  #  print "overlap n_invariant: $overlap_n_invariant,  length: $overlap_length\n";
  #  print "overlap pinv: ", $overlap_n_invariant/$overlap_length, "\n";

  $self->{id_overlapseq} = \%id_overlapseq;
  $self->{id_overlapnongapcount} = \%id_overlapnongapcount;
  # die "overlap length inconsistency??? $overlap_length \n" if($overlap_length > 0 and $overlap_length != length $id_overlapseq{$ids[0]});

  if ($overlap_length > 0 and $overlap_length != length $id_overlapseq{$ids[0]}) {
    for (@ids) {
      print $_, "  ", length $id_overlapseq{$_}, "\n";
    }
    die "$overlap_length  ", length $id_overlapseq{$ids[0]}, "\n";
  }
  $self->{overlap_length} = $overlap_length;
  #	$self->{ids} = \@ids;

  return $self;
}


sub  weed_sequences{
  # weed out sequences which have poor overlap with others
  my $self = shift;
  my $fraction = shift || 0.3;
  my $min_nongapcount = int($fraction * $self->{overlap_length});
  #my %id_overlap_count = (); # 
  my @ids = $self->{id_overlapnongapcount};
  foreach (@ids) {
    if ($self->{id_overlapnongapcount}->{$_} < $min_nongapcount) {
      # delete this sequence
      delete $self->{id_overlapseq}->{$_};
      delete $self->{id_overlapnongapcount}
    }
  }
}


sub align_fasta_string{
  my $self = shift;
  my $spacer = shift || '';
  my $align_fasta = '';
  foreach my $id (@{$self->{ids}}) {
    my $sequence = $self->{id_seq}->{$id};
    $align_fasta .= ">$spacer$id\n$sequence\n";
  }
  chomp $align_fasta;
  return $align_fasta;
}

sub overlap_fasta_string{
  my $self = shift;
  my $spacer = shift || '';
  my $overlap_fasta = '';
  foreach my $id (@{$self->{ids}}) {
    my $sequence = $self->{id_overlapseq}->{$id};
    $overlap_fasta .= ">$spacer$id\n$sequence\n";
  }
  return $overlap_fasta;
}

sub get_overlap_hashref{
  my $self = shift;
  return $self->{id_overlapseq};
}

sub interleaved_overlap_fasta_string{
  my $self= shift;
  my $max_length = shift || 1000;
  my $overlap_fasta = '';
  foreach my $id (@{$self->{ids}}) {
    $overlap_fasta .= "$id \n";
  }
  foreach my $id (@{$self->{ids}}) {
    my $sequence = $self->{id_overlapseq}->{$id};
    $sequence = substr($sequence, 0, 120);
    $overlap_fasta .= "$sequence\n";
  }
  chomp $overlap_fasta;
  return $overlap_fasta;
}


sub overlap_nexus_string{ # basic nexus format string for use by MrBayes.
  my $self = shift;
  my $n_leaves = scalar @{$self->{ids}}; 
  my $overlap_length = length ($self->{id_overlapseq}->{$self->{ids}->[0]});
  my $nexus_string = "#NEXUS\n" . "begin data;\n";
  $nexus_string .= "dimensions ntax=$n_leaves nchar=$overlap_length;\n";
  $nexus_string .= "format datatype=protein interleave=no gap=-;\n";
  $nexus_string .= "matrix\n";

  foreach my $id (@{$self->{ids}}) {
    my $sequence = $self->{id_overlapseq}->{$id};
    $id =~ s/[|].*//;
    my $id_padded =  $id . "                                                  ";
    my $target_l = 50;
    my $l = length $id;
    while ($target_l < $l+5) {
      $target_l += 10;
    }
    $id_padded = substr($id_padded, 0, $target_l);
    $nexus_string .= "$id_padded$sequence\n";
  }
  $nexus_string .= "\n;\n\n" . "end;\n";
  return $nexus_string;
}

sub bootstrap_overlap_fasta_string{
  my $self = shift;
  my $spacer = shift || '';
  my %id_bootstrapoverlapseq = ();
  my $overlap_length = $self->{overlap_length};

  my @indices = ();
  for (1..$overlap_length) {
    my $index = int( rand($overlap_length) );
    push @indices, $index;
    #	$index_count{$index}++;
  }

  for my $id (@{$self->{ids}}) {
    my $std_overlap = $self->{id_overlapseq}->{$id};
    my $string = '';
    foreach my $index (@indices) {
      $string .=  substr($std_overlap, $index, 1);
    }
    $id_bootstrapoverlapseq{$id} = $string;
  }

  my $bofstring = '';
  foreach my $id (@{$self->{ids}}) {
    #	print "id, seq: $id; $sequence\n";
    my $sequence = $id_bootstrapoverlapseq{$id};
    $bofstring .= ">$spacer$id\n$sequence\n";
  }
  chomp $bofstring;
  return $bofstring;
}

sub get_n_sequences{
  my $self = shift;
  return scalar @{$self->{ids}};
}

sub get_overlap_length{
  my $self = shift;
  return $self->{overlap_length};
}
sub set_overlap_length{
  my $self = shift;
  $self->{overlap_length} = shift;
}

sub array_of_cols{
   my $self = shift;
   my @sorted_ids = sort keys %{$self->{id_seq}};
   my $align_length = $self->{align_length};
   my @array_of_col_arrays = ();
   my @array_of_nongap_counts = ();
for(my $i=0; $i < $align_length; $i++){ # initialize to array of $align_length empty array refs.
   push @array_of_col_arrays, [];
}

   for my $id (@sorted_ids){
      my $the_sequence = $self->{id_seq}->{$id};
      my @row_chars = split("", $the_sequence);
      for (my $i=0; $i < $align_length; $i++){ # push each char in row (gapped seq) onto appropriate col
         my $the_char = $row_chars[$i];
         push @{$array_of_col_arrays[$i]}, $the_char;
         $array_of_nongap_counts[$i]++ if($the_char ne '-');
      }
   }
   return \@array_of_col_arrays, \@array_of_nongap_counts;

}

sub array_of_col_objects{
   my $self = shift;
#   print STDERR "XXX: ", $self->{cols}->[0]->{col_str}, "\n";
   return $self->{cols};
}


sub array_of_overlap_cols{
   my $self = shift;
   my @sorted_ids = sort keys %{$self->{id_overlapseq}};
   my $overlap_length = $self->{overlap_length};
   my @array_of_col_arrays = ();
   my @array_of_col_strings = ();
   my @array_of_nongap_counts = ();
my @array_of_plurality_counts = ();
for(my $i=0; $i < $overlap_length; $i++){ # initialize to array of $align_length empty array refs.
   push @array_of_col_arrays, [];
}
 my @col__AA_count = ();
   for (my $i=0; $i < $overlap_length; $i++){ $col__AA_count[$i] = {}; }
   for my $id (@sorted_ids){
      my $the_sequence = $self->{id_overlapseq}->{$id};
      my @row_chars = split("", $the_sequence);
      for (my $i=0; $i < $overlap_length; $i++){ # push each char in row (gapped seq) onto appropriate col
         my $the_char = $row_chars[$i];
         push @{$array_of_col_arrays[$i]}, $the_char;
         $array_of_col_strings[$i] .= $the_char;
         $array_of_nongap_counts[$i]++ if($the_char ne '-');
         $col__AA_count[$i]->{$the_char}++ if($the_char ne '-');
      }
   }
   for my $AA_count (@col__AA_count){
      my $max_AA_count = -1;
     for my $AA (keys %$AA_count){
        if($AA_count->{$AA} > $max_AA_count){
           $max_AA_count = $AA_count->{$AA};
        }
     }
      push @array_of_plurality_counts, $max_AA_count;
   }
   return (\@array_of_col_strings, \@array_of_col_arrays, \@array_of_nongap_counts, \@array_of_plurality_counts);

}

1;
