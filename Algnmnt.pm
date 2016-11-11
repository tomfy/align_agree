package Algnmnt;
use strict;
use List::Util qw ( min max sum );

my $NO_ID = 'NOT_AN_ACTUAL_ID';

sub  new {
   my $class = shift;
   my $fasta_str_or_filename = shift; # either filename or string with contents of fasta file
   my $args = shift; # e.g. {min_nongap_fraction => 0.3, seed => 134573} 
   my $default_args = {
                       min_nongap_fraction => 0,
                       min_nongap_chars => 0,
                       seed => undef
                      };
   my $self = bless $default_args, $class;

   #  print "Algnmnt object fields: \n";
   for my $option ( keys %$args ) {
      warn "Unknown option: $option in Algnmnt constructor.\n" if ( !exists $self->{$option} );
      if ( defined $args->{$option} ) { # if arg is undef, leaves default in effect
         $self->{$option} = $args->{$option};
         #        print "$option  ", $self->{$option}, "\n";
      }
   }
   #  print "************************************************************\n";

   if (defined $self->{seed}) {
      srand($self->{seed});
   } else {
      srand();
   }
   my @ids = ();
   my %id_overlapseq = ();
   my %id_sequence = ();
   my %id_chars = ();
   my @lines = ();
   if (!($fasta_str_or_filename =~ /\n/) and -f $fasta_str_or_filename) { # $arg is filename
      open my $fhin, "<$fasta_str_or_filename";
      @lines = <$fhin>;
      close $fhin;
   } else {			# treat arg as string
      @lines = split("\n", $fasta_str_or_filename);
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
   $self->{ids} = \@ids;        # sorted ids
   $self->{cols} = [];          # array of column objects.

   my $n_sequences = scalar @ids;
   my $n_required = ($self->{min_nongap_fraction} >= 1)? $n_sequences: int ($self->{min_nongap_fraction} * $n_sequences) + 1;
   $n_required = $self->{min_nongap_chars} if($self->{min_nongap_chars} > $n_required);
   $self->{n_required} = $n_required;

   my $seq_length = length $id_sequence{$ids[0]};
   for (my $i=0; $i < $seq_length; $i++) {
      my $col_obj = Column->new($self->{id_seq}, $self->{ids}, $i); # , $self->{id_chars});
      push $self->{cols}, $col_obj if($col_obj->get_nongap_count('-') >= $n_required);
   }

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
   for (my $i=0; $i < $align_length; $i++) { # initialize to array of $align_length empty array refs.
      push @array_of_col_arrays, [];
   }

   for my $id (@sorted_ids) {
      my $the_sequence = $self->{id_seq}->{$id};
      my @row_chars = split("", $the_sequence);
      for (my $i=0; $i < $align_length; $i++) { # push each char in row (gapped seq) onto appropriate col
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
   for (my $i=0; $i < $overlap_length; $i++) { # initialize to array of $align_length empty array refs.
      push @array_of_col_arrays, [];
   }
   my @col__AA_count = ();
   for (my $i=0; $i < $overlap_length; $i++) {
      $col__AA_count[$i] = {};
   }
   for my $id (@sorted_ids) {
      my $the_sequence = $self->{id_overlapseq}->{$id};
      my @row_chars = split("", $the_sequence);
      for (my $i=0; $i < $overlap_length; $i++) { # push each char in row (gapped seq) onto appropriate col
         my $the_char = $row_chars[$i];
         push @{$array_of_col_arrays[$i]}, $the_char;
         $array_of_col_strings[$i] .= $the_char;
         $array_of_nongap_counts[$i]++ if($the_char ne '-');
         $col__AA_count[$i]->{$the_char}++ if($the_char ne '-');
      }
   }
   for my $AA_count (@col__AA_count) {
      my $max_AA_count = -1;
      for my $AA (keys %$AA_count) {
         if ($AA_count->{$AA} > $max_AA_count) {
            $max_AA_count = $AA_count->{$AA};
         }
      }
      push @array_of_plurality_counts, $max_AA_count;
   }
   return (\@array_of_col_strings, \@array_of_col_arrays, \@array_of_nongap_counts, \@array_of_plurality_counts);

}

1;
