#!/usr/bin/perl -w

package UTIL;

use strict;
use File::Basename;

sub read_fasta_sequences {
  my $input_file = shift;
  my @sequences = ();
  open(IN, $input_file) or die "Cannot open $input_file\n";
  my $curr_sequence = '';
  while(my $line = <IN>) {
    if($line =~ /^>/) {
      #print "Header = $line\n";
      if(length($curr_sequence) > 0) {
        push(@sequences, $curr_sequence);
        $curr_sequence = '';
      }
    } else {
      chomp($line);
      $line =~ s/\W//g;
      $curr_sequence .= $line;
    }
  }
  push(@sequences, $curr_sequence);
  close(IN);
  return @sequences;
}

sub trim_extension {
  my $str = shift;
  $str =~ s/\.[^.]+$//;
  return $str;
}

sub pMHC_name_signature {
  my $str = shift;
  my $dirname = basename(dirname($str));
  my $signature = $dirname. "_" . UTIL::trim_extension(basename($str));
  if($dirname eq ".") { $signature = UTIL::trim_extension(basename($str)); }
  return $signature;
}

1;
