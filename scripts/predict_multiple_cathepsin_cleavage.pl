#!/usr/bin/perl -w

use strict;
use Getopt::Long;

# File: predict_multiple_cathepsin_cleavage.pl
# Created by: Natalia Khuri
# Date: May 27, 2013

if($#ARGV < 1) {
    print "Usage:\n";
    print "\tperl predict_multiple_cathepsin_cleavage.pl fasta_file matrix_file1 matrix_file2 ...\n";
    print "\tSee for example matrix.txt and p01308.txt\n";
    print "--cleavage_threshold = cleavage threshold\n";
    exit(0);
}

my $cleavage_threshold = 3.0;
my $nterm = 0;
my $cterm = 0;

my $input_file = $ARGV[0]; # for example, p01308.fasta;

my $result = GetOptions("cleavage_threshold=f" => \$cleavage_threshold,
                        "nterm=i" => \$nterm,
                        "cterm=i" => \$cterm);

my @sequences = read_fasta_sequences($input_file);
# pad sequences:
for(my $i=0; $i <= $#sequences; $i++) {
  my $pad_seq = "XXX" . $sequences[$i] . "XXX";
  $sequences[$i] = $pad_seq;
}

my %aa_to_index = (
        'A' => 0,
        'R' => 1,
        'N' => 2,
        'D' => 3,
        'C' => 4,
        'Q' => 5,
        'E' => 6,
        'G' => 7,
        'H' => 8,
        'I' => 9,
        'L' => 10,
        'K' => 11,
        'M' => 12,
        'F' => 13,
        'P' => 14,
        'S' => 15,
        'T' => 16,
        'W' => 17,
        'Y' => 18,
        'V' => 19
);

my @index_to_aa = ();
foreach my $key ( keys %aa_to_index) {
	$index_to_aa[$aa_to_index{$key}] = $key;
}

# human proteome frequences are from http://www.pasteur.fr/~tekaia/aafreq.html
# note that the total does not add to 1 (.9986)!
my %bg_frequency = (
	'A' => .0697,
	'R' => .0568,
	'N' => .0358,
	'D' => .0468,
	'C' => .0231,
	'Q' => .0471,
	'E' => .0695,
	'G' => .0670,
	'H' => .0264,
	'I' => .0435,
	'L' => .0987,
	'K' => .0572,
	'M' => .0220,
	'F' => .0364,
	'P' => .0637,
	'S' => .0825,
	'T' => .0534,
	'W' => .0127,
	'Y' => .0263,
	'V' => .0600
);

my $sum = 0;
foreach my $aa (keys %bg_frequency)
{
    $sum += $bg_frequency{$aa};
}

my @tables = ();

for(my $i=1; $i<$#ARGV+1; $i++) {
    my $matrix_file = $ARGV[$i]; # for example, matrix.txt;
    my @table = read_frequency_matrix_file($matrix_file);
    push(@tables, \@table);
}


my $window = 8;
my @out_sequences;

# predict cleavage sites
foreach my $seq (@sequences) {
  # skip short sequences
  if(length($seq) < 8) {
      my $subsequence = substr($seq, 3, -3);
      push(@out_sequences, $subsequence);
      next;
  }

  # determine search range
  my $out_sequence = '';
  my $start = 0;
  my $end = length($seq) - $window;
  if($nterm > 0) { $end = $nterm-1; }
  if($cterm > 0) {
     $start = $end - $cterm + 1;
     $out_sequence = substr($seq, 3, $start);
  }

  for (my $i = $start; $i <= $end; $i++) {
    my $subsequence = substr($seq, $i, $window);
    $out_sequence .= substr($seq, $i+3, 1);
    my $cleaved=0;
    for(my $j=1; $j<$#ARGV+1; $j++) {
      my $score = score_sequence_window($subsequence, $tables[$j-1]);
      if($score > $cleavage_threshold) {
        if($cleaved == 0) {
          #print "\n>";
          $cleaved++;
          if($nterm > 0) { $end++; }
          push(@out_sequences, $out_sequence);
          $out_sequence = '';
        }
      }
    }
  }
  my $end_seq = substr($seq, length($seq) - $window/2, 1);
  if($nterm > 0) {
      $end_seq = substr($seq, $end+4, -3);
  }
  $out_sequence .= $end_seq;
  push(@out_sequences, $out_sequence);
}

#print "out_sequences: $#out_sequences+1 \n";
my $total = 0;
my $atotal = 0;
foreach my $seq (@out_sequences) {
    my $seq_length = length($seq);
    $atotal += $seq_length;
    if($seq_length >= 12) {
      $total += ($seq_length - 11);
    }
    print ">\n$seq\n";
}
$atotal -=11;
#print ">$total / $atotal\n";

sub log2 {
  my $n = shift;
  return log($n)/log(2);
}

sub read_frequency_matrix_file {
  my $matrix_file = shift;
  my @table = ();

  my @counts = (0, 0, 0, 0, 0, 0, 0, 0);

  open(MX, $matrix_file);
  while(my $line = <MX>) {
    chomp($line);
    my @record = split("\t", $line);
    my $aa = shift(@record);
    $aa =~ s/\s+//g;
    if(length($aa) > 1) {
      print "Only one-letter AA code is accepted\n";
      exit(0);
    }
    for(my $i = 0; $i < 8; $i++) {
      # print "$record[$i] ";
      if($record[$i] ne 'X') {
        $counts[$i] += $record[$i];
      }
      $table[$aa_to_index{$aa}][$i] = $record[$i];
    }
  }
  close(MX);

  my @sorted = sort {$a <=> $b} @counts;

  #do Laplacian correction
  my $nsamples = pop(@sorted) + 20;
  for(my $i = 0; $i < 20; $i++) {
    for(my $j = 0; $j < 8; $j++) {
      if($table[$i][$j] ne 'X') {
        $table[$i][$j] = ($table[$i][$j] + 1) / $nsamples;
        $table[$i][$j] = log2($table[$i][$j] / $bg_frequency{$index_to_aa[$i]});
      } else {
        $table[$i][$j] = 0;
      }
      #print "$table[$i][$j] ";
    }
    #print " $index_to_aa[$i]\n";
  }
  return @table;
}

sub read_fasta_sequence {
  my $input_file = shift;
  my $seq = '';
  open(IN, $input_file) or die "Cannot open $input_file\n";
  while(my $line = <IN>) {
    if($line =~ /^>/) {
      #print "Header = $line\n";
    } else {
      chomp($line);
      $line =~ s/\W//g;
      $seq .= $line;
    }
  }
  close(IN);
  return $seq;
}

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

sub score_sequence_window {
  my $subsequence = shift;
  my @table = @{(shift)};
  #print @table;
  my $score = 0;
  my @letters = split('', $subsequence);
  for (my $j = 0; $j < scalar(@letters); $j++) {
    if($letters[$j] ne 'X') {
      my $k = $aa_to_index{$letters[$j]};
      #print $letters[$j];
      $score += $table[$k][$j];
    }
  }
  return $score;
}
