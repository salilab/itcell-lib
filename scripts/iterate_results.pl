#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;

my $home = "$FindBin::Bin";
require "$home/UTIL.pm";

if ($#ARGV != 1 and $#ARGV != 2 and $#ARGV != 3) {
  print "Usage: iterate_results.pl <peptide_sequence_file> <pMHC_template> [cleaved_sequence_file] [pMHC score file]\n";
  print "Assumes that template MHC chains are AB, and peptide chain is C\n";
  exit;
}

my $sequence_file = $ARGV[0];
my $pMHC = $ARGV[1];
my $pMHC_signature = UTIL::pMHC_name_signature($pMHC);

my $cleaved_sequence_file = $sequence_file;
if($#ARGV == 2) {
  $cleaved_sequence_file = $ARGV[2];
}

my $soap_file =  UTIL::trim_extension($sequence_file) . "_" . $pMHC_signature . ".soap_pep";
if($#ARGV == 3) {
  $soap_file = $ARGV[3];
}

my @sequences = UTIL::read_fasta_sequences($sequence_file);
my @cleaved_sequences =  UTIL::read_fasta_sequences($cleaved_sequence_file);
my @scores = read_soap_file($soap_file);

my $peptide = `$home/pdb2fasta $pMHC | tail -n1`;
chomp $peptide;
my $peptide_length = length($peptide);

my $average_pMHC_score = 0.0;
my $std_pMHC_score = 0.0;
my $average_TCR_score = 0.0;
my $std_TCR_score = 0.0;
my $counter = 0;

my $window_size = $peptide_length;

my $filtered_file_name = UTIL::trim_extension($sequence_file) . "_" . $pMHC_signature . "_filtered.txt";
open OUT, ">$filtered_file_name";
foreach my $sequence (@sequences) {
  for(my $i=0; $i<length($sequence)-($window_size-1); $i++) {
    my $peptide_sequence = substr($sequence, $i, $window_size);
    my $name_signature = $pMHC_signature . "_" . $peptide_sequence;
    my $TCR_score_file = $name_signature . "_TCR_score_SOAP.txt";

    #print "$TCR_score_file\n";

    if(-e $TCR_score_file) {
      my $pMHC_score = $scores[$counter];
      chomp $pMHC_score;  $pMHC_score += 0.0;
      my $TCR_score = `tail -n1 $TCR_score_file | cut -d '|' -f2`;
      chomp $TCR_score;  $TCR_score += 0.0;
      my $total_score = $pMHC_score + $TCR_score;

      $average_pMHC_score += $pMHC_score;
      $std_pMHC_score += ($pMHC_score*$pMHC_score);
      $average_TCR_score += $TCR_score;
      $std_TCR_score += ($TCR_score*$TCR_score);
      $counter++;
      #print "$peptide_sequence | $total_score | $pMHC_score | $TCR_score\n";
    }
  }
}

$average_pMHC_score /= $counter;
$average_TCR_score /= $counter;

$std_pMHC_score /= $counter;
$std_pMHC_score -= ($average_pMHC_score*$average_pMHC_score);
$std_pMHC_score = sqrt($std_pMHC_score);

$std_TCR_score /= $counter;
$std_TCR_score -= ($average_TCR_score*$average_TCR_score);
$std_TCR_score = sqrt($std_TCR_score);

#print "$average_pMHC_score $std_pMHC_score\n";

$counter = 0;
# conversion to Z-scores
foreach my $sequence (@sequences) {
  for(my $i=0; $i<length($sequence)-($window_size-1); $i++) {
    my $peptide_sequence = substr($sequence, $i, $window_size);
    my $name_signature = $pMHC_signature . "_" . $peptide_sequence;
    my $TCR_score_file = $name_signature . "_TCR_score_SOAP.txt";

    if(-e $TCR_score_file) {
      my $pMHC_score = $scores[$counter];
      chomp $pMHC_score;  $pMHC_score += 0.0;
      my $TCR_score = `tail -n1 $TCR_score_file | cut -d '|' -f2`;
      chomp $TCR_score;  $TCR_score += 0.0;
      my $total_score = $pMHC_score + $TCR_score;
      # normalize
      my $pMHC_Zscore  = sprintf("%.2f", ($pMHC_score - $average_pMHC_score)/$std_pMHC_score);
      my $TCR_Zscore = sprintf("%.2f", ($TCR_score - $average_TCR_score)/$std_TCR_score);
      my $Zscore = sprintf("%.2f", $pMHC_Zscore + $TCR_Zscore);
      $counter++;
      print "$counter | $peptide_sequence | $Zscore | $pMHC_Zscore | $TCR_Zscore | $total_score | $pMHC_score | $TCR_score \n";

      if($pMHC_Zscore < -0.5 and $TCR_Zscore < -0.5) {
        # check if cleaved
        my $cleaved = 1;
        foreach my $seq (@cleaved_sequences) {
          if(index($seq, $peptide_sequence) != -1) {
            #print $peptide_sequence . " in " . $seq . "\n";
            $cleaved = 0;
          }
        }

        if($cleaved == 0) {
          print OUT "$counter | $peptide_sequence | $Zscore | $pMHC_Zscore | $TCR_Zscore | $total_score | $pMHC_score | $TCR_score\n";
        }
      }
    }
  }
}


sub read_soap_file {
  my $input_file = shift;
  my @scores = ();
  open(IN, $input_file) or die "Cannot open $input_file\n";
  while(my $line = <IN>) {
    if($line =~ /^filename/) {
      #print "Header = $line\n";
      # read the next line too
      $line = <IN>;
      #print "Header2 = $line\n";
    } else {
      chomp($line);
      my @tmp = split('\|', $line);
      my $score = $tmp[1];
      #print $score;
      push (@scores, $score);
    }
  }
  close(IN);
  return @scores;
}
