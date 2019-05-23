#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;

my $home = "$FindBin::Bin";
my $modeller_home = "/salilab/diva1/home/modeller";

if ($#ARGV != 0) {
  print "Usage: iterate_model_pMHC.pl <TCRs_fasta_file>\n";
  exit;
}

my $sequence_file = $ARGV[0];
my @sequences = read_fasta_sequences($sequence_file);

my $seq_counter = 1;
for(my $i=0; $i<@sequences-1; $i+=2, $seq_counter++) {
  print $i;
  my $header = $sequences[$i];
  my $sequence = $sequences[$i+1];
  my $dirname = "seq_$seq_counter";
  print "$dirname $header $sequence\n";
  mkdir $dirname; # or die "Can't make directory $dirname\n";
  chdir $dirname;

  open OUT, ">seq";
  print OUT "$header$sequence\n";
  close OUT;

  open OUT2, ">TCR.ali";
  print OUT2 ">P1;TCR\n";
  print OUT2 "sequence:TCR:::::::0.00: 0.00\n";
  print OUT2 "$sequence*\n";
  close OUT2;

  my $cmd = "$modeller_home/modpy.sh python $home/modelTCR.py seq >& model.log";
  print "$cmd\n";
  `$cmd`;
  chdir "..";

}

sub read_fasta_sequences {
  my $input_file = shift;
  my @sequences = ();
  open(IN, $input_file) or die "Cannot open $input_file\n";
  my $curr_sequence = '';
  my $curr_header;
  while(my $line = <IN>) {
    if($line =~ /^>/) {
      $curr_header = $line;
      #print "Header = $line\n";
      if(length($curr_sequence) > 0) {
        push(@sequences, ($curr_header, $curr_sequence));
        $curr_sequence = '';
      }
    } else {
      chomp($line);
      $curr_sequence .= $line;
    }
  }
  if(length($curr_sequence) > 0) {
    push(@sequences, ($curr_header, $curr_sequence));
  }
  close(IN);
  return @sequences;
}
