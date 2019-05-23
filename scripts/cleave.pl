#!/usr/bin/perl -w

use strict;
use FindBin;

my $home = "$FindBin::Bin";

if($#ARGV < 0) {
  print "Usage: cleave.pl <fasta_file1> <fasta_file2> ...\n";
  exit;
}

for(my $i=0; $i<$#ARGV+1; $i++) {
  my $cmd = "$home/predict_multiple_cathepsin_cleavage.pl $ARGV[$i] $home/../libs/cleavage/catS_15_count.txt $home/../libs/cleavage/catS_60_count.txt $home/../libs/cleavage/catS_240_count.txt $home/../libs/cleavage/catS_choe.txt $home/../libs/cleavage/cathepsinS_matrix.txt  $home/../libs/cleavage/catB_15_count.txt $home/../libs/cleavage/catB_60_count.txt  $home/../libs/cleavage/catB_240_count.txt $home/../libs/cleavage/cathepsinB_matrix.txt";
  print "$cmd\n";
  `$home/predict_multiple_cathepsin_cleavage.pl $ARGV[$i] $home/../libs/cleavage/catS_240_count.txt $home/../libs/cleavage/catS_240_count.txt $home/../libs/cleavage/catS_240_count.txt $home/../libs/cleavage/catS_choe.txt $home/../libs/cleavage/cathepsinS_matrix.txt  $home/../libs/cleavage/catB_15_count.txt $home/../libs/cleavage/catB_60_count.txt  $home/../libs/cleavage/catB_240_count.txt $home/../libs/cleavage/cathepsinB_matrix.txt > tmp_file`;
  $cmd = "$home/predict_multiple_cathepsin_cleavage.pl tmp_file $home/../libs/cleavage/catH_15_count.txt $home/../libs/cleavage/catH_60_count.txt $home/../libs/cleavage/catH_240_count.txt --cleavage_threshold 2 -nterm 1";
    print "$cmd\n";

  `$home/predict_multiple_cathepsin_cleavage.pl tmp_file $home/../libs/cleavage/catH_15_count.txt $home/../libs/cleavage/catH_60_count.txt $home/../libs/cleavage/catH_240_count.txt --cleavage_threshold 2 -nterm 1 > $ARGV[$i].out`
}
