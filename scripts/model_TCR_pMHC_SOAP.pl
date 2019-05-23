#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;
use File::Copy;

my $home = "$FindBin::Bin";
require "$home/UTIL.pm";
require "$home/Transformation.pm";

if ($#ARGV != 2) {
  print "Usage: model_TCR_pMHC.pl <complex_template_DB> <TCR_model_list> <pMHC_pdb>\n";
  print "Assumes that in pMHC_pdb, MHC chains are AB, and peptide chain is C\n";
  exit;
}

my $imp_home = "~/imp_fast/";

my $template_DB = $ARGV[0];
my $TCR_DB = $ARGV[1];
my $pMHC = $ARGV[2];
my $pMHC_aligned = 1; # TODO: add as an option

# remove unnecessary parts from the filename and copy it
my $pMHC_signature = UTIL::pMHC_name_signature($pMHC);
print "pMHC = $pMHC $pMHC_signature\n";

open FILE, $template_DB or die "Can't open file template DB file $template_DB\n";
my @template_files = <FILE>;
close FILE;
open FILE, $TCR_DB or die "Can't open file TCR models file $TCR_DB\n";
my @tcr_files = <FILE>;
close FILE;

my $best_score = 10000.0;
my $best_structure = '';



#open SCORES_OUT, ">$pMHC_signature\_TCR_score_best_scores_SOAP.txt";
#my $scores_files = "$pMHC_signature\_TCR_score_all_scores_SOAP.txt";
#`rm -f $scores_files`;

my $alignment_dirname = $pMHC_signature .".tmp";
mkdir $alignment_dirname;
`cp $home/params.txt $alignment_dirname/`;

# iterate TCRs
foreach my $tcr (@tcr_files) {
  chomp $tcr;
  print "TCR = $tcr\n";
  # copy a
  my $tcr_name = basename(dirname($tcr)) . "_" . basename($tcr);
  if(basename(dirname($tcr)) eq ".") { $tcr_name = basename($tcr) }
  else {
    copy("$tcr", "$tcr_name") or die "Copy failed from $tcr to $tcr_name\n";
  }

  # iterate templates and prepare transformation file for SOAP
  my $trans_file = $pMHC_signature.".trans";
  open OUT, ">$trans_file" or die "Can't open file $trans_file\n";
  my $counter = 0;
  foreach my $template (@template_files) {
    chomp $template;
    print "$counter Template = $template\n";

    chdir $alignment_dirname;

    if ($pMHC_aligned == 1) { # prealigned pMHC and complex template
      # compute T = TCR transformation to the template
      my $template_TCR = $template . "/TCR.pdb";
      my ($rx, $ry, $rz, $tx, $ty, $tz) = get_alignment_transformation("$template_TCR", "../$tcr_name");
      print OUT "$counter $rx $ry $rz $tx $ty $tz\n";
    } else {
      # compute T1 = pMHC transformation from template
      my $template_MHC = $template . "/pMHC.pdb";
      print "get_alignment_transformation( $pMHC, $template_MHC)\n";
      my ($rx1, $ry1, $rz1, $tx1, $ty1, $tz1) = get_alignment_transformation("../$pMHC", "$template_MHC");
      print "$rx1, $ry1, $rz1, $tx1, $ty1, $tz1\n";

      # compute T2 = TCR transformation to the template
      my $template_TCR = $template . "/TCR.pdb";
      my ($rx2, $ry2, $rz2, $tx2, $ty2, $tz2) = get_alignment_transformation("$template_TCR", "../$tcr_name");
      print "$rx2, $ry2, $rz2, $tx2, $ty2, $tz2\n";

      my ($rx, $ry, $rz, $tx, $ty, $tz) =
          Transformation::multiply_transforms($rx1, $ry1, $rz1, $tx1, $ty1, $tz1,
                                              $rx2, $ry2, $rz2, $tx2, $ty2, $tz2);
      print OUT "$counter $rx $ry $rz $tx $ty $tz\n";
    }

    chdir "..";

    $counter++;
  } # end templates iteration
  close (OUT);

  # run SOAP
  my $run_name = $pMHC_signature . "_" . $tcr_name;
  #my($filename, $directories, $suffix) = fileparse($run_name, qr/\.[^.]*/);
  #$run_name = $filename;

  my $cmd = "$imp_home/setup_environment.sh $imp_home/bin/soap_score $pMHC $tcr_name $trans_file -o $run_name.res";
  print "$cmd\n";
  `$cmd`;

  # pick the best template
  my $curr_score = `tail -n$counter $run_name.res | sort -nk3 | head -n1 | cut -d '|' -f2`;
  print "$counter $curr_score\n";
  chomp $curr_score;
  $curr_score += 0.0;
  print "Curr score = $curr_score\n";
  #`tail -n$counter $run_name.res >> $scores_files`;
  my $average_score = `tail -n$counter $run_name.res | awk '{sum+=\$3} END {print sum/NR}'`;
  chomp $average_score;
  my $stdev = `tail -n$counter $run_name.res | awk '{sum+=\$3; sumsq+=\$3*\$3} END {print sqrt(sumsq/NR - (sum/NR)**2)}'`;
  chomp $stdev;
  my $sum = $curr_score + $average_score;
  my $norm_score = ($curr_score-$average_score)/$stdev;
  #print SCORES_OUT "$tcr_name | $curr_score | $average_score | $sum | $norm_score\n";

  if(length($best_structure) == 0 || $best_score > $curr_score)  {
      my $best_num = `tail -n$counter $run_name.res | sort -nk3 | head -n1 | cut -d '|' -f1`;
      chomp $best_num;
      $best_num =~ s/\s+$//; $best_num =~ s/^\s+//;
      $best_structure = $run_name . "_" . $best_num . ".ref.pdb";
      #$best_structure = $tcr . "_" . $counter;
      print "best_structure = $best_structure\n";
      $best_score = $curr_score;
  }

  #`rm -f $run_name.res $trans_file`;
  #`rm -f log_multiprot.txt $trans_file modeledABC.pdb modeledDE.pdb soap_score.res`;

}  # end TCR iteration

# record best scoring model for the peptide
open OUT, ">>$pMHC_signature\_TCR_score_SOAP.txt";
print OUT "$pMHC_signature | $best_score | $best_structure\n";

`rm -rf $alignment_dirname`;

sub get_alignment_transformation {
 `$home/multiprot.Linux $_[0] $_[1]`;
 my $p = `pwd`;
 print $p;
  my $grep_line = `grep Trans 2_sol.res`;
  my @tmp = split(' ', $grep_line);
  unlink "2_sol.res";
  return ($tmp[$#tmp-5], $tmp[$#tmp-4], $tmp[$#tmp-3], $tmp[$#tmp-2], $tmp[$#tmp-1], $tmp[$#tmp]);
}

sub get_alignment_transformation2 {
  `match.linux $_[0] $_[1] | head -n2 > match`;
  #print "match.linux $_[0] $_[1]\n";
  my $grep_line = `grep RESULT match`;
  my @tmp = split(' ', $grep_line);
  unlink "match";
  return ($tmp[$#tmp-8], $tmp[$#tmp-7], $tmp[$#tmp-6], $tmp[$#tmp-5], $tmp[$#tmp-4], $tmp[$#tmp-3]);
}
