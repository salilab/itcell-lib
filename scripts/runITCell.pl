#!/usr/bin/perl -w

use strict;
use FindBin;
use Getopt::Long;

if($#ARGV < 2) {
    print "runITcell.pl <pMHC> <TCR> <seq>\n";
    exit;
}

# software directories
my $home = "$FindBin::Bin";
my $libs = "$home/../libs";

my $pMHC = $ARGV[0];
if(-e $ARGV[0]) { # uploaded file
# do something
} else {
  $pMHC =~ tr/\*/_/;
  $pMHC = "${libs}/MHCII_models/$pMHC/model.pdb";
}

my $TCR = $ARGV[1];
my $seq = $ARGV[2];


`$home/cleave.pl $seq`;
`$home/iterate_model_pMHC_SOAP.pl $seq $pMHC >& iterate_model_pMHC_SOAP.out`;
`echo TCR.pdb > TCR_list`;
`$home/iterate_model_TCR_pMHC_SOAP.pl ${libs}/templatesIItest/templates_DB TCR_list $seq $pMHC >& iterate_model_TCR_pMHC_SOAP.out`;
`$home/iterate_results.pl $seq $pMHC $seq.out | sort -gk5 | cat -n > scores.txt`;
