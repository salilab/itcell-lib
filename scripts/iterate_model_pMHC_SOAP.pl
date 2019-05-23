#!/usr/bin/perl -w

use strict;
use FindBin;
use File::Basename;

my $home = "$FindBin::Bin";
require "$home/UTIL.pm";

if ($#ARGV != 1) {
  print "Usage: iterate_model_pMHC.pl <peptide_sequence_file> <pMHC_template>\n";
  print "Assumes that template MHC chains are AB, and peptide chain is C\n";
  exit;
}

my $imp_home = "/netapp/sali/dina/imp_server_build2/";

my $sequence_file = $ARGV[0];
my $pMHC = $ARGV[1];
my $pMHC_signature = UTIL::pMHC_name_signature($pMHC);

my @sequences = UTIL::read_fasta_sequences($sequence_file);

my $peptide = `$home/pdb2fasta $pMHC | tail -n1`;
chomp $peptide;
my $peptide_length = length($peptide);
print "MHC peptide_length = $peptide_length\n";

# we assume MHC chains A and B
my $MHC = "pMHC_" . $pMHC_signature . ".pdb";
`$home/getChain.Linux AB $pMHC > $MHC`;
my $mhc_fasta_file = $pMHC_signature . ".seq";
`$home/pdb2fasta $MHC | grep -v \">\" | perl -CSD -ne \'print lc\' > $mhc_fasta_file`; # lowercase

my $filenames_file = UTIL::trim_extension($sequence_file) . "_" . $pMHC_signature . "_filenames.txt";
open OUT, ">$filenames_file";
my $window_size = $peptide_length;
foreach my $sequence (@sequences) {
  for(my $i=0; $i<length($sequence)-($window_size-1); $i++) {
    my $peptide_sequence = substr($sequence, $i, $window_size);
    my $name_signature = $pMHC_signature . "_" . $peptide_sequence;

    # 1.0 scwrl
    #create sequence file for scwrl
    my $seq_file = $name_signature . ".seq";
    `cat $mhc_fasta_file > $seq_file`;
    `echo $peptide_sequence >> $seq_file`; # upper case
    my $out_file = "pMHC_" . $name_signature . ".pdb";
    print "$home/scwrl4/Scwrl4 -i $pMHC -s $seq_file -o $out_file\n";
    `$home/scwrl4/Scwrl4 -i $pMHC -s $seq_file -o $out_file`;

    my $peptide_file = $name_signature . ".pdb";
    `$home/getChain.Linux C $out_file > $peptide_file`;
    unlink $seq_file;

    print OUT "$MHC $peptide_file\n";
  }
}
close OUT;

my $out_file_name = UTIL::trim_extension($sequence_file) . "_" . $pMHC_signature;
my $out_file = $out_file_name . ".soap_mhc";

my $cmd = "$imp_home/setup_environment.sh $imp_home/bin/soap_score $filenames_file -r -f $home/../libs/soap/cs1.hdf5 --ligand_potentials_file $home/../libs/soap/bs2.hdf5 -o $out_file -n";
print "$cmd\n";
`$cmd`;

$out_file = $out_file_name . ".soap_pep";
$cmd = "$imp_home/setup_environment.sh $imp_home/bin/soap_score $filenames_file -r -f $home/../libs/soap/soap_pep.hdf5 -o $out_file";
print "$cmd\n";
`$cmd`;

my $resfile = $out_file_name . ".soap";
`cp $out_file $resfile`;

#unlink $MHC;
