from modeller import *
from modeller.automodel import *
from modeller import soap_protein_od
from modeller.scripts import complete_pdb

import subprocess
import os,sys

current_home =  os.path.dirname(sys.argv[0])
print current_home;

blast_home = "/salilab/diva1/programs/linux/blast/blast-latest/bin/";

def run_blast(filename):
    out_file = filename + ".blast"
    blast_file = current_home + "/../libs/blastdb/pdbaa ";
    cmd = blast_home + "blastpgp -d " + blast_file + " -i " + filename + " >& " + out_file
    print cmd
    # Put stderr and stdout into pipes
    proc = subprocess.Popen(cmd, shell=True, stdout=sys.stdout, stderr=sys.stderr)
    return_code = proc.wait()
    if return_code != 0:
        print "blast subprocess failed with exit code " , return_code

if(len(sys.argv) != 2) :
    print "Please provide input fasta file";
    sys.exit(0)

fasta_file_name = sys.argv[1]
run_blast(fasta_file_name)
blast_file_name = fasta_file_name + ".blast"
count = 0
template_list = []
for line in open(blast_file_name):
    if "|" in line and count <10:
        fields = line.split(' ')
        fields2 = fields[0].split('|')
        #print fields2[3],fields2[4]
        template = fields2[3].lower() + fields2[4]
        template_list.append(template)
        count+=1


print template_list;

log.verbose()
env = environ()
env.io.atom_files_directory = current_home + '/../libs/pdbs/'
env.libs.topology.read(file='$(LIB)/top_heav.lib')
env.libs.parameters.read(file='$(LIB)/par.lib') # read parameters

aln = alignment(env)
for (template_code) in (template_list):
    code = template_code[:4]
    chain = template_code[-1:]
    #print code;
    #print chain;
    mdl = model(env, file=code, model_segment=('FIRST:'+chain, '+120:'+chain))
    aln.append_model(mdl, atom_files=code, align_codes=code+chain)

for (weights, write_fit, whole) in (((1., 0., 0., 0., 1., 0.), False, True),
                                    ((1., 0.5, 1., 1., 1., 0.), False, True),
                                    ((1., 1., 1., 1., 1., 0.), True, False)):
    aln.salign(rms_cutoff=3.5, normalize_pp_scores=False,
               rr_file='$(LIB)/as1.sim.mat', overhang=30,
               gap_penalties_1d=(-450, -50),
               gap_penalties_3d=(0, 3), gap_gap_score=0, gap_residue_score=0,
               alignment_type='tree', # If 'progresive', the tree is not
                                      # computed and all structues will be
                                      # aligned sequentially to the first
               feature_weights=weights, # For a multiple sequence alignment only
                                        # the first feature needs to be non-zero
               improve_alignment=True, fit=True, write_fit=write_fit,
               write_whole_pdb=whole, output='ALIGNMENT QUALITY')

aln.write(file='templates.ali', alignment_format='PIR')
aln_block = len(aln)

# Read aligned sequence(s):
aln.append(file='TCR.ali', align_codes='TCR')

# Structure sensitive variable gap penalty sequence-sequence alignment:
aln.salign(output='', max_gap_length=10,
           gap_function=True,   # to use structure-dependent gap penalty
           alignment_type='PAIRWISE', align_block=aln_block,
           feature_weights=(1., 0., 0., 0., 0., 0.), overhang=99,
           gap_penalties_1d=(-450, 0),
           gap_penalties_2d=(0.35, 1.2, 0.9, 1.2, 0.6, 8.6, 1.2, 0., 0.),
           similarity_flag=True, local_alignment=True)

aln.write(file='TCR-mult.ali', alignment_format='PIR')

f = open('gnuplotfile', 'w')
f.write('plot ')

sp = soap_protein_od.Scorer()
for (template_code) in (template_list):
    code = template_code[:4]
    mdl = complete_pdb(env, code + '_fit.pdb')
    s = selection(mdl)   # all atom selection
    score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT', file=code + '.profile',
                     normalize_profile=True, smoothing_window=15)
    f.write('\''+ code + '.profile\' u 1:42 w l lc rgb \'#61a0e2\',')


a = automodel(env, alnfile='TCR-mult.ali',
              knowns=template_list,
              sequence='TCR',assess_methods=(assess.DOPE,sp))
a.starting_model = 1
a.ending_model = 10
a.make()

for n in range(1, 9):
    model_name = 'TCR.B9999000' + str(n) + '.pdb'
    mdl = complete_pdb(env, model_name)
    s = selection(mdl)   # all atom selection
    score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT', file='TCR.B9999000' + str(n) + '.profile',
                     normalize_profile=True, smoothing_window=15)
    f.write('\'TCR.B9999000' + str(n) + '.profile\' u 1:42 w l lc rgb \'#e26261\',')


model_name = 'TCR.B99990010.pdb'
mdl = complete_pdb(env, model_name)
s = selection(mdl)   # all atom selection
score = s.assess(sp, output='ENERGY_PROFILE NO_REPORT', file='TCR.B99990010.profile',
                     normalize_profile=True, smoothing_window=15)
f.write('\'TCR.B99990010.profile\' u 1:42 w l lc rgb \'#e26261\'')
