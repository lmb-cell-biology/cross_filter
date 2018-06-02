#!/usr/bin/python

import os
import sys
current_path = os.path.realpath(__file__)
current_path = os.path.dirname(current_path) + '/cell_bio_util'
sys.path.append(current_path)
import cell_bio_util as util
import cross_fil_exe as exe

for name in exe.EXE:
  # print(exe.EXE)
  util.check_exe(exe.EXE[name])

PROG_NAME   = 'cross_fil_genotype'
DESCRIPTION = 'CrossFil Python script to generate genotype VCF files'

util.init_app('cf') # Redefine variables from cross_fil_util.py

# Below should be available in exe.EXE
CALLER_FREEBAYES = 'freebayes'
CALLER_GATK      = 'gatk'

VAR_CALLERS = (CALLER_FREEBAYES, CALLER_GATK) # First is the default


def gatk_haplotype_job(bam_file_path, genome_fasta_path, sub_dir_name):
  
  path_root, file_ext = os.path.splitext(bam_file_path)
  dir_name, file_root = os.path.split(path_root)
  
  vcf_dir_name = os.path.join(dir_name, sub_dir_name)
  # vcf_file_path  = os.path.join(vcf_dir_name, '%s_hap.vcf' % (file_root))
  gvcf_file_path  = os.path.join(vcf_dir_name, '%s_hap.g.vcf' % (file_root))
    
  if os.path.exists(gvcf_file_path):
    util.info("VCF file %s already exists. Skipping haplotype calling for %s" % (gvcf_file_path, file_root))
    return gvcf_file_path
  
  util.makedirs(vcf_dir_name, exist_ok=True)
  
  util.info('Creating GVCF file for %s using GATK' % file_root)

  cmd_args = list(util.JAVA) + ['-jar', exe.EXE[CALLER_GATK],
                                '-T', 'HaplotypeCaller',
                                '-R', genome_fasta_path,
                                '-I', bam_file_path,
                                '-o', gvcf_file_path,
                                '-ERC', 'GVCF',
                                '-variant_index_type', 'LINEAR',             # Deprecated for GATK 4.0
                                '-variant_index_parameter', '128000']        # Deprecated for GATK 4.0
  util.call(cmd_args)
  
  return gvcf_file_path


def _get_merged_vcf_path(dir_name, strain_names, tag):
  
  sort_strains = sorted(strain_names)
  
  if len(strain_names) < 7:
    strain_text = '_'.join(sort_strains)
    merge_file_path = os.path.join(dir_name, 'merged_%s_comb_%s.vcf' % (strain_text, tag))
  
  else:
    num_strains = len(sort_strains)
    merge_file_path = os.path.join(dir_name, 'merged_%d_comb_%s.vcf' % (num_strains, tag))
           
  util.info('Strains to be combined:\n%s\n' % ' '.join(sort_strains))

  
  """
  if os.path.exists(merge_file_path):
    util.warn('%s already exists and won\'t be overwritten...' % merge_file_path)
    
    i = 0
    merge_file_path = '%s_%03d.vcf' % (merge_file_path[:-4], i)
    
    while os.path.exists(merge_file_path):
      i += 1
      merge_file_path = '%s_%03d.vcf' % (merge_file_path[:-8], i)
      
    util.info('Results will be saved in %s' % merge_file_path)
 
  """
    
  return merge_file_path
    
    
def gatk_merge_vcfs(dir_name, strain_vcf_paths, genome_fasta_path, num_cpu=util.MAX_CORES):

  merge_file_path = _get_merged_vcf_path(dir_name, strain_vcf_paths.keys(), CALLER_GATK)

  cmd_args = list(util.JAVA) + ['-jar', exe.EXE[CALLER_GATK],
                                '-T', 'GenotypeGVCFs',
                                '-R', genome_fasta_path,
                                #'-nt', str(min(8, num_cpu)), # Seems to fail with multiple CPU threads...
                                '-o', merge_file_path]

  for strain in strain_vcf_paths:
    cmd_args += ['-V', strain_vcf_paths[strain]]
  
  util.call(cmd_args)

  return merge_file_path


def gatk_select_homozygous_vars(strain_name, merged_vcf_path, genome_fasta_path, file_tag='extracted'):

  dir_name, file_name = os.path.split(merged_vcf_path)
  
  out_vcf_path = os.path.join(dir_name, '%s_%s.vcf' % (strain_name, file_tag))
  # Original naming: "%s_sorted_f3_F4_q1_mark_dups_w_mate_cig_gatk_hap_call_extracted.vcf" % strain_name
 
  util.info('Creating VCF file for %s' % strain_name)

  cmd_args = list(util.JAVA) + ['-jar', exe.EXE[CALLER_GATK],
                               '-T', 'SelectVariants',
                               '-R', genome_fasta_path,
                               '-V', merged_vcf_path,
                               '-o', out_vcf_path,
                               '-sn', 'sample_%s' % strain_name,
                               '-select', "vc.getGenotype('sample_%s').isHomVar()" % strain_name] # Check quotes

  util.call(cmd_args)
  util.info('All done for strain %s. VCF file can be found in %s' % (strain_name, out_vcf_path))
  

def call_genotype_gatk(strain_bam_paths, genome_fasta_path, num_cpu, out_dir, sub_dir_name):
  # GATK pipeline - parallelise strains in python
  
  genome_index_file = genome_fasta_path + '.fai'
  genome_dict_file  = os.path.splitext(genome_fasta_path)[0] + '.dict'
  
  if not os.path.exists(genome_index_file):
    util.info('Making index for genome FASTA file %s' % genome_fasta_path)
    cmd_args = [exe.EXE['samtools'], 'faidx', genome_fasta_path]
    util.call(cmd_args) 

  if not os.path.exists(genome_dict_file):
    util.info('Making dict filr for genome FASTA file %s' % genome_fasta_path)
    cmd_args = [exe.EXE['samtools'], 'dict', genome_fasta_path, '-o', genome_dict_file]
    util.call(cmd_args) 
  
  strains = sorted(strain_bam_paths)
  
  bam_paths    = [strain_bam_paths[s] for s in strains] # Each parallel call will be sent one of these
  common_args  = [genome_fasta_path, sub_dir_name]      # All tasks share this
  
  vcf_paths = util.parallel_split_job(gatk_haplotype_job, bam_paths, common_args, num_cpu, collect_output=True)
  # BAM and VCF path will be in corresponding order
    
  # Multi-sample
  strain_vcf_paths = {strains[i]:p for i, p in enumerate(vcf_paths)}
  merged_vcf_path = gatk_merge_vcfs(out_dir, strain_vcf_paths, genome_fasta_path, num_cpu=num_cpu)    
  
  return merged_vcf_path


def freebayes_genotype_job(region, genome_fasta_path, bam_paths):
  
  out_vcf_path = 'temp_%s_freebayes.vcf' % region
  
  if not os.path.exists(out_vcf_path):
 
    cmd_args = [exe.EXE['freebayes'],
                '--no-mnps', # make this optional
                '--no-complex', # make this optional
                '-f', genome_fasta_path,
                '-r', region,
                '-v', out_vcf_path] #, '--ploidy', '2']
 
    cmd_args += bam_paths
 
    util.call(cmd_args)

  return out_vcf_path
  

def call_genotype_freebayes(strain_bam_paths, genome_fasta_path, num_cpu, out_dir, sub_dir_name):
  # FreeBayes pipeline
  
  
  strain_names, bam_file_paths = zip(*list(strain_bam_paths.items()))
  
  merge_file_path = _get_merged_vcf_path(out_dir, bam_file_paths, CALLER_FREEBAYES)
  temp_file_path_a = util.get_temp_path(merge_file_path)
  temp_file_path_b = util.get_temp_path(merge_file_path)

  # Make regions for parallelisation, splitting all chromos according to number of CPUs
  
  chromo_sizes = util.get_bam_chromo_sizes(bam_file_paths[0])
      
  regions = []
  region_fmt = '%s:%d-%d'
  
  for chromo, size in chromo_sizes:
    step = int(size/num_cpu) + 1 # will be rounded up
    
    i = 0
    j = step
    
    while j < size:
      regions.append(region_fmt % (chromo, i, j))
      i = j
      j += step
    
    regions.append(region_fmt % (chromo, i, size))
  
  # Call haplotype for all strains at once, split into parallel regions
  
  common_args = [genome_fasta_path, bam_file_paths]
  region_vcf_paths = util.parallel_split_job(freebayes_genotype_job, regions, common_args,
                                             num_cpu, collect_output=True)
  
  # Combine the regions which were run in parallel
  
  util.info('Combining freebayes regions')
  out_file_obj = open(temp_file_path_a, 'w')
  write = out_file_obj.write
  
  for i, region_vcf in enumerate(region_vcf_paths):
    with open(region_vcf) as file_obj:
      for line in file_obj:
        if line[0] == '#':
          if i == 0:
            write(line)
        
        else:
          write(line)
  
  out_file_obj.close()
  cmd_args = [exe.EXE['vcfuniq']]
  util.call(cmd_args, stdin=temp_file_path_a, stdout=merge_file_path)
  
  # Cleanup temp files
  
  os.unlink(temp_file_path_a)

  for file_path in region_vcf_paths:
    os.unlink(file_path)
 
  return merge_file_path


def cross_fil_genotype(bam_file_paths, genome_fasta_path, var_caller=CALLER_FREEBAYES,
                       num_cpu=util.MAX_CORES, out_dir=None, sub_dir_name=None):
  
  import subprocess
  
  # Main function to call genotypes given (cleaned) input BAM files with the option of different caller programs
  # Requires filtered, duplicate marked BAM file for each sample/strain
  if not sub_dir_name:
    sub_dir_name = 'vcf_%s' % var_caller
  
  if not out_dir:
    out_dir = os.getcwd()
  
  for file_path in bam_file_paths + [genome_fasta_path]:
    is_ok, msg = util.check_regular_file(file_path)
  
    if not is_ok:
      util.critical(msg)
  
  strain_bam_paths = {}
  for bam_file_path in sorted(bam_file_paths):
    file_name = os.path.basename(bam_file_path)
    
    if util.FILE_TAG not in file_name:
      msg = 'BAM file name %s does not contain %s; it does not appear to have been created by cross_fil_map' % (file_name, util.FILE_TAG)
      util.critical(msg)
        
    sample_name = file_name.split(util.FILE_TAG)[0]
    strain_bam_paths[sample_name] = bam_file_path
  
  strain_names = list(strain_bam_paths.keys())
  num_strains = len(strain_names)  
  # Code for different varient callers splits early as the pipelines differ somewhat
  
  if var_caller == CALLER_GATK:
    # Create separate gVCFS, combine and call genotype
    merged_vcf_path = call_genotype_gatk(strain_bam_paths, genome_fasta_path, num_cpu, out_dir, sub_dir_name)
  
  else:
    # Call genotype directly on multiple BAM files
    merged_vcf_path = call_genotype_freebayes(strain_bam_paths, genome_fasta_path, num_cpu, out_dir, sub_dir_name)
  
  
  # Select homozygous variants for each strain from combined genotype VCF file - parallelise strains in python
  
  file_tag = 'homozy_%s' % var_caller
  common_args = [merged_vcf_path, genome_fasta_path, file_tag]
  util.parallel_split_job(gatk_select_homozygous_vars, strain_names, common_args, num_cpu)
    
  util.info('%s finished for %d input strains' % (PROG_NAME, num_strains))  



if __name__ == '__main__':

  from argparse import ArgumentParser
   
  epilog = 'For further help on running this program please email tjs23@cam.ac.uk.\n\n'
  epilog += 'Example use:\n\n'
  epilog += 'python3 cross_fil_genotype.py /data/SLX-12506/trimmed/strain/*/*._sorted_f3_F4_q1_mdwm.bam '
  epilog += '/home/paulafp/Documents/temp/WS255_WBcel235/uncompressed_fa/c_elegans.PRJNA13758.WS255.genomic.fa '
  epilog += '-outdir /data/SLX-12506/trimmed/strain/'
  epilog += ' -vc freebayes'
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('bam_paths', nargs='+', metavar='BAM_FILES',
                         help='Input file paths of cleaned BAM files as output by cross_fil_map (may contain wildcards) ') 

  arg_parse.add_argument('genome_fasta', metavar='GENOME_FASTA',
                         help='File path of genome sequence FASTA file') 
  
  default_var_caller = VAR_CALLERS[0]
  other_var_callers  = ', '.join(VAR_CALLERS[1:])
  
  arg_parse.add_argument('-vc', metavar='CALLER_NAME', default=None,
                         help='Name of the program to perform variant calling: Default: %s Other options: %s' % (default_var_caller, other_var_callers)) 

  arg_parse.add_argument('-outdir', metavar='DIR_NAME', default=None,
                         help='Optional name of directory for output VCF files. Default is current working directory') 

  arg_parse.add_argument('-cpu', metavar='NUM_CORES', default=util.MAX_CORES, type=int,
                         help='Number of parallel CPU cores to use. Default: All available (%d)' % util.MAX_CORES) 

  arg_parse.add_argument('-q', default=False, action='store_true',
                         help='Sets quiet mode to supress on-screen reporting.')
  
  arg_parse.add_argument('-log', default=False, action='store_true',
                         help='Log all reported output to a file.')
  
  args = vars(arg_parse.parse_args())

  genome_fasta   = args['genome_fasta']
  bam_file_paths = args['bam_paths']
  num_cpu        = args['cpu'] or None # May not be zero
  out_dir        = args['outdir']
  var_caller     = args['vc']

  # Reporting handled by cross_fil_util
  util.QUIET   = args['q']
  util.LOGGING = args['log']  
  
  if var_caller:
    var_caller = var_caller.lower()
  
    if var_caller not in VAR_CALLERS:
      util.critical('Variant caller option must be one of: %s' % ', '.join(VAR_CALLERS)) 
  
  else:
    var_caller = default_var_caller
  
  
  if out_dir and not os.path.exists(out_dir):
    util.critical('Output directory "%s" does not exist' % out_dir)  
  
  cross_fil_genotype(bam_file_paths, genome_fasta, var_caller, num_cpu, out_dir)
