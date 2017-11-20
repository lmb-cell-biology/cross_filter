#!/usr/bin/python

# Make a single background VCF file for combined strains/samples which can be subtracted from any individual sample

#!/usr/bin/python

import os
sys.path.append('./cell_bio_util')
import cross_fil_util as util

PROG_NAME   = 'cross_fil_background'
DESCRIPTION = 'CrossFil Python script to make a single background VCF file for combined strains/samples which may then be subtracted from any individual sample using cross_fil_subtract'
DEFAULT_MIN_NUM_OBS = 3

def cross_fil_background(strain_vcf_files, out_vcf_path=None, min_num_obs=3):
    
  for file_path in strain_vcf_files:
    is_ok, msg = util.check_regular_file(file_path)
  
    if not is_ok:
      util.critical(msg)
  
  if not out_vcf_path:
    out_vcf_path = 'bg_s%d_m%d.vcf' % (len(strain_vcf_files), min_num_obs)
  
  file_root, file_ext = os.path.splitext(out_vcf_path)
  
  comb_vcf_path = '%s_comb_input.vcf' % file_root
  temp_comb_vcf_path = util.get_temp_path(comb_vcf_path)
  
  util.info('Creating background VCF file for %d input files' % (len(strain_vcf_files)))
  
  # Combine each strain's diploid variants into a combined VCF
  # vcfcombine: Combine multiple VCF files together, handling samples when alternate allele descriptions are identical
  # vcfintersect -u any better?
   
  cmd_args = [util.EXE['vcfcombine']]
  cmd_args += strain_vcf_files
  util.call(cmd_args, stdout=temp_comb_vcf_path)
  
  # Check chromosome sorting 
  
  cmd_args = [util.EXE['vcfstreamsort'], '-a']
  util.call(cmd_args, stdin=temp_comb_vcf_path, stdout=comb_vcf_path)
  
  os.unlink(temp_comb_vcf_path)
  
  # Filter on number of samples represented in the genotype fields
  # i.e. for vars that occur at least a given number of times across strainss
  
  out_file_obj = open(out_vcf_path, 'w')
  write = out_file_obj.write
  num_samples = None # Filled from header info
  
  with open(comb_vcf_path) as file_obj:
    for line in file_obj:
      if line[0] == '#':
        if line[1:6] == 'CHROM':
          header = line.split()
          
          if len(header) < 10:
            util.critical('Cannot filter sample genotypes in a VCF file without FORMAT and sample/genotype information')
           
          else: 
            sample_names = header[9:] #Could use in future to track which strains are selected
            num_samples = len(sample_names)
            
        write(line)
      
      else:
        data = line.split()      
        genotypes = data[9:]
        num_obs = num_samples - genotypes.count('.')
        
        if num_obs >= min_num_obs:
          write(line)      
  
  util.info('Background VCF file output at "%s"' % (out_vcf_path, ))
  
  
if __name__ == '__main__':

  from argparse import ArgumentParser
   
  epilog = 'For further help on running this program please email tjs23@cam.ac.uk.\n\n'
  epilog += 'Example use:\n\n'
  epilog += 'python3 cross_fil_background.py /data/SLX-12506/trimmed/strain/*_extracted.vcf '
  epilog += '-o /data/SLX-12506/trimmed/strain/background_1.vcf -m 3'
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)

  arg_parse.add_argument('vcf_paths', nargs='+', metavar='VCF_FILES',
                         help='Input file paths of extracted strain VCF files as created by cross_fil_genotype (may contain wildcards) ') 

  arg_parse.add_argument('-o', metavar='VCF_FILE', default=None,
                         help='Optional output path for combined backgrouund VCF file. Default is "bg_s{NUM_STRAINS}_m{MIN_VAR_COUNT}.vcf" in the current working directory') 

  arg_parse.add_argument('-m', metavar='MIN_VAR_COUNT', default=DEFAULT_MIN_NUM_OBS, type=int,
                         help='Minimum number of variant occurrences (in input strains) required for acceptance into background file. Default: %d' % DEFAULT_MIN_NUM_OBS) 
  
  arg_parse.add_argument('-q', default=False, action='store_true',
                         help='Sets quiet mode to supress on-screen reporting.')
  
  arg_parse.add_argument('-log', default=False, action='store_true',
                         help='Log all reported output to a file.')
  
  args = vars(arg_parse.parse_args())

  vcf_paths = args['vcf_paths']
  out_path  = args['o']
  min_num_obs = args['m']

  # Reporting handled by cross_fil_util
  util.QUIET   = args['q']
  util.LOGGING = args['log']  
  
  cross_fil_background(vcf_paths, out_path, min_num_obs)
