#!/usr/bin/python

import csv
import os
import shutil
import uuid
sys.path.append('./cell_bio_util')
import cross_fil_util as util

PROG_NAME   = 'cross_fil_map'
DESCRIPTION = 'CrossFil Python script to map FASTQ files then sort, quality filter and mark duplicate BAM output'
PICARD_TAG  = 'mdwm'
CLEAN_TAG   = 'f3_F4_q1'

ALIGNERS = ('bbmap', 'bwa', 'bt2')
ALIGNER_BBMAP, ALIGNER_BWA, ALIGNER_BT2 = ALIGNERS
DEFAULT_ALIGNER = ALIGNER_BT2
OTHER_ALIGNERS = [ALIGNER_BBMAP, ALIGNER_BWA]

def genome_map(aligner, strain_name, strain_num, fastq_paths, genome_index_path, num_cpu=util.MAX_CORES):  
  
  dir_name, base_name = os.path.split(fastq_paths[0])
  
  path_root         = os.path.join(dir_name, strain_name)
  sam_file_path     = '%s.sam' % path_root  
  
  if os.path.exists(sam_file_path):  
    util.info("SAM file %s already exists. Skipping genome mapping" % sam_file_path)
  
  else:
    util.info("Running aligner %s on %s..." % (aligner, strain_name))
 
    if aligner == ALIGNER_BWA:
      rg_header = "@RG\\tID:%s\\tSM:sample_%s\\tPL:illumina\\tLB:lib%d\\tPU:unit%d" % (strain_name, strain_name, strain_num, strain_num)
      cmd_args = [util.EXE[ALIGNER_BWA], 'mem',
                  '-t', str(num_cpu),
                  '-M',
                  '-R', rg_header,
                  genome_index_path] + list(fastq_paths)
      util.call(cmd_args, stdout=open(sam_file_path, 'w'))
   
    elif aligner == ALIGNER_BT2:
      cmd_args = [util.EXE[ALIGNER_BT2], '--sensitive', 
                  '-x', genome_index_path,
                  '-p', str(num_cpu),
                  '-q', # FASTQ input
                  '--rg-id', strain_name,
                  '--rg', "SM:sample_%s\tPL:illumina\tLB:lib%d\tPU:unit%d" % (strain_name, strain_num, strain_num),
                  '-S', sam_file_path]
      
      if len(fastq_paths) > 1:
        cmd_args += ['-1', fastq_paths[0], '-2', fastq_paths[1]]
      else:
        cmd_args += ['-U', fastq_paths[0]]
      
      util.call(cmd_args)
      
    else: # bbmap
      cmd_args = [util.EXE[ALIGNER_BBMAP],
                  'ref=%s' % genome_index_path,
                  'sam=1.3',
                  'in=%s' % fastq_paths[0],
                  'out=%s' % sam_file_path,
                  't=%d' % num_cpu,
                  'rgid=%s' % strain_name, 
                  'rgsm=sample_%s' % strain_name,
                  'rgpl=illumina', 
                  'rglb=lib%d' % strain_num, 
                  'rgpu=unit%d' % strain_num]
 
      if len(fastq_paths) > 1:
        cmd_args += ['in2=%s' % fastq_paths[1]]
 
      util.call(cmd_args)

  util.info('Done %s genome alignment for strain %s' % (aligner, strain_name))
  
  return sam_file_path
    

def sam_cleanup(sam_file_path, num_cpu=2):
  
  file_tag = util.FILE_TAG
  
  path_root, file_ext = os.path.splitext(sam_file_path)
  strain_name = os.path.basename(path_root)
  
  bam_file_path     = '%s%ssrt.bam' % (path_root, file_tag)
  clean_bam_path    = '%s%ssrt_%s.bam' % (path_root, file_tag, CLEAN_TAG)
  out_bam_path      = '%s%ssrt_%s_%s.bam' % (path_root, file_tag, CLEAN_TAG, PICARD_TAG)
  metrics_file_path = '%s%ssrt_%s_%s_metrics.txt' % (path_root, file_tag, CLEAN_TAG, PICARD_TAG)
  
  if os.path.exists(out_bam_path):
    util.info("BAM file %s already exists. Skipping SAM cleanup" % out_bam_path)
    return out_bam_path
   
  util.info("Converting SAM file from genome aligner output into sorted BAM...")
  
  cmd_args = [util.EXE['samtools'], 'sort',
              '-O', 'bam',
 #             '-@', str(num_cpu),
              '-o', bam_file_path,
              sam_file_path]
              
  util.call(cmd_args)
  
  util.info('Removing unmapped reads, PCR duplicates and low quality ones (MAPQ smaller than 1) keeping only paired reads which are properly mapped...') # Log strains individually
  
  cmd_args = ['samtools','view','-b','-f','3','-F','4','-q','1', bam_file_path]
  util.call(cmd_args, stdout=open(clean_bam_path, 'wb'))  
 
  util.info("Marking duplicate reads using Picard")
   
  cwd = os.getcwd()
  os.chdir('/') # Picard picky about relative paths
 
  cmd_args = list(util.JAVA)
  cmd_args += ['-jar', util.EXE['picard'],
               'MarkDuplicatesWithMateCigar',
               'I=%s' % clean_bam_path,
               'O=%s' % out_bam_path,
               'M=%s' % metrics_file_path]
 
  util.call(cmd_args)
 
  os.chdir(cwd)
 
  util.info("Indexing %s" % out_bam_path)
  util.call([util.EXE['samtools'],'index', out_bam_path])  

  util.info('Done BAM clean-up for strain %s' % strain_name)
  
  return out_bam_path


def bedtools_coverage(bam_file_path, genome_fasta_path, exon_gff_file_path):
  
  dir_name, base_name = os.path.split(bam_file_path)
  file_root = os.path.splitext(base_name)[0]
  dir_name = os.path.join(dir_name, 'coverage')
  
  util.makedirs(dir_name, exist_ok=True) # Not the os version to be Python 2 and 3 compatible
  
  genome_cvr_file_path = os.path.join(dir_name, base_name + '.genomecov')
  exon_cvr_file_path = os.path.join(dir_name, base_name + '_exon.coverage')
  exon_cvr_temp_file_path = os.path.join(dir_name, base_name + '_exon.coverage.temp')
  R_cvr_file_path = os.path.join(dir_name, base_name + '_R_coverage.out')
   
  if os.path.exists(R_cvr_file_path):
    util.info("Coverage file %s already exists. Skipping coverage calculations" % (R_cvr_file_path,))
    return
  
  bedtools_exe = util.EXE['bedtools']
    
  util.info("Running bedtools genomecov...")
  
  cmd_args = [bedtools_exe, 'genomecov', '-ibam', bam_file_path, '-g', genome_fasta_path]              
  util.call(cmd_args, stdout=genome_cvr_file_path)  
                     
  util.info("Done... Results saved in: %s" % genome_cvr_file_path)
  util.info("Converting %s into a sorted bed file..." % bam_file_path)
  
  temp_dir = os.path.join(dir_name, 'TEMP_%s' % uuid.uuid4())
  os.makedirs(temp_dir)
  
  temp_bed_file1 = os.path.join(temp_dir, '%s.bed' % file_root)
  temp_bed_file2 = os.path.join(temp_dir, '%s_sortBed.bed' % file_root)
  
  cmd_args = [bedtools_exe, 'bamtobed', '-i', bam_file_path]
  util.call(cmd_args, stdout=temp_bed_file1)   
             
  cmd_args = [bedtools_exe, 'sort', '-i', temp_bed_file1]
  util.call(cmd_args, stdout=temp_bed_file2)              
  
  util.info("Done... Results saved in temporary directory: %s" % temp_dir)
  util.info("Running bedtools coverage...")
  
  cmd_args = [bedtools_exe, 'coverage', '-hist' ,'-a', exon_gff_file_path,'-b', temp_bed_file2]
  util.call(cmd_args, stdout=exon_cvr_file_path)   
  
  util.info("Done... Results saved in: %s" % exon_cvr_file_path)

  # In order to calculate exon coverage in R, we need to extract all lines starting with "all" from exon.coverage file
  # This file is saved in a temporary directory
  
  cmd_args = ['grep', 'all', exon_cvr_file_path]
  util.call(cmd_args, stdout=exon_cvr_temp_file_path)   

  util.info("Running R to compute mean genome coverage and mean exon coverage..")
  
  cmd_args = ['Rscript', '--vanilla', util.EXE['mgcr'],
              genome_cvr_file_path, exon_cvr_temp_file_path]
  util.call(cmd_args, stdout=R_cvr_file_path)   

  util.info("Delete temporary directory and files...")
  shutil.rmtree(temp_dir)   

     
def cross_fil_map(barcode_csv, genome_fasta_path, exon_gff_path, fastq_paths_r1,
                  fastq_paths_r2=None, out_top_dir=None,  aligner=ALIGNER_BWA, bowtie2_index=None,
                  num_cpu=util.MAX_CORES, sub_dir_name=None, file_ext=None):
  
  if not sub_dir_name:
    sub_dir_name = 'strain_%s' % aligner

  if fastq_paths_r2:
    fastq_paths_r2_list = fastq_paths_r2
  else:
    fastq_paths_r2_list = []
  
  #for file_path in [barcode_csv, genome_fasta_path, exon_gff_path] + fastq_paths_r1 + fastq_paths_r2 or []:
  for file_path in [barcode_csv] + fastq_paths_r1 + fastq_paths_r2_list or []:
    is_ok, msg = util.check_regular_file(file_path)
  
    if not is_ok:
      util.critical(msg)
      
  if not file_ext:
    file_ext = util.get_file_ext(fastq_paths_r1[0])
    
  if not out_top_dir:
    file_path = fastq_paths_r1[0]
    out_top_dir = os.path.dirname(file_path)
    
  # # Concatenate FASTQ files 
  
  sample_barcodes = {}
  barcode_samples = {}
  
  # Read CVS
  with open(barcode_csv, 'rU') as file_obj:
    csv_data = csv.reader(file_obj)
    header = next(csv_data)
    
    for seq_run_id, barcode_name, barcode_seq, sample_name in csv_data:
      barcode_name = barcode_name.replace('-','_')
      sample_name = sample_name.replace(' ','_')
      
      if sample_name in barcode_samples:
        util.critical('Multiple samples/strains with name "%s" present in CSV file %s' % (sample_name, barcode_csv))
      
      barcode_samples[sample_name] = barcode_name
      sample_barcodes[barcode_name] = (seq_run_id, sample_name)
  
  # Make subdirs
  for sample_name in barcode_samples:
    dir_name = os.path.join(out_top_dir, sub_dir_name, sample_name)
    util.makedirs(dir_name, exist_ok=True)
  
  # Concatenate FASTQ files with the same barcode for each read 
  # and save results in corresponding strain folder
  # - now makes a symbolic link of only one file for a strain and not gzipped
  
  strain_fastq_paths = {}
  
  for barcode_name in sample_barcodes:
    seq_run_id, sample_name = sample_barcodes[barcode_name]
    file_pattern = '%s*%s*' % (seq_run_id, barcode_name)
    
    if fastq_paths_r2:
      out_file_name_1 = '%s_r_1%s' % (sample_name, file_ext)
      out_file_name_2 = '%s_r_2%s' % (sample_name, file_ext)
      in_fastq_paths_1 = util.match_files(fastq_paths_r1, file_pattern) # Read pairs files already separated
      in_fastq_paths_2 = util.match_files(fastq_paths_r2, file_pattern)
      
      io_paths =  [(in_fastq_paths_1, out_file_name_1),
                   (in_fastq_paths_2, out_file_name_2)] 
  
    else:
      out_file_name  = '%s%s' % (sample_name, file_ext)
      in_fastq_paths = util.match_files(fastq_paths_r1, file_pattern)
      
      io_paths =  [(in_fastq_paths, out_file_name)]
    
    fastq_paths = []
       
    for in_fastq_paths, out_file_name in io_paths:
      if not in_fastq_paths:
        util.critical('No FASTQ read files found for run %s barcode %s' % (seq_run_id, barcode_name))
      
      out_fastq_path = os.path.join(out_top_dir, sub_dir_name, sample_name, out_file_name)

      if os.path.exists(out_fastq_path):
        util.warn('FASTQ file %s already exists and won\'t be overwritten...' % out_fastq_path)
      
      else:
        # Concatenate or sym link
        
        if len(in_fastq_paths) == 1 and not in_fastq_paths[0].endswith('.gz'):
          util.info('Sym linking %s reads to %s' % (barcode_name, out_fastq_path))
          os.symlink(in_fastq_paths[0], out_fastq_path)
        
        else:
          with open(out_fastq_path, 'wb') as out_file_obj:
            util.info('Concatenating %s reads to %s' % (barcode_name, out_fastq_path))
            
            for fastq_path in in_fastq_paths:
              shutil.copyfileobj(util.open_file(fastq_path, 'rb'), out_file_obj) # Accepts GZIP input
 
      fastq_paths.append(out_fastq_path)
    
    strain_fastq_paths[sample_name] = fastq_paths  
  
  # Genome alignment/mapping
  if aligner == ALIGNER_BT2:
    genome_index = bowtie2_index
  else:
    genome_index = genome_fasta_path
  
  sam_paths = []
  for i, strain_name in enumerate(strain_fastq_paths):
    sam_path = genome_map(aligner, strain_name, i+1, strain_fastq_paths[strain_name], genome_index, num_cpu)
    sam_paths.append(sam_path)  
  
  # Parallel BAM cleanup
  # samtools sort can use multiple cores but is generally I/O cound in any case
  bam_paths = util.parallel_split_job(sam_cleanup, sam_paths, [], num_cpu, collect_output=True)
  
  # Parallel coverage with BEDtools
  common_args = [genome_fasta_path, exon_gff_path]
  util.parallel_split_job(bedtools_coverage, bam_paths, common_args, num_cpu, collect_output=False)
      
  util.info('%s finished for %d input strains' % (PROG_NAME, len(strain_fastq_paths)))  
  
  
   
if __name__ == '__main__':

  from argparse import ArgumentParser
   
  epilog = 'For further help on running this program please email tjs23@cam.ac.uk.\n\n'
  epilog += 'Example use:\n\n'
  epilog += 'python3 cross_fil_map.py /data/SLX-12506/SLX-12506.valperga.contents.csv '
  epilog += '/home/paulafp/Documents/temp/WS255_WBcel235/uncompressed_fa/c_elegans.PRJNA13758.WS255.genomic.fa  '
  epilog += '/home/paulafp/Documents/temp/WS255_WBcel235/c_elegans.PRJNA13758.WS255.WormbaseExons_sorted.gff3 '
  epilog += '/data/SLX-12506/trimmed/SLX-12506*.fq -pe r_1_val_1 r_2_val_2 '
  epilog += '-bt2_index /data/genome_builds/c_elegans_WS259/c_elegans_WS259'
  
  arg_parse = ArgumentParser(prog=PROG_NAME, description=DESCRIPTION,
                             epilog=epilog, prefix_chars='-', add_help=True)
  
  arg_parse.add_argument('barcode_csv', metavar='BARCODE_CSV_FILE',
                         help='CSV format file containing barcode strain/sample names') 

  arg_parse.add_argument('genome_fasta', metavar='GENOME_FASTA',
                         help='File path of genome sequence FASTA file (for use by genome aligner)') 

  arg_parse.add_argument('exon_gff_path', metavar='GFF3_FILE',
                         help='File path of GFF3 format file containing exon information') 

  arg_parse.add_argument('fastq_paths', nargs='+', metavar='FASTQ_FILES',
                         help='File paths of FASTQ sequence read files (may contain wildcards) ') 
                         
  arg_parse.add_argument('-bt2_index', metavar='BOWTIE2_GENOME_INDEX', default=None,
                         help='File path of genome index files for Bowtie2 without any file extension (.1.bt2 etc.)') 

  arg_parse.add_argument('-outdir', metavar='DIR_NAME', default=None,
                         help='Name of directory for output files and sub-directories (need not exist). Defaults to where first FASTQ file is located') 

  arg_parse.add_argument('-cpu', metavar='NUM_CORES', default=util.MAX_CORES, type=int,
                         help='Number of parallel CPU cores to use. Default: All available (%d)' % util.MAX_CORES) 

  arg_parse.add_argument('-al', metavar='ALIGNER_NAME', default=DEFAULT_ALIGNER,
                         help='Name of the program to perform the genome alignment/mapping: Default: %s Other options: %s' % (DEFAULT_ALIGNER, OTHER_ALIGNERS)) 
 
  arg_parse.add_argument('-pe', nargs=2, metavar='PAIRED_READ_TAGS', default=['r_1','r_2'],
                         help='The subtrings/tags which are the only differences between paired FASTQ file paths. Default: r_1 r_2') 
                         
  arg_parse.add_argument('-se', default=False, action='store_true',
                         help='Input reads are single-end data, otherwise defaults to paired-end.')

  arg_parse.add_argument('-q', default=False, action='store_true',
                         help='Sets quiet mode to supress on-screen reporting.')
  
  arg_parse.add_argument('-log', default=False, action='store_true',
                         help='Log all reported output to a file.')

  args = vars(arg_parse.parse_args())

  barcode_csv   = args['barcode_csv']
  fastq_paths   = args['fastq_paths']
  genome_fasta  = args['genome_fasta']
  bowtie2_index = args['bt2_index']
  exon_gff_path = args['exon_gff_path']
  aligner       = args['al']
  pair_tags     = args['pe']
  is_single_end = args['se']
  out_top_dir   = args['outdir']
  num_cpu       = args['cpu'] or None # May not be zero

  # Reporting handled by cross_fil_util
  util.QUIET   = args['q']
  util.LOGGING = args['log']
  
  if aligner:
    aligner = aligner.lower()
    
    if aligner not in ALIGNERS:
      util.critical('Aligner option must be one of: %s' % ', '.join(ALIGNERS))
 
    if bowtie2_index and aligner != ALIGNER_BT2:
      util.critical('Aligner option must be "%s" (or left unset) if a Bowtie2 genome index is specified' % ALIGNER_BT2)
        
  else:
    if bowtie2_index:
      aligner = ALIGNER_BT2
    else:
      aligner = ALIGNER_BWA
  
  if len(pair_tags) != 2:
    util.critical('When specified, exactly two paired-end filename tags must be given.')
  
  if is_single_end:
    fastq_paths_r1 = fastq_paths
    fastq_paths_r2 = None
  
  else:
    if not fastq_paths:
      util.critical('No FASTQ files found')

    elif len(fastq_paths) % 2 != 0:
      util.critical('When using paired-end data an even number of FASTQ files must be specified.')
  
    fastq_paths_r1, fastq_paths_r2 = util.pair_fastq_files(fastq_paths, pair_tags)
     
  cross_fil_map(barcode_csv, genome_fasta, exon_gff_path, fastq_paths_r1, fastq_paths_r2,
                out_top_dir, aligner, bowtie2_index, num_cpu)
