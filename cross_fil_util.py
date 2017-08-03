#!/usr/bin/python

import gzip
import multiprocessing
import os
import random
import string
import subprocess
import sys
import uuid

QUIET   = False
LOGGING = False

FILE_TAG = '_cf_' # This tag is used for formatting file names so they can be passed between the various cross_fil programs 

JAVA = ['java', '-d64', '-Xmx4g'] # 4 gigabyte heap size, for GATK, Picard etc.

TEMP_ID = '%s' % uuid.uuid4()
LOG_FILE_PATH = 'cf-out-%s.log' % TEMP_ID
LOG_FILE_OBJ = None # Created when needed
MAX_CORES = multiprocessing.cpu_count()

VCF_LIB_DIR = '/home/tjs23/apps/freebayes/vcflib/bin/'
       
EXE = {#'bbmap'     :'bbmap',
       'bwa'           :'bwa',
       'bt2'           :'bowtie2',
       'bedtools'      :'bedtools',
       'samtools'      :'samtools',
       'freebayes'     :'freebayes',
       'vcfcombine'    :'vcfcombine',
       'vcfuniq'       :'vcfuniq',
       'vcfstreamsort' :'vcfstreamsort',
       'vcfoverlay'    :'vcfoverlay',
       'vcffirstheader':'vcffirstheader',
       'picard'        :os.environ["PICARD"],
       'gatk'          :os.environ["GATK"],
       'snpeff'        :os.environ["SnpEff"],
       'snpsift'       :os.environ["SnpSift"],
       'mgcr'          :os.environ["meanGenomeCoverage"],
       }


def get_temp_path(file_path):
  '''Get a temporary path based on some other path or directory'''
  
  path_root, file_ext = os.path.splitext(file_path)
  
  return '%s_temp_%s%s' % (path_root, get_rand_string(8), file_ext)


def makedirs(dir_path, exist_ok=False):
  
  try: # Python 3
    os.makedirs(dir_path, exist_ok=exist_ok)
  
  except TypeError:
    if not (os.path.exists(dir_path) and os.path.isdir(dir_path)):
      os.makedirs(dir_path)
      

def report(msg):

  global LOG_FILE_OBJ
  
  if LOGGING:
    if not LOG_FILE_OBJ:
      LOG_FILE_OBJ = open(LOG_FILE_PATH, 'a')
      
    LOG_FILE_OBJ.write(msg)
  
  if not QUIET:
    print(msg)

  
def get_file_ext(file_path):
  
  file_root, file_ext = os.path.splitext(file_path)
   
  if file_ext.lower() in ('.gz','.gzip'):
    file_root, file_ext = os.path.splitext(file_root)
   
  return file_ext
   

def warn(msg, prefix='WARNING'):

  report('%s: %s' % (prefix, msg))

 
def critical(msg, prefix='FAILURE'):

  report('%s: %s' % (prefix, msg))
  sys.exit(0)


def info(msg, prefix='INFO'):

  report('%s: %s' % (prefix, msg))


def check_exe(file_name):
  
  if not os.path.exists(file_name):
    if not locate_exe(file_name):
      critical('Could not locate command exacutable "%s" in system $PATH' % file_name)
 

def locate_exe(file_name):
 
  for path in os.environ["PATH"].split(os.pathsep):
    if os.path.exists(os.path.join(path, file_name)):
      return os.path.join(path, file_name)

  return None
    
    
def call(cmd_args, stdin=None, stdout=None, stderr=None, verbose=True, wait=True, path=None):
  """
  Wrapper for external calls to log and report commands,
  open stdin, stderr and stdout etc.
  """
  
  if verbose:
    info(' '.join(cmd_args))
  
  if path:
    env = dict(os.environ)
    prev = env.get('PATH', '')
    
    if path not in prev.split(':'):
      env['PATH'] = prev + ':' + path
  
  else:
    env = None # Current environment variables 
      
  if stdin and isinstance(stdin, str):
    stdin = open(stdin)

  if stdout and isinstance(stdout, str):
    stdout = open(stdout, 'w')

  if stderr and isinstance(stderr, str):
    stderr = open(stderr, 'a')
  
  if wait:
    subprocess.call(cmd_args, stdin=stdin, stdout=stdout, stderr=stderr, env=env)
      
  else:
    subprocess.Popen(cmd_args, stdin=stdin, stdout=stdout, stderr=stderr, env=env)
  

def open_file(file_path, mode='rU', gzip_exts=('.gz','.gzip')):
  """
  GZIP agnostic file opening
  """
  
  if os.path.splitext(file_path)[1].lower() in gzip_exts:
    file_obj = gzip.open(file_path, mode)
  else:
    file_obj = open(file_path, mode)
  
  return file_obj
 
 
def check_regular_file(file_path):

  msg = ''
  
  if not os.path.exists(file_path):
    msg = 'File "%s" does not exist'
    return False, msg % file_path
  
  if not os.path.isfile(file_path):
    msg = 'Location "%s" is not a regular file'
    return False, msg % file_path
  
  if os.stat(file_path).st_size == 0:
    msg = 'File "%s" is of zero size '
    return False, msg % file_path
    
  if not os.access(file_path, os.R_OK):
    msg = 'File "%s" is not readable'
    return False, msg % file_path
  
  return True, msg

 
def get_rand_string(size):
  
  return ''.join(random.choice(string.ascii_uppercase + string.ascii_lowercase + string.digits) for i in range(size))


def get_safe_file_path(path_name, file_name=None):
   
  if file_name:
    file_path = os.path.join(path_name, file_name)
  else:
    file_path = path_name
 
  if os.path.exists(file_path):
    warn("%s already exists and won't be overwritten..." % file_path)
    
    path_root, file_ext = os.path.splitext(file_path)
    file_path = '%s_%s%s' % (path_root, get_rand_string(8), file_ext)
    
    info('Results will be saved in %s' % file_path)
  
  return file_path


def pair_fastq_files(fastq_paths, pair_tags=('r_1','r_2'), err_msg='Could not pair FASTQ read files.'):
  
  if len(fastq_paths) != len(set(fastq_paths)):
    msg = '%s Repeat file path present.'
    critical(msg % (err_msg))
      
  t1, t2 = pair_tags
  
  paths_1 = []
  paths_2 = []
  
  for path in fastq_paths:
    dir_name, base_name = os.path.split(path)
    
    if (t1 in base_name) and (t2 in base_name):
      msg = '%s Tags "%s" and "%s" are ambiguous in file %s'
      critical(msg % (err_msg, t1, t2, base_name))
    
    elif t1 in base_name:
      paths_1.append((path, dir_name, base_name))
    
    elif t2 in base_name:
      paths_2.append((path, dir_name, base_name))
     
    else:
      msg = '%s File name %s does not contain tag "%s" or "%s"'
      critical(msg % (err_msg, base_name, t1, t2))
  
  n1 = len(paths_1)
  n2 = len(paths_2)
  
  if n1 != n2:
    msg = '%s Number of read 1 (%d) and read 2 (%d) files do not match'
    critical(msg % (err_msg, n1, n2))
  
  fastq_paths_r1 = []
  fastq_paths_r2 = []
  
  for path_1, dir_name_1, file_1 in paths_1:
    seek_file = file_1.replace(t1, t2)
    found = []
    
    for path_2, dir_name_2, file_2 in paths_2:
      if dir_name_1 != dir_name_2:
        continue
    
      if file_2 == seek_file:
        found.append(path_2)
    
    if len(found) == 0:
      # No need to check unpaired read 2 files as these always result in an isolated read 1
      msg = '%s No read 2 file "%s" found to pair with %s'
      critical(msg % (err_msg, seek_file, path_1))
         
    else: 
      # Repeat Read 2 files not possible as repeats checked earlier
      fastq_paths_r1.append(path_1)
      fastq_paths_r2.append(found[0])
  
  return fastq_paths_r1, fastq_paths_r2


def match_files(file_paths, file_pattern):

  from fnmatch import fnmatch
  from os.path import basename
  
  # Like glob, but on a list of strings
  
  return [fp for fp in file_paths if fnmatch(basename(fp), file_pattern)]


def get_bam_chromo_sizes(bam_file_path):
  # Looks in header of BAM file to get chromosome/contig names and thier lengths  
    
  cmd_args = ['samtools', 'idxstats', bam_file_path]
  
  proc = subprocess.Popen(cmd_args, shell=False,
                          stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE)
                          
  std_out_data, std_err_data = proc.communicate()
  chromos_sizes = []
  
  for line in std_out_data.decode('ascii').split('\n'):
    if line:
      ref_name, seq_len, n_mapped, n_unmapped = line.split()
      seq_len = int(seq_len)
      
      if seq_len:
        chromos_sizes.append((ref_name, seq_len))
        
  return chromos_sizes
  

def _parallel_func_wrapper(queue, target_func, proc_data, common_args):
  
  for t, data_item in proc_data:
    result = target_func(data_item, *common_args)
    
    if queue:
      queue.put((t, result))
  

def parallel_split_job(target_func, split_data, common_args, num_cpu=MAX_CORES, collect_output=True):
  
  num_tasks   = len(split_data)
  num_process = min(num_cpu, num_tasks)
  processes   = []
  
  if collect_output:
    queue = multiprocessing.Queue() # Queue will collect parallel process output
  
  else:
    queue = None
    
  for p in range(num_process):
    # Task IDs and data for each task
    # Each process can have multiple tasks if there are more tasks than processes/cpus
    proc_data = [(t, data_item) for t, data_item in enumerate(split_data) if t % num_cpu == p]
    args = (queue, target_func, proc_data, common_args)

    proc = multiprocessing.Process(target=_parallel_func_wrapper, args=args)
    processes.append(proc)
  
  for proc in processes:
    proc.start()
  
  
  if queue:
    results = [None] * num_tasks
    
    for i in range(num_tasks):
      t, result = queue.get() # Asynchronous fetch output: whichever process completes a task first
      results[t] = result
 
    queue.close()
 
    return results
  
  else:
    for proc in processes: # Asynchromous wait and no output captured
      proc.join()

  
for name in EXE:
  check_exe(EXE[name])
