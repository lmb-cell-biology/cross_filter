#!/usr/bin/python

import os

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


