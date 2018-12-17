from django import forms

#class FileFieldForm(forms.Form):
#    file_field = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True}))
    
  
class LoadForm(forms.Form):
  csv = forms.FileField(required=False)

  def __init__(self, num_files, *args, **kw):
    forms.Form.__init__(self, *args, **kw)
  
    for i in range(num_files):
      self.fields['upload_fastq%d' % (i+1)] = forms.FileField(widget=forms.ClearableFileInput(attrs={'multiple': True}),required=False)
#      self.fields['upload_fastq%d' % (i+1)] = FileFieldForm()
  
  
  genome_fasta = forms.FileField(required=False)
  
  genome_gff3 = forms.FileField(required=False)
  
  snpeff_gv = forms.CharField(label="SnpEff genome version",max_length = 100, required=False)
  
  snpeff_ivl = forms.IntegerField(label="SnpEff interval length",required=False,min_value=0)

  pe = forms.BooleanField(label="Paired-end data", required=False)

  pe_tags = forms.CharField(label="Paired-end tags",max_length = 100, required=False)
  
  al = (("bt2","Bowtie2"),("bwa","BWA"),("bbmap","Bbmap"))
  aligner = forms.ChoiceField(choices=al,required=False) 
  
  bt2_index = forms.FileField(required=False,label="Bowtie2 index")
  
  bbmap_index = forms.FileField(required=False,label="Bbmap index")
  
  vc = (("fb","FreeBayes"),("gatk","GATK"))
  variant_caller = forms.ChoiceField(choices=vc,required=False) 
  
  heterozygous = forms.BooleanField(label="Include heterozygous variants", required=False)
  
  min_var_count = forms.IntegerField(label="Minimum number of common variants",required=False,min_value=0)
  
  ref_allele = forms.CharField(label="Reference allele",max_length = 100, required=False)
  
  q = forms.BooleanField(label="Silence progress", required=False)
  
  log = forms.BooleanField(label="Log progress", required=False)
  
  cpu = forms.IntegerField(label="Number of CPUs",required=False,min_value=1)
  
  outdir = forms.FileField(required=False,label="Output directory")
  
