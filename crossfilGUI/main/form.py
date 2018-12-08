from django import forms
  
class LoadForm(forms.Form):
  csv = forms.FileField(required=False)
  
  def __init__(self, num_files, *args, **kw):
    forms.Form.__init__(self, args, kw)
  
    for i in range(num_files):
      self.fields['upload_fastq%d' % (i+1)] = forms.FileField(required=False)
  
  genome_build = forms.FileField(required=False)

