#from django.shortcuts import render

from django.shortcuts import render
from main.form import LoadForm


#def welcome(request):
#    return render(request,"welcome.html")
  
def set_crossfil_page(request, num_files=1):
    
  if request.method == "POST":
    form_dict = dict(request.POST)

    if 'addfq' in form_dict:
      num_files = int(form_dict['num_files'][0]) + 1
    
    elif 'reset' in form_dict:
      num_files = 1
    
    form = LoadForm(num_files, request.POST, request.FILES)
    
  else:
    form = LoadForm(num_files)

  my_vars = {'form':form, 'num_files':num_files}

  return render(request,'main/crossfiltemplate.html',my_vars)
