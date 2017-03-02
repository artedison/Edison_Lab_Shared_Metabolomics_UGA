#!/usr/bin/python

# Requires an fid_*.com and ft_*.com template file in the folder "template"
# AND
# Requires a path_*.txt giving the path to the raw data
# will run template/fid.com for each folder in the raw data directory


import os
import glob
import shutil
import re


# Get directory listing, not including certain directories
raw_path_filename = glob.glob('path_*.txt')[0]
with open(raw_path_filename, 'r') as raw_data_path:
    raw_data = raw_data_path.read()
    raw_data = raw_data.replace("\n","")
path_orig = os.getcwd()
os.chdir(raw_data)
dirList = os.listdir(os.getcwd())
dirList = [f for f in dirList if os.path.isdir(f) and f not in ['fid','ft','proc_file','template']]
os.chdir(path_orig)

# Delete directories if they exist
if os.path.exists('fid'):
   shutil.rmtree('fid')
if os.path.exists('ft'):
    shutil.rmtree('ft')
if os.path.exists('proc_file'):
    shutil.rmtree('proc_file')
os.mkdir('ft')
os.mkdir('fid')
os.mkdir('proc_file')

# excecute fid.com for all dirs
os.chdir('template')
for f in dirList:
    exit_status = os.system("fid_*.com " + f)
    if exit_status:
        quit()
os.chdir('..')
fidList=os.listdir('fid')
fidList=[ os.path.splitext(f)[0] for f in fidList] #remove extension

# excecute ft.com for all dirs
os.chdir('template')
for f in fidList:
    exit_status = os.system("ft_*.com " + f)
    if exit_status:
        quit()
os.chdir('..')
#read in the template ft script
#replace in and out with each file's name
ft_name = glob.glob("template/ft_*.com")[0]
with open(ft_name,'r') as fid:
    txt = fid.read()

for f in fidList:
    new_txt = re.sub('\$1', f, txt)
    with open('./proc_file/' + 'ft_'+f+'.com','w') as fid:
        fid.writelines(new_txt)

os.system('chmod +x ./proc_file/*')
    

