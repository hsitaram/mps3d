import yt
from sys import argv
import matplotlib.pyplot as plt
import numpy as np
import glob

#===================================================================
def getspecinfo():

    str1="#NUMBER DENSITIES"
    infile=open('plasma_inputs','r')

    all_lines=infile.readlines()

    gapd=float(all_lines[1].split()[1])

    specstartlinenum=0
    for line in all_lines:
        specstartlinenum+=1
        if(line.strip() == str1):
            break

    specnames=[]

    nspecnum=0
    for l in range(specstartlinenum,len(all_lines)):
        if(all_lines[l].strip()=="#PROBLEM SPECIFIC PARAMETERS"):
            break
        nspecnum+=1

    for l in range(specstartlinenum,specstartlinenum+nspecnum):
        specnames.append(all_lines[l].split()[0])
    

    nbgspc=0
    for i in range(len(specnames)):
        if(specnames[i] == 'E'):
            break
    nbgspc=i
    especnum=i+1

    nions=0
    for i in range(len(specnames)):
        if(specnames[i] == 'E'):
            especnum=i
        if(len(specnames[i].split('+')) > 1):
            nions=nions+1
        if(len(specnames[i].split('-')) > 1):
            nions=nions+1

    infile.close()

    nneutrals=len(specnames)-nions-nbgspc-1

    return(specnames,especnum,nbgspc,nions,nneutrals,gapd)
#===================================================================
def readfile(filename,enum,nbgspecies,nions,nneutrals,nspecies):

    infile=open(filename,'r')

    arr=np.loadtxt(filename)

    offset=0
    x        = arr[:,offset]
    offset+=1
    potential= arr[:,offset]
    offset+=1
    efield   = arr[:,offset]
    offset+=1

    nden     = arr[:,offset:offset+nbgspecies]
    offset += nbgspecies
    elecden  = arr[:,offset]
    offset+=1
    ionden    = arr[:,offset:offset+nions]
    offset+= nions
    elecenrg = arr[:,offset]
    offset+=1
    electemp  = arr[:,offset]
    offset+=1
    jheating = arr[:,offset]
    offset+=1
    elasticol = arr[:,offset]
    offset+=1
    inelasticol = arr[:,offset]
    offset+=1
    electroncurr = arr[:,offset]
    offset+=1
    ioncurr = arr[:,offset]

    return(x,elecden,ionden,nden,potential,efield,electemp,jheating,elasticol,inelasticol,electroncurr,ioncurr)
#===================================================================

fn_pattern= argv[1]
fn_list=[]
try:
    fn_list = sorted(glob.glob(fn_pattern), key=lambda f: int(f.split("soln_")[1].split('.')[0]))
except:
    if(fn_list==[]):
        print("using file of plotfiles..")
        infile=open(argv[1],'r')
        for line in infile:
            fn_list.append(line.split()[0])
        infile.close()

print(fn_list)
first_res=1.0
resarr=np.zeros((len(fn_list)-1,3))
for i in range(1,len(fn_list)):
    (specnames,enum,nbgspecies,nions,nneut,gapd)=getspecinfo()
    (x,eden,iden1,nden,pot,efield,etemp,ejheat,elcol,inelcol,electroncurr,ioncurr)=readfile(fn_list[i-1],enum,nbgspecies,nions,nneut,len(specnames))
    (x,eden,iden2,nden,pot,efield,etemp,ejheat,elcol,inelcol,electroncurr,ioncurr)=readfile(fn_list[i],enum,nbgspecies,nions,nneut,len(specnames))
    residual=np.abs(np.max(iden1)-np.max(iden2))
    if(i==1):
        first_res=residual
    print("abs/rel residual:",residual,residual/first_res)
    resarr[i-1,0]=i
    resarr[i-1,1]=residual
    resarr[i-1,2]=residual/first_res


np.savetxt("resnorms",resarr,delimiter="  ")
