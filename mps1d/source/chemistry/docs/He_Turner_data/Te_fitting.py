import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from sys import argv
import sys
from io import StringIO

def Te_fit_arrh1(x,A,B,C):
    f = A*np.log(x) + B + C/x
    return(f)

def Te_fit_arrh2(x,A,B,C):
    f = B*(x**A)*np.exp(-C/x)
    return(f)

def is_number(word):
    try:
        float(word)
    except:
        return(False)
    return(True)

def is_line_of_numbers(string):
    splt=string.split()
    for i in range(len(splt)):
        if(is_number(splt[i])==False):
            return(False)
    return(True)


bolsigoutfile=argv[1]
search=argv[2]
infile=open(bolsigoutfile,'r')
lines=infile.readlines()
nlines=len(lines)

for l,lineno in  zip(lines,range(nlines)):
    if(search in l):
        break

all_lines=""
for i in range(lineno,nlines):
    line=lines[i]
    splt=line.split()
    if(len(splt)>0):
        if(is_line_of_numbers(line)):
            all_lines=all_lines+line
            #print(np.loadtxt(StringIO(line)))
    else:
        break

arr=np.loadtxt(StringIO(all_lines))
print(arr.shape)

minval=np.min(arr[np.nonzero(arr[:,1]),1])
for i in range(len(arr[:,1])):
    if(arr[i,1]==0.0):
        arr[i,1]=minval

avgval=np.mean(arr[:,1])

#parameters, cov = curve_fit(Te_fit_arrh1,arr[:,0],np.log(arr[:,1]/avgval))
parameters, cov = curve_fit(Te_fit_arrh2,arr[:,0],arr[:,1])
A,B,C=parameters
fitfunc=np.exp(Te_fit_arrh2(arr[:,0],A,B,C))
print(parameters)
np.savetxt("fitting_"+search.replace(" ","")+".dat",np.transpose(np.vstack((np.transpose(arr),fitfunc))))
