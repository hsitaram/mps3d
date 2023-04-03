import matplotlib.pyplot as plt
import numpy as np
from sys import argv
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

fname=argv[1]

(specnames,enum,nbgspecies,nions,nneut,gapd)=getspecinfo()

(x,eden,iden,nden,pot,efield,etemp,ejheat,elcol,inelcol,electroncurr,ioncurr)=readfile(fname,enum,nbgspecies,nions,nneut,len(specnames))
x=x/gapd

xlen=6.0
ylen=3.0
me = 5
m=2
n=3
lt=3


npoints=len(x)

(fig,ax)=plt.subplots(m,n,figsize=(n*xlen,m*ylen))
#fig.tight_layout()

i=0
j=0
ax[i][j].set_title("Electron/ion densities")
ax[i][j].plot(x,eden,'k*-',markevery=me,markeredgecolor='black',linewidth=lt,label="E")
offset=nbgspecies+1
for n in range(nions):
    ax[i][j].plot(x,iden[:,n],markevery=me+7,linewidth=lt,label=specnames[offset+n])
ax[i][j].set_ylabel("number density (#/m3)")
lg=ax[i][j].legend(loc="best")
lg.draw_frame(False)

i=0
j=1
ax[i][j].set_title("Potential/Electric fields")
ax[i][j].plot(x,pot,'r*-',markevery=me,markeredgecolor='red',linewidth=lt,label="Potential")
ax[i][j].set_ylabel("Voltage (V)")
ax_twin=ax[i][j].twinx()
ax_twin.plot(x,efield,'bo-',markevery=me+7,markeredgecolor='blue',linewidth=lt,label="Electric field")
ax_twin.set_ylabel("Electric field (V/m)")
lg=ax[i][j].legend(loc=1)
lg.draw_frame(False)
lg=ax_twin.legend(loc=2)
lg.draw_frame(False)

i=1
j=0
ax[i][j].set_title("Potential/Electron temperature")
ax[i][j].plot(x,pot,'r*-',markevery=me,markeredgecolor='red',linewidth=lt,label="Potential")
ax[i][j].set_ylabel("Voltage (V)")
ax[i][j].set_xlabel("gap distance")
ax_twin=ax[i][j].twinx()
ax_twin.plot(x,etemp,'bo-',markevery=me+7,markeredgecolor='blue',linewidth=lt,label="Electron temperature")
ax_twin.set_ylabel("electron temperature (eV)")
lg=ax[i][j].legend(loc=1)
lg.draw_frame(False)
lg=ax_twin.legend(loc=2)
lg.draw_frame(False)

i=1
j=1
ax[i][j].set_title("Electron source terms")
ax[i][j].plot(x,ejheat,'r*-',markevery=me,markeredgecolor='red',linewidth=lt,label="Joule heating")
ax[i][j].plot(x,elcol,'bo-',markevery=me,markeredgecolor='blue',linewidth=lt,label="Elastic collisions")
ax[i][j].plot(x,inelcol,'gs-',markevery=me,markeredgecolor='green',linewidth=lt,label="Inelastic collisions")
ax[i][j].set_xlabel("gap distance")
ax[i][j].set_ylabel("power density (W/m3)")
lg=ax[i][j].legend(loc="best")
lg.draw_frame(False)

i=0
j=2
ax[i][j].set_title("Current densities")
ax[i][j].plot(x,electroncurr,'r*-',markevery=me,markeredgecolor='red',linewidth=lt,label="electron current")
ax[i][j].plot(x,ioncurr,'bo-',markevery=me,markeredgecolor='blue',linewidth=lt,label="Ion current")
ax[i][j].plot(x,ioncurr+electroncurr,'ks-',markevery=me,markeredgecolor='black',linewidth=lt,label="Total current")
ax[i][j].set_xlabel("gap distance")
ax[i][j].set_ylabel("current density (A/m2)")
lg=ax[i][j].legend(loc="best")
lg.draw_frame(False)

i=1
j=2
ax[i][j].set_title("Neutral densities")
if(nneut > 0):
    offset=nbgspecies+1+nions
    for n in range(nneut):
        ax[i][j].plot(x,nden[:,n],markevery=me+7,linewidth=lt,label=specnames[offset+n])
    ax[i][j].set_ylabel("number density (#/m3)")
    lg=ax[i][j].legend(loc="best")
    lg.draw_frame(False)
    print(np.trapz(ndens[0],x))

fig.tight_layout()
#plt.savefig("plasmaparams.pdf")
plt.show()
