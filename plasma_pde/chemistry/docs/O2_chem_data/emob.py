import matplotlib.pyplot as plt
import numpy as np
from sys import argv

x=np.array([])
data=np.array([])

fname=argv[1]

infile=open(fname,'r')

for line in infile:
	x    = np.append(x,float(line.split()[1]))
	data = np.append(data,float(line.split()[2]))

npoints=100

#convert x to electron temperature, scale by 1.5 (3/2)
x=x/1.5

#convert data to electron mobility
e=1.602176e-19
m=9.109382e-31
kB=1.380650e-23
eVtoK=11604.5
pi=3.14159265
A2tom2=1e-20

cbar=np.sqrt(8.0*kB*x*eVtoK/pi/m)

#print cbar

mob_times_N=(e/m)/cbar
#print mob_times_N
mob_times_N=mob_times_N/(data*A2tom2)
#print mob_times_N

Tmin=min(x)
Tmax=max(x)*4

#print Tmin,Tmax

datafit=np.zeros(npoints)
Temp=np.zeros(npoints)
dTemp=(Tmax-Tmin)/(float(npoints)-1)

p=np.polyfit(np.log(x),np.log(mob_times_N),2)
print p
print "test value at 20 eV",np.exp(p[0]*np.log(20.)**2 + p[1]*np.log(20.) + p[2])

print np.exp(0.061*2.995**2 - 0.342*2.995 + 55.89)

for i in range(npoints):
	Temp[i]=Tmin+i*dTemp
	#datafit[i]=np.exp(p[0]*np.log(Temp[i])**4 + p[1]*np.log(Temp[i])**3 + p[2]*np.log(Temp[i])**2 + p[3]*np.log(Temp[i]) + p[4])
	#datafit[i]=np.exp(p[0]*np.log(Temp[i])**3+p[1]*np.log(Temp[i])**2+p[2]*np.log(Temp[i])+p[3])
	datafit[i]=np.exp(p[0]*np.log(Temp[i])**2 + p[1]*np.log(Temp[i]) + p[2])
	


(fig,ax)=plt.subplots(1,2)
ax[0].plot(np.log(x),np.log(mob_times_N),marker='o',label="data")
ax[0].plot(np.log(Temp),np.log(datafit),label="datafit")
ax[0].legend()

ax[1].plot(x,mob_times_N,marker='o',label="data")
ax[1].plot(Temp,datafit,label="datafit")
ax[1].set_yscale('log')
ax[1].legend()

plt.tight_layout()
plt.show()
