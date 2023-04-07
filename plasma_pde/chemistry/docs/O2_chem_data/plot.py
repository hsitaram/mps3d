import matplotlib.pyplot as plt
import numpy as np
from sys import argv

x=np.array([])
data=np.array([])

fname=argv[1]

infile=open(fname,'r')

for line in infile:
	x    = np.append(x,float(line.split()[0]))
	data = np.append(data,float(line.split()[1]))

npoints=100
Tmin=min(x)
Tmax=max(x)*2

print Tmin,Tmax

datafit=np.zeros(npoints)
Temp=np.zeros(npoints)
dTemp=(Tmax-Tmin)/(float(npoints)-1)

p=np.polyfit(1.0/x,np.log(data),2)
print p

for i in range(npoints):
	Temp[i]=Tmin+i*dTemp
	datafit[i]=np.exp(p[2])*np.exp(p[0]/(Temp[i]**2)+(p[1]/Temp[i]))
	
	#O2+E->O2+ + 2E
	#datafit[i]=1.043e-13*np.exp(-2.74e5/Temp[i] - 2.0015e8/(Temp[i]**2) + 5.33e13/(Temp[i]**3) - 4.001e17/(Temp[i]**4))
	
	#O2+E->2O + E
	#datafit[i]=9.577e-16*np.exp(1.246e5/Temp[i] - 8.647e9/(Temp[i]**2) + 1.381e14/(Temp[i]**3) - 6.953e17/(Temp[i]**4))
	
	#O2+E->O + O-
	#datafit[i]=1.46e-14*np.exp(6.21e4/Temp[i] - 7.27e9/(Temp[i]**2) + 1.25e14/(Temp[i]**3) - 6.57e17/(Temp[i]**4))
	

(fig,ax)=plt.subplots(1,2)
ax[0].plot(1.0/x,np.log(data),label="data")
ax[0].plot(1.0/Temp,np.log(datafit),label="datafit")
ax[0].legend()

ax[1].plot(x,data,label="data")
ax[1].plot(Temp,datafit,label="datafit")
ax[1].set_yscale('log')
ax[1].legend()

plt.tight_layout()
plt.show()
