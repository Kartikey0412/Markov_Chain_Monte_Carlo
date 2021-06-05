# -*- coding: utf-8 -*-


import numpy as np
from pandas import DataFrame, Series
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import scipy
import numpy as np
import pymc as mc


xls = pd.ExcelFile('data worksheet1.xlsx')
df=xls.parse('Sheet1',   index_col=None,na_values=['NA'])
a=df['MOI20-30'].values
b=df['Coculture_Mon'].values

mean=df.mean(axis =0)
mean=mean[1:]
mean=np.array_split(mean,4)
mean=np.array(mean)
meanmon=[]
meanfri=[]
meanclus=[]
for i in range(len(mean)):
    meanmon1=mean[i][0]
    meanmon.append(meanmon1)
    meanfri1=mean[i][1]
    meanfri.append(meanfri1)
    meanclus1=mean[i][2]
    meanclus.append(meanclus1)
    #print mean
    
std=df.std(axis=0)
std=std[1:]
std=np.array_split(std,4)
std=np.array(std)
stdmon=[]
stdfri=[]
stdclus=[]
for i in range(len(std)):
    stdmon1=std[i][0]
    stdmon.append(stdmon1)
    stdfri1=std[i][1]
    stdfri.append(stdfri1)
    stdclus1=std[i][2]
    stdclus.append(stdclus1)
    

t=[0.09,0.05,0.02,0.01]

def sd(t,stdmon,irh,it):
    return (irh)/(1-np.exp(-(t/((it)*0.693))))
    
def fitmon(a,meanmon,t,stdmon):
    irh,it=a
    chisqmon=(meanmon-sd(t,stdmon,irh,it))/(stdmon)

    return chisqmon
    
a0=[0.33,0.01]   
minmon=scipy.optimize.leastsq(fitmon,a0,args=(meanmon,t,stdmon))
print minmon 

z=np.array(df)
z1=np.transpose(z)
s1=scipy.stats.skew(z1[1])
k1=scipy.stats.kurtosis(z1[1])
s4=scipy.stats.skew(z1[4])
k4=scipy.stats.kurtosis(z1[4])
s7=scipy.stats.skew(z1[7])
k7=scipy.stats.kurtosis(z1[7])
s10=scipy.stats.skew(z1[10])
k10=scipy.stats.kurtosis(z1[10])

#skew=np.array(s1,s4,s7,s10)
y=[s1,s4,s7,s10]
skew=np.array(y)
plt.figure(0)      
t1=np.linspace(0.009,0.1,100)
yfit=  ((minmon[0][0])/(1-np.exp(-(t1/(minmon[0][1]*0.693)))))
plt.plot(t,meanmon,"or",label='mean values')
plt.title('hek-cell line')

t2=[]
for i in range (len(z1[1])):
    ta=t[0]
    t2.append(ta)
t3=[]
for i in range (len(z1[4])):
    ta=t[1]
    t3.append(ta)
t4=[]
for i in range (len(z1[7])):
    ta=t[2]
    t4.append(ta)
t5=[]
for i in range (len(z1[10])):
    ta=t[3]
    t5.append(ta)



plt.plot(t2,z1[1],"or",color='k',label='90ms')

plt.plot(t3,z1[4],"or",color='g',label='50ms')

plt.plot(t4,z1[7],"or",color='y',label='20ms')

plt.plot(t5,z1[10],"or",color='m',label='10ms')

plt.errorbar(t, meanmon, yerr=stdmon, fmt='o',color='r')

plt.plot(t1,yfit,label='model')
legend(loc='upper right')
plt.show()

p = pd.read_excel('cookiesdata.xlsx','Sheet1',index_col=None,na_values=['NA'])

mean1= np.mean(p)
mean_co=mean1.values

std1=np.std(p)
sig_co=std1.values
tn=[0.09,0.05,0.03,0.01]


def sdn(tn,std1,irhn,itn):
    return (irhn)/(1-np.exp(-(t/((itn)*0.693))))
    
def fitmonn(an,mean1,tn,std1):
    irhn,itn=an
    chisqmon=(mean1-sdn(tn,std1,irhn,itn))/(std1)

    return chisqmon
    
a0n=[0.33,0.01]   
minmonn=scipy.optimize.leastsq(fitmonn,a0n,args=(mean1,tn,std1))
print(minmonn)

zn=np.array(p)
z1n=np.transpose(zn)
s1n=scipy.stats.skew(z1n[0])
k1n=scipy.stats.kurtosis(z1n[0])
s4n=scipy.stats.skew(z1n[1])
k4n=scipy.stats.kurtosis(z1n[1])
s7n=scipy.stats.skew(z1n[2])
k7n=scipy.stats.kurtosis(z1n[2])
s10n=scipy.stats.skew(z1n[3])
k10n=scipy.stats.kurtosis(z1n[3])

#skew1=np.array(s1n,s4n,s7n,s10n)
y1=(s1n,s4n,s7n,s10n)
skew1=np.array(y1)
plt.figure(1)      
t1n=np.linspace(0.009,0.1,100)
yfitn=  ((minmonn[0][0])/(1-np.exp(-(t1n/(minmonn[0][1]*0.693)))))
plt.plot(tn,mean1,"or",label='mean values')
plt.title('fibroblast')
t2n=[]
for i in range (len(z1n[0])):
    ta=tn[0]
    t2n.append(ta)
t3n=[]
for i in range (len(z1n[1])):
    ta=tn[1]
    t3n.append(ta)
t4n=[]
for i in range (len(z1n[2])):
    ta=tn[2]
    t4n.append(ta)
t5n=[]
for i in range (len(z1n[3])):
    ta=tn[3]
    t5n.append(ta)



plt.plot(t2n,z1n[0],"or",color='k',label='90ms')

plt.plot(t3n,z1n[1],"or",color='g',label='50ms')

plt.plot(t4n,z1n[2],"or",color='y',label='30ms')

plt.plot(t5n,z1n[3],"or",color='m',label='10ms')
plt.errorbar(tn, mean1, yerr=std1, fmt='o',color='r')

plt.plot(t1n,yfitn,label='model')
legend(loc='upper right')
plt.show()
 

#prior distribution of our parameters
Irh1 = mc.Uniform("Irh1",0,0.5) #A uniform prior model is chosen for Irh
tch1 = mc.Uniform("tch1",0,5.0) # A uniform prior model is also chosen for tch

#
stdmon=np.array(stdmon)
meanmon=np.array(meanmon)

# using the deterministic decorator to make our deterministic function
@mc.deterministic(plot=False)  
def model1(Irh1=Irh1,tch1=tch1):   #the deterministic function model1
    out=np.empty(len(stdmon))
    out=(Irh1)/(1-np.exp(-(tn/((tch1)*0.693))))  #out variable-which is for intensity uses the stochastic variables to be computed
     
    return out

D=mc.Normal('D', mu=model1, tau=1/sig_co**2, value=mean_co, observed=True)

m1 = mc.MAP([Irh1,tch1,D])   #using the MAP class for setting the vlaues of the stochastic variables to their maximum posteriori values

m1.fit()
#running the Morkov Chain Monte Carlo algorithm
m1=mc.MCMC(m1)
#giving 50,000 iterations
m1.sample(50000)
#summary to see the results
m1.summary()
#plot is made Matplot is conjuction with PyMC
mc.Matplot.plot(m1)




#getting the mean value of Irh and tch and making a fit
psfit1=(Irh1.stats()['mean'])/(1-np.exp(-(t/((tch1.stats()['mean'])*0.693))))
plt.figure()
plt.plot(t,psfit1)
plt.show()