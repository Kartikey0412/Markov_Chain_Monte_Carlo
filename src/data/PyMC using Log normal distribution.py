# -*- coding: utf-8 -*-

import numpy as np
from pandas import DataFrame, Series
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import scipy
import numpy as np
import pymc as mc
import statsmodels.sandbox.distributions.extras as extras

xls = pd.ExcelFile('data worksheet1.xlsx')
df=xls.parse('Sheet1',   index_col=None,na_values=['NA'])
a=df['MOI20-30'].values
b=df['Coculture_Mon'].values

c=np.array(df)
d=c.T
d_1=d[1:]
d_2=np.log(d_1)
d_2_0=np.mean(d_2[0])
d_2_1=np.mean(d_2[3])
d_2_2=np.mean(d_2[6])
d_2_3=np.mean(d_2[9])

me=(d_2_0,d_2_1,d_2_2,d_2_3)
mean_ch=np.array(me)


st_0=np.std(d_2[0])
st_1=np.std(d_2[3])
st_2=np.std(d_2[6])
st_3=np.std(d_2[9])

st=(st_0,st_1,st_2,st_3)
std_ch=np.array(st)

mean=df.mean(axis =0)
mean_z=mean[1:]
mean=np.array_split(mean_z,4)
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


y=[s1,s4,s7,s10]
skew=np.array(y)
#plt.figure(0)      
t1=np.linspace(0.009,0.1,100)
yfit=  ((minmon[0][0])/(1-np.exp(-(t1/(minmon[0][1]*0.693)))))


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





p = pd.read_excel('cookiesdata.xlsx','Sheet1',index_col=None,na_values=['NA'])
q=np.array(p)
q1=np.log(q)
q2=q1.T



mean0=np.mean(q2[0])
meana=np.mean(q2[1])
mean2=np.mean(q2[2])
mean3=np.mean(q2[3])

mean_1=(mean0,meana,mean2,mean3)
mean_co=np.array(mean_1)




std0=np.std(q2[0])
stda=np.std(q2[1])
std2=np.std(q2[2])
std3=np.std(q2[3])

std_1=(std0,stda,std2,std3)
std_co=np.array(std_1)


tn=[0.09,0.05,0.03,0.01]


#
#BAYESIAN FITTING

Irh_co = mc.Uniform("Irh_co",0.02,0.5)
tch_co = mc.Uniform("tch_co",0.02,5.0)

Irh_ch = mc.Uniform("Irh_ch",0.02,0.5)
tch_ch = mc.Uniform("tch_ch",0.02,5.0)

stdmon=np.array(stdmon)
meanmon=np.array(meanmon)

@mc.deterministic(plot=False)

def modelco(Irh_co=Irh_co,tch_co=tch_co):
    out=np.empty(len(stdmon))
    out=(Irh_co)/(1-np.exp(-(tn/((tch_co)*0.693))))
    out1=np.log(out)
    
#    a1=np.exp(out+(stdmon**2)/2)
     
    return out1
    
@mc.deterministic(plot=False)

def modelch(Irh_ch=Irh_ch,tch_ch=tch_ch):
    out=np.empty(len(stdmon))
    out=(Irh_ch)/(1-np.exp(-(t/((tch_ch)*0.693))))
    out1=np.log(out)
     
    return out1
    
@mc.deterministic(plot=False)
def diff_irh(Irh_ch=Irh_ch,Irh_co=Irh_co):     
    return Irh_ch-Irh_co
    
@mc.deterministic(plot=False)
def diff_tch(tch_ch=tch_ch,tch_co=tch_co):     
    return tch_ch-tch_co
    
#D=mc.Normal('D', mu=model1, tau=1/(np.exp(2*meanmon+stdmon**2)*(np.exp(stdmon**2)-1)), value=np.exp(meanmon+(stdmon**2)/2), observed=True)
Dco=mc.Normal('D', mu=modelco, tau=1/std_co**2, value=mean_co, observed=True)
Dch=mc.Normal('D', mu=modelch, tau=1/std_ch**2, value=mean_ch, observed=True)

m1 = mc.MAP([Irh_co,tch_co,Dco])
m2 = mc.MAP([Irh_ch,tch_ch,Dch])
m3 = mc.MAP([diff_irh,Dco,Dch])

m1.fit()
m2.fit()
m3.fit()
#m4.fit()
#mz.fit()
m1=mc.MCMC(m1)
m2=mc.MCMC(m2)
m3=mc.MCMC(m3)
#m4=mc.MCMC(m4)
#mz=mc.MCMC(mz)
m1.sample(50000)
m2.sample(50000)
m3.sample(50000)
#m4.sample(50000)
#mz.sample(50000)
m1.summary()
m2.summary()
m3.summary()
#m4.summary()
#mz.summary()

mc.Matplot.plot(m1)
mc.Matplot.plot(m2)
mc.Matplot.plot(m3)
#mc.Matplot.plot(m4)
#mc.Matplot.plot(mz)

#psfit=np.log(Irh.stats()['mean'])-np.log((1-np.exp(-(t/((tch.stats()['mean'])*0.693)))))
psfit_co=(Irh_co.stats()['mean'])/(1-np.exp(-(tn/((tch_co.stats()['mean'])*0.693))))
psfit_ch=(Irh_ch.stats()['mean'])/(1-np.exp(-(tn/((tch_ch.stats()['mean'])*0.693))))

plt.figure()
plt.plot(tn,psfit_co)
plt.figure()
plt.plot(t,psfit_ch)
#plt.show()
y=m3.trace("diff_irh")[:]
y1=m1.trace("Irh_co")[:]
y2=m2.trace("Irh_ch")[:]
y3=((y2-y1)/y2)*1000
binwidth=0.005
plt.figure()
plt.hist(y3,histtype='stepfilled',bins=25,alpha=0.85,normed=True)
#plt.hist(y)
plt.show()