# -*- coding: utf-8 -*-
# Running MCMC simulation with using Weibull Distribution of final Stochastic Variable
# Finding trace of difference in rheobases of 2 cell lines

import numpy as np
from pandas import DataFrame, Series
import pandas as pd
import matplotlib.pyplot as plt
from lmfit import minimize, Parameters
import scipy
import numpy as np
import pymc as mc

xls = pd.ExcelFile('data worksheet1.xlsx')
df = xls.parse('Sheet1', index_col=None, na_values=['NA'])
a = df['MOI20-30'].values
b = df['Coculture_Mon'].values

c = np.array(df)
d = c.T
d_1 = d[1:]

data_temp = (d_1[0], d_1[3], d_1[6], d_1[9])
data = np.array(data_temp)

mean = df.mean(axis=0)
mean_z = mean[1:]
mean = np.array_split(mean_z, 4)
mean = np.array(mean)
meanmon = []
meanfri = []
meanclus = []
for i in range(len(mean)):
    meanmon1 = mean[i][0]
    meanmon.append(meanmon1)
    meanfri1 = mean[i][1]
    meanfri.append(meanfri1)
    meanclus1 = mean[i][2]
    meanclus.append(meanclus1)
    # print mean

std = df.std(axis=0)
std = std[1:]
std = np.array_split(std, 4)
std = np.array(std)
stdmon = []
stdfri = []
stdclus = []

for i in range(len(std)):
    stdmon1 = std[i][0]
    stdmon.append(stdmon1)
    stdfri1 = std[i][1]
    stdfri.append(stdfri1)
    stdclus1 = std[i][2]
    stdclus.append(stdclus1)

t = [0.09, 0.05, 0.02, 0.01]

z = np.array(df)
z1 = np.transpose(z)
s1 = scipy.stats.skew(z1[1])
k1 = scipy.stats.kurtosis(z1[1])
s4 = scipy.stats.skew(z1[4])
k4 = scipy.stats.kurtosis(z1[4])
s7 = scipy.stats.skew(z1[7])
k7 = scipy.stats.kurtosis(z1[7])
s10 = scipy.stats.skew(z1[10])
k10 = scipy.stats.kurtosis(z1[10])

# skew=np.array(s1,s4,s7,s10)
y = [s1, s4, s7, s10]
skew = np.array(y)
# plt.figure(0)
t1 = np.linspace(0.009, 0.1, 100)

t2 = []
for i in range(len(z1[1])):
    ta = t[0]
    t2.append(ta)
t3 = []
for i in range(len(z1[4])):
    ta = t[1]
    t3.append(ta)
t4 = []
for i in range(len(z1[7])):
    ta = t[2]
    t4.append(ta)
t5 = []
for i in range(len(z1[10])):
    ta = t[3]
    t5.append(ta)

p = pd.read_excel('cookiesdata.xlsx', 'Sheet1', index_col=None, na_values=['NA'])

tn = [0.09, 0.05, 0.03, 0.01]
q = np.array(p)
q1 = q.T
n1 = []
for i in range(len(q1)):
    #    q2=0
    q2 = q1[i]

    # q2=q1[3]
    q2s = np.sort(q2)

    z = plt.hist(q2s, len(q2s), normed=True)

    z1 = z[0]  # y-values
    z2 = z[1]  # x-values
    k = 1.5


    # plt.plot(z)
    def wb(q2s, lam):
        #    return ((k/lam)*((q2s/lam)**(k-1)))*(np.exp(-(q2s/lam)**k))
        min = ((k / lam) * ((q2s / lam) ** (k - 1))) * (np.exp(-(q2s / lam) ** k))
        #    min=((a[0]/a[1])*((q2s/a[1])**(a[0]-1)))*(np.exp(-((q2s/a[1])**a[0])))
        return min


    def fitmon(a, z1, q2s):
        lam = a
        chisq = (z1 - wb(q2s, lam)) / (0.1)
        #    min=np.sqrt(z1-wb(q2s,k,lam))
        #
        return chisq


    # k=1.5

    a = [1]

    n = scipy.optimize.leastsq(fitmon, a, args=(z1, q2s))
    n1 = np.append(n[0][0], n1)

n1c = []
for i in range(len(data)):
    #    q2=0
    q2c = data[i]

    # q2=q1[3]
    q2sc = np.sort(q2c)

    zc = plt.hist(q2sc, len(q2sc), normed=True)

    z1c = zc[0]  # y-values
    z2c = zc[1]  # x-values
    kc = 1.5


    # plt.plot(z)
    def wbc(q2sc, lamc):
        #    return ((k/lam)*((q2s/lam)**(k-1)))*(np.exp(-(q2s/lam)**k))
        min = ((kc / lamc) * ((q2sc / lamc) ** (kc - 1))) * (np.exp(-(q2sc / lamc) ** kc))
        #    min=((a[0]/a[1])*((q2s/a[1])**(a[0]-1)))*(np.exp(-((q2s/a[1])**a[0])))
        return min


    def fitmonc(ac, z1c, q2sc):
        lamc = ac
        chisqc = (z1c - wbc(q2sc, lamc)) / (0.1)
        #    min=np.sqrt(z1-wb(q2s,k,lam))
        #
        return chisqc


    # k=1.5

    ac = [1]

    nc = scipy.optimize.leastsq(fitmonc, ac, args=(z1c, q2sc))
    n1c = np.append(nc[0][0], n1c)

# BAYESIAN FITTING

Irh_co = mc.Uniform('Irh_co', 0.02, 0.5)
tch_co = mc.Uniform('tch_co', 0.02, 5.0)
a = 1.5
# b = mc.Uniform('b', lower=0, upper=10, value=5, doc='Weibull beta parameter')

Irh_ch = mc.Uniform('Irh_ch', 0.02, 5.0)
tch_ch = mc.Uniform('tch_ch', 0.02, 5.0)
a1 = 1.5


@mc.deterministic(plot=False)
#
def modelco(Irh_co=Irh_co, tch_co=tch_co):
    out = np.empty(len(stdmon))
    out = (Irh_co) / (1 - np.exp(-(tn / ((tch_co) * 0.693))))
    #    out1=np.log(out)
    #
    ##    a1=np.exp(out+(stdmon**2)/2)
    #
    return out


#
@mc.deterministic(plot=False)
def modelch(Irh_ch=Irh_ch, tch_ch=tch_ch):
    out = np.empty(len(stdmon))
    out = (Irh_ch) / (1 - np.exp(-(t / ((tch_ch) * 0.693))))
    #    out1=np.log(out)

    #    a1=np.exp(out+(stdmon**2)/2)

    return out


@mc.deterministic(plot=False)
def diff_irh(Irh_ch=Irh_ch, Irh_co=Irh_co):
    return Irh_ch - Irh_co


@mc.deterministic(plot=False)
def diff_tch(tch_ch=tch_ch, tch_co=tch_co):
    return tch_ch - tch_co


Dco = mc.Weibull('Dco', alpha=a, beta=modelco / 0.48, value=n1, observed=True)
Dch = mc.Weibull('Dch', alpha=a1, beta=modelch / 0.48, value=n1c, observed=True)
#
#
m1 = mc.MAP([Irh_co, tch_co, Dco])
m2 = mc.MAP([Irh_ch, tch_ch, Dch])
m3 = mc.MAP([diff_irh, Dco, Dch])

m1.fit()
m2.fit()
m3.fit()
##m4.fit()
##mz.fit()
m1 = mc.MCMC(m1)
m2 = mc.MCMC(m2)
m3 = mc.MCMC(m3)
##m4=mc.MCMC(m4)
##mz=mc.MCMC(mz)
m1.sample(50000)
m2.sample(50000)
m3.sample(50000)
##m4.sample(50000)
##mz.sample(50000)
m1.summary()
m2.summary()
m3.summary()
##m4.summary()
##mz.summary()
#
mc.Matplot.plot(m1)
mc.Matplot.plot(m2)
mc.Matplot.plot(m3)

psfit_co = (Irh_co.stats()['mean']) / (1 - np.exp(-(tn / ((tch_co.stats()['mean']) * 0.693))))
psfit_ch = (Irh_ch.stats()['mean']) / (1 - np.exp(-(tn / ((tch_ch.stats()['mean']) * 0.693))))

plt.figure()
plt.plot(tn, psfit_co)
plt.figure()
plt.plot(t, psfit_ch)
plt.show()
y = m3.trace("diff_irh")[:]
y1 = m1.trace("Irh_co")[:]
y2 = m2.trace("Irh_ch")[:]
y3 = ((y2 - y1) / y2) * 1000
binwidth = 0.005
plt.figure()
plt.hist(y3, histtype='stepfilled', bins=25, alpha=0.85, normed=True)
plt.hist(y)
plt.show()