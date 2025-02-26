"""
Author: Lydia A. Kanari-Naish, Imperial College London
Last update: September 2023

EXAMPLE 2: Photon-subtracted/added two-mode squeezed vacuum state


Imports functions from NPyT
Code for optimizing over all suitable NPT criteria for example of a 
photon-subtracted/added two-mode squeezed vacuum state

"""

"Import packages"
from itertools import cycle
import numpy as np
from qutip import *
from scipy.special import comb as comb
import math
from matplotlib import pyplot as plt
import matplotlib.ticker as mtick
from matplotlib.ticker import FormatStrFormatter
import matplotlib.style as style 
from matplotlib.ticker import (AutoMinorLocator, MultipleLocator)
style.use('seaborn-v0_8-colorblind')
from win32com.client import Dispatch
speak = Dispatch("SAPI.SpVoice").Speak
from datetime import datetime
start = datetime.now()


"Import the functions from the NPyT"
from NPyT import *



"Define the photon-subtracted/added TMSV"


def add(fock_dims, z, k, l):
    a1 = tensor(destroy(fock_dims),qeye(fock_dims))
    a2 = tensor(qeye(fock_dims),destroy(fock_dims))

    ground_state = tensor(basis(fock_dims, 0),basis(fock_dims,0))
    TMSV = squeezing(a1, a2, z)*ground_state
    #phonon added state
    add_state = ((a1.dag())**k)*((a2.dag())**l)*TMSV
    return add_state.unit()



def sub(fock_dims, z, k, l):
    a1 = tensor(destroy(fock_dims),qeye(fock_dims))
    a2 = tensor(qeye(fock_dims),destroy(fock_dims))

    ground_state = tensor(basis(fock_dims, 0),basis(fock_dims,0))
    TMSV = squeezing(a1, a2, z)*ground_state
    #phonon subtracted state
    subbed_state = ((a1)**k)*((a2)**l)*TMSV
    return subbed_state.unit()



add_vec=np.frompyfunc(add,4,1)



sub_vec=np.frompyfunc(sub,4,1)




# Squeezing parameter
zs = np.linspace(0.0001,2,10)
#zs = np.linspace(0.0001,2,100)
sub_states = sub_vec(N_fock,zs,1,1)



eta = 0.8
N=0
M=1000



"Preliminary search of successful determinants"
vals_22,combs_22=my_state(2,2,N_fock,np.array([sub_states],dtype=object))
#Note EIV and EV are not unique as noted in the manuscript so let's removes these from combs_22
combs_22=np.array([[ 2,  4],[ 6, 11],[ 7, 12],[ 7, 14]])




## This generates the data for value of determinant plus error bars. 
# as a function of ks
# fixed M, t, kappa, N
data_y=[]
data_err=[]
for ls in combs_22:
    sh = len(zs)
    y=np.zeros((sh))
    err=np.zeros((sh))
    for i in range(len(zs)):
        #y[i]=(TD_det(20,sub_states[i],2,ls,eta,N))
        y[i]=(TD_det(N_fock,sub_states[i],2,ls,eta,N))
        #err[i]=(error(20,sub_states[i],2,ls,eta,N,M))
        err[i]=(error(N_fock,sub_states[i],2,ls,eta,N,M))
 
    data_y.append(y)
    data_err.append(err)
data_y=np.array([data_y])
data_err=np.array([data_err])




"FIG 1 PHOTON SUBTRACTED TWO-MODE SQUEEZED STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

#plt.plot(zs,np.zeros(len(zs)), 'black', alpha = 0.5)

for i in range(len(combs_22)):
    y=data_y[0][i]
    err=data_err[0][i]
    c1 = ax.plot(zs, y, '-', label = str(combs_22[i]))
    c2 = ax.fill_between(zs, y-err, y+err,alpha=0.2)
    

plt.xlabel(r'$\zeta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Determinant',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(loc=(0.2,0.2))
plt.title(r'Determinant for M={}, $\eta$={}, N={}'.format(M,eta,N))


ax.set_xlim([zs[0],zs[-1]])

# symmetric log scale for both positive and negative y values
plt.yscale("symlog")

plt.grid()

plt.savefig('det_error_bars_sub.svg')

plt.show()


"FIG 2 PHOTON SUBTRACTED TWO-MODE SQUEEZED STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(zs,0.95*np.ones(len(zs)),'k--',label='95%')

ys = data_y
errs = data_err

for i in range(len(combs_22)):   
    c1 = ax.plot(zs, confidence_level(ys[0][i],errs[0][i]),label=str(combs_22[i]))


plt.xlabel(r'$\zeta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.xticks(fontsize=12)
fig.legend(loc=(0.4,0.2))
plt.title(r'M={}, $\eta$={}, N={}'.format(M,eta,N))
ax.set_xlim([zs[0],zs[-1]])
ax.set_ylim([0.5,1.01])
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

plt.grid()

plt.savefig('confidence_sub.svg')

plt.show()




# Data for how confidence depends on measurement number
Ms = [10**i for i in np.linspace(2,4,10)]
#Ms = [10**i for i in np.linspace(2,4,100)]
eta = 0.8
zeta = 1
N = 0



sub_state=sub_vec(N_fock,zeta,1,1)



## This generates the data for value of determinant plus error bars. 
## for 2x2 matrices
# as a function of M
# fixed alpha, t, kappa, N
data_ms_y=[]
data_ms_err=[]
for ls in combs_22:
    L = len(Ms)
    y=np.zeros((L))
    err=np.zeros((L))
    for i in range(L):
        #y[i]=TD_det(20,sub_state,2,ls,eta,N)
        y[i]=TD_det(N_fock,sub_state,2,ls,eta,N)
        #err[i]=error(20,sub_state,2,ls,eta,N,Ms[i])
        err[i]=error(N_fock,sub_state,2,ls,eta,N,Ms[i])
        #print(y,err)
    data_ms_y.append(y)
    data_ms_err.append(err)



"FIG 3 PHOTON SUBTRACTED TWO-MODE SQUEEZED STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(Ms,0.95*np.ones(len(Ms)),'k--',label='95%')

ys = data_ms_y
errs = data_ms_err

for i in range(len(combs_22)):
    c1 = ax.plot(Ms, confidence_level(ys[i],errs[i]),label=str(combs_22[i]))

plt.xlabel(r'$M_{tot}$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.xticks(fontsize=12)
fig.legend(loc=(0.7,0.2))
plt.title(r'$\zeta$={}, $\eta$={}, N={}'.format(zeta,eta,N))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xscale('log')

plt.grid(which='both')
ax.set_xlim([Ms[0],Ms[-1]])
plt.ylim(0.62,1.02)

plt.show()

plt.savefig('confidence_ms_sub.svg')



# Data for how confidence depends on measurement number
etas = np.linspace(1, 0, 10)
#etas = np.linspace(1, 0, 100)
M = 1000
zeta = 1
N = 0



## This generates the data for value of determinant plus error bars. 
## for 2x2 matrices
# as a function of kappa
# fixed alpha, t, M, N

data_ks_y=[]
data_ks_err=[]
for ls in combs_22:
    L = len(etas)
    y=np.zeros((L))
    err=np.zeros((L))
    for i in range(L):
        #y[i]=TD_det(20,sub_state,2,ls,etas[i],N)
        y[i]=TD_det(N_fock,sub_state,2,ls,etas[i],N)
        #err[i]=error(20,sub_state,2,ls,etas[i],N,M)
        err[i]=error(N_fock,sub_state,2,ls,etas[i],N,M)
        #print(y,err)
    data_ks_y.append(y)
    data_ks_err.append(err)
    
    

"FIG 4 PHOTON SUBTRACTED TWO-MODE SQUEEZED STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(etas,0.95*np.ones(len(etas)),'k--',label='95%')

ys = data_ks_y
errs = data_ks_err

for i in range(len(combs_22)):
    c1 = ax.plot(etas, confidence_level(ys[i],errs[i]),label=str(combs_22[i]))


plt.xlabel(r'$\eta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(loc=(0.3,0.2))
plt.title(r'M={}, $\zeta$={}, N={}'.format(M,zeta,N))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xlim([etas[0],etas[-1]])
ax.set_ylim([0.5,1.01])

plt.grid(which="both")

plt.savefig('confidence_ks_sub.svg')

plt.show()

#%%

print(datetime.now() - start)
# The mission was succesful and you can now find optimal NPT criteria :) 
speak("Mission complete. N P T optimised.")