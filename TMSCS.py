"""
Author: Lydia A. Kanari-Naish, Imperial College London
Last update: September 2023

EXAMPLE 3: Two-mode Schrodinger cat state
"""

#     /\_____/\            /\_____/\
#    /  O   x  \          /  X   o  \
#   ( ==  ^  == )   +    ( ==  ^  == )
#    )         (          )         (
#   (           )        (           )
#  ( (  )   (  ) )      ( (  )   (  ) )
# (__(__)___(__)__)    (__(__)___(__)__) 
   
"""  
Imports functions from NPyT
Code for optimizing over all suitable NPT criteria for example of a 
two-mode SCS

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

"Define the two-mode SCS"


def state_generator(fock_dims, alpha, beta, gamma, delta, phi):
    #this gives a general 2 mode cat state of the form |alpha>|beta>+exp(i\phi)|gamma>|delta>
    state=(tensor(coherent(fock_dims,alpha),coherent(fock_dims,beta))+np.exp(1j*phi)*tensor(coherent(fock_dims,gamma),coherent(fock_dims,delta))).unit()
    return state

#this creates a vectorised version
vec_states=np.frompyfunc(state_generator,6,1)



# create range of alphas to consider
alphas = np.linspace(0.0001,2,10) # Smaller vector for shorter total run time
#alphas = np.linspace(0.0001,2,100) # This vector of alphas was used for the results in the paper.
#cat_alphas = vec_states(20,alphas,0,0,alphas,np.pi)
cat_alphas = vec_states(N_fock,alphas,0,0,alphas,np.pi)


"Preliminary search of successful determinants"
# For order 2, matrix dimensions 2 we can find suitable determinants with a preliminary search
#vals_22,combs_22 = my_state(2,2,20,cat_alphas)
vals_22,combs_22 = my_state(2,2,N_fock,cat_alphas)


# Here, we notice that some of these are repeats because of the symmetry of the cat state.
# i.e. we could change the labels of the operators from subsystem A to B and some determinants are the same.
# So here we identify the set of unique determinants from this preliminary search
combs_22=np.array([[ 0, 12],[1,12] ,[ 5, 12], [ 8, 12], ])



# Total number of measurments
M=1e3
# Phonon occupation number of bath
N=0.01
# Optical losses
eta = 0.95



TD_det_vec = np.frompyfunc(TD_det,7,1)



## This generates the data for value of determinant plus error bars. 
## for 2x2 matrices
# as a function of alpha
# fixed M, eta, N
data_22_alphas_y=[]
data_22_alphas_err=[]
for ls in combs_22:
    L = len(cat_alphas)
    y=np.zeros((L))
    err=np.zeros((L))
    for i in range(L):
        y[i]=TD_det(N_fock,cat_alphas[i],2,ls,eta,N)
        err[i]=error(N_fock,cat_alphas[i],2,ls,eta,N,M)
        #print(y,err)
    data_22_alphas_y.append(y)
    data_22_alphas_err.append(err)
data_22_alphas_y=np.array([data_22_alphas_y])
data_22_alphas_err=np.array([data_22_alphas_err])




sh = len(cat_alphas)
data_alphas_S3_y=np.zeros(sh)
data_alphas_S3_err=np.zeros(sh)
for i in range(sh):
    data_alphas_S3_y[i]=TD_det(N_fock,cat_alphas[i],2,[0,4,12],eta,N)
    data_alphas_S3_err[i]=error(N_fock,cat_alphas[i],2,[0,4,12],eta,N,M)
    
data_alphas_S3_y = np.array([data_alphas_S3_y])
data_alphas_S3_err = np.array([data_alphas_S3_err])


"FIG 1 TWO MODE SCHRODINGER CAT STATE"

fig = plt.figure(figsize=(8,15))
fig, ax = plt.subplots()

ax.axhline(0,color='black',alpha=0.8)

for i in range(len(combs_22)):
    y=data_22_alphas_y[0][i]
    err=data_22_alphas_err[0][i]
    c1 = ax.plot(alphas, y, '-',label = str(combs_22[i]))
    c2 = ax.fill_between(alphas, y-err, y+err,alpha=0.2)
    
c3 = ax.plot(alphas, data_alphas_S3_y[0],label='S3')
c4 = ax.fill_between(alphas,data_alphas_S3_y[0]-data_alphas_S3_err[0],data_alphas_S3_y[0]+data_alphas_S3_err[0],alpha=0.2)

plt.xlabel(r'$\alpha$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Determinant',fontsize=14)
plt.yticks(fontsize=12)
plt.ylim([-0.25, 0.25])

fig.legend(loc=(0.2,0.6))
plt.title(r'Determinant for M={}, $\eta$={}, N={}'.format(M,eta,N))

ax.set_xlim([alphas[0],alphas[-1]])

plt.grid()

plt.savefig('det_error_bars.svg')

plt.show()


"FIG 2 TWO MODE SCHRODINGER CAT STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(alphas,0.95*np.ones(len(alphas)),'k--',label='95%')

ys = data_22_alphas_y
errs = data_22_alphas_err

for i in range(len(combs_22)):
    c1 = ax.plot(alphas, confidence_level(ys[0][i],errs[0][i]),label=str(combs_22[i]))
    

c2 = ax.plot(alphas, confidence_level(data_alphas_S3_y[0],data_alphas_S3_err[0]),label='S3' )

plt.ylim([0.35, 1.02])

plt.xlabel(r'$\alpha$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.yticks(fontsize=12)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol='%', is_latex=False))
fig.legend(loc=(0.4,0.2))
plt.title(r'M={}, $\eta$={}, N={}'.format(M,eta,N))
ax.set_xlim([alphas[0],alphas[-1]])



plt.grid()

plt.savefig('confidence.svg')

plt.show()



# Data for how confidence depends on measurement number

# Plot as a function of M, number of measurement
Ms = [10**i for i in np.linspace(2,4,10)]
#Ms = [10**i for i in np.linspace(2,4,100)]
# Keeping optical coupling fixed
eta = 0.95
# alpha is cat state parameter
alpha = 1
# phonon occupation number of bath
N = 0.01




#cat_state=state_generator(20,alpha,0,0,alpha,np.pi)
cat_state=state_generator(N_fock,alpha,0,0,alpha,np.pi)



## This generates the data for value of determinant plus error bars. 
## for 2x2 matrices
# as a function of M
# fixed alpha, t, kappa, N
data_22_ms_y=[]
data_22_ms_err=[]
for ls in combs_22:
    L = len(Ms)
    y=np.zeros((L))
    err=np.zeros((L))
    for i in range(L):
        #y[i]=TD_det(20,cat_state,2,ls,eta,N)
        y[i]=TD_det(N_fock,cat_state,2,ls,eta,N)
        # err[i]=error(20,cat_state,2,ls,eta,N,Ms[i])
        err[i]=error(N_fock,cat_state,2,ls,eta,N,Ms[i])
        #print(y,err)
    data_22_ms_y.append(y)
    data_22_ms_err.append(err)




sh = len(Ms)
data_ms_S3_y=np.zeros(sh)
data_ms_S3_err=np.zeros(sh)
for i in range(sh):
    data_ms_S3_y[i]=TD_det(N_fock,cat_state,2,[0,4,12],eta,N)
    data_ms_S3_err[i]=error(N_fock,cat_state,2,[0,4,12],eta,N,Ms[i])


"FIG 3 TWO MODE SCHRODINGER CAT STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(Ms,0.95*np.ones(len(Ms)),'k--',label='95%')

ys = data_22_ms_y
errs = data_22_ms_err

for i in range(len(combs_22)):
    c1 = ax.plot(Ms, confidence_level(ys[i],errs[i]),label=str(combs_22[i]))
c2 = ax.plot(Ms, confidence_level(data_ms_S3_y,data_ms_S3_err),label='S3' )



plt.xlabel(r'$M_{tot}$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(loc=(0.7,0.2))
plt.title(r'$\alpha$={}, $\eta$={}, N={}'.format(alpha,eta,N))

ax.set_xscale('log')
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))


plt.grid(which='both')
ax.set_xlim([Ms[0],Ms[-1]])
#ax.set_ylim(0.5,1.02)
ax.set_ylim(0.62,1.02)

plt.show()

plt.savefig('confidence_ms.svg')



# Data for how confidence depends on measurement number
etas = np.linspace(1,0.00,10)
#etas = np.linspace(1,0.00,100)
M = 1000
alpha = 1
N = 0.001




cat_state=state_generator(N_fock,alpha,0,0,alpha,np.pi)



## This generates the data for value of determinant plus error bars. 
## for 2x2 matrices
# as a function of kappa
# fixed alpha, t, M, N


data_22_ks_y=[]
data_22_ks_err=[]
for ls in combs_22:
    L = len(etas)
    y=np.zeros((L))
    err=np.zeros((L))
    for i in range(L):
        #y[i]=TD_det(20,cat_state,2,ls,etas[i],N)
        y[i]=TD_det(N_fock,cat_state,2,ls,etas[i],N)
        #err[i]=error(20,cat_state,2,ls,etas[i],N,M)
        err[i]=error(N_fock,cat_state,2,ls,etas[i],N,M)
        #print(y,err)
    data_22_ks_y.append(y)
    data_22_ks_err.append(err)



sh = len(etas)
data_ks_S3_y=np.zeros(sh)
data_ks_S3_err=np.zeros(sh)
for i in range(sh):
    data_ks_S3_y[i]=TD_det(N_fock,cat_state,2,[0,4,12],etas[i],N)
    data_ks_S3_err[i]=error(N_fock,cat_state,2,[0,4,12],etas[i],N,M)


"FIG 4 TWO MODE SCHRODINGER CAT STATE"


fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(etas,0.95*np.ones(len(etas)),'k--',label='95%')

ys = data_22_ks_y
errs = data_22_ks_err

for i in range(len(combs_22)):
    c1 = ax.plot(etas, confidence_level(ys[i],errs[i]),label=str(combs_22[i]))

c2 = ax.plot(etas, confidence_level(data_ks_S3_y,data_ks_S3_err),label='S3' )


plt.xlabel(r'$\eta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(loc=(0.3,0.2))
plt.title(r'M={}, $\alpha$={}, N={}'.format(M,alpha,N))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xlim([etas[0],etas[-1]])
ax.set_ylim([0.35,1.05])

plt.grid(which="both")

plt.savefig('confidence_ks.svg')

plt.show()


#%%

print(datetime.now() - start)
# The mission was succesful and you can now find optimal NPT criteria :) 
speak("Mission complete. N P T optimised.")
