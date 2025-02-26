"""
Author: Lydia A. Kanari-Naish, Imperial College London
Last update: September 2023

EXAMPLE 1: Two-mode squeezed vacuum state
                                                       

Imports functions from NPyT 
Code for optimizing over all suitable NPT criteria for example of a 
two-mode squeezed vacuum state

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



"Define the displaced TMSV"
def sq(fock_dims, z, amp):
    a1 = tensor(destroy(fock_dims),qeye(fock_dims))
    a2 = tensor(qeye(fock_dims),destroy(fock_dims))
    # Displace each mode of TMSV with amp 
    disp_2=tensor(displace(fock_dims, amp),displace(fock_dims, amp))

    ground_state = tensor(basis(fock_dims, 0),basis(fock_dims,0))
    TMSV = disp_2*squeezing(a1, a2, z)*ground_state
    return TMSV.unit()


TMSV_vec=np.frompyfunc(sq,3,1)



# Parameters for fig 1 and 2 (dep on measurement sq parameter)
zs = np.linspace(0.0001,2,10) #turn number of samples much higher for actual plotting
amp=0.0
eta = 1.0 #0.8
N=0
M=200
tms_states = TMSV_vec(N_fock,zs,amp)


"Preliminary search of successful determinants"
"Submatrices of order=2 from 2x2 to 5x5 (Simon's criterion') that can see entanglement"
vals_12,combs_12=my_state(1,2,N_fock,np.array([tms_states],dtype=object))
vals_13,combs_13=my_state(1,3,N_fock,np.array([tms_states],dtype=object))
vals_14,combs_14=my_state(1,4,N_fock,np.array([tms_states],dtype=object))
# Simon's criterion
vals_15,combs_15=my_state(1,5,N_fock,np.array([tms_states],dtype=object))

# As the TMSV is symmetric under the swapping of the modes, we find that
# for the 3x3 matrix determinants [1,2,4] is equivalent to [2,3,4] and so
# remove [2, 3, 4] just for convenience
combs_13=np.array([[0,  2,  4],[1, 2, 4]]) 
#for the 4x4 matrix determinants [0,1,2,4] is equivalent to [0,2,3,4] so
# remove [0, 2, 3, 4] just for convenience
combs_14=np.array([[0,  1, 2,  4],[1, 2, 3, 4]]) 



## This generates the data for value of determinant plus error bars. 
# as a function of ks
# fixed M, t, kappa, N

#2x2
data_y_12=[]
data_err_12=[]
for ls in combs_12:
    sh = len(zs)
    y_12=np.zeros((sh))
    err_12=np.zeros((sh))
    for i in range(len(zs)):
        y_12[i]=(TD_det(N_fock,tms_states[i],2,ls,eta,N))
        err_12[i]=(error(N_fock,tms_states[i],2,ls,eta,N,M))
 
    data_y_12.append(y_12)
    data_err_12.append(err_12)
data_y_12=np.array([data_y_12])
data_err_12=np.array([data_err_12])

#3x3
data_y_13=[]
data_err_13=[]
for ls in combs_13:
    sh = len(zs)
    y_13=np.zeros((sh))
    err_13=np.zeros((sh))
    for i in range(len(zs)):
        y_13[i]=(TD_det(N_fock,tms_states[i],2,ls,eta,N))
        err_13[i]=(error(N_fock,tms_states[i],2,ls,eta,N,M))
 
    data_y_13.append(y_13)
    data_err_13.append(err_13)
data_y_13=np.array([data_y_13])
data_err_13=np.array([data_err_13])

#4x4
data_y_14=[]
data_err_14=[]
for ls in combs_14:
    sh = len(zs)
    y_14=np.zeros((sh))
    err_14=np.zeros((sh))
    for i in range(len(zs)):
        y_14[i]=(TD_det(N_fock,tms_states[i],2,ls,eta,N))
        err_14[i]=(error(N_fock,tms_states[i],2,ls,eta,N,M))
 
    data_y_14.append(y_14)
    data_err_14.append(err_14)
data_y_14=np.array([data_y_14])
data_err_14=np.array([data_err_14])

#5x5
data_y_15=[]
data_err_15=[]
for ls in combs_15:
    sh = len(zs)
    y_15=np.zeros((sh))
    err_15=np.zeros((sh))
    for i in range(len(zs)):
        y_15[i]=(TD_det(N_fock,tms_states[i],2,ls,eta,N))
        err_15[i]=(error(N_fock,tms_states[i],2,ls,eta,N,M))
 
    data_y_15.append(y_15)
    data_err_15.append(err_15)
data_y_15=np.array([data_y_15])
data_err_15=np.array([data_err_15])




"FIG 1 TWO-MODE SQUEEZED STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

## Different linestyles to cycle through
#lines = ["-","--",":","-."] 
lines = ["-","--"]
linecycler = cycle(lines)

for i in range(len(combs_12)):
    y=data_y_12[0][i]
    err=data_err_12[0][i]
    c1 = ax.plot(zs, y, next(linecycler), label = str(combs_12[i]),)
    c2 = ax.fill_between(zs, y-err, y+err,alpha=0.2)
    
for i in range(len(combs_13)):
    y=data_y_13[0][i]
    err=data_err_13[0][i]
    c1 = ax.plot(zs, y, next(linecycler), label = str(combs_13[i]),)
    c2 = ax.fill_between(zs, y-err, y+err,alpha=0.2)

for i in range(len(combs_14)):
    y=data_y_14[0][i]
    err=data_err_14[0][i]
    c1 = ax.plot(zs, y, next(linecycler), label = str(combs_14[i]),)
    c2 = ax.fill_between(zs, y-err, y+err,alpha=0.2)
    
for i in range(len(combs_15)):
    y=data_y_15[0][i]
    err=data_err_15[0][i]
    c1 = ax.plot(zs, y, next(linecycler), label = str(combs_15[i]),)
    c2 = ax.fill_between(zs, y-err, y+err,alpha=0.2)

#Analytic expression for DI in the absence of loss  
ax.plot(zs,-np.sinh(0.5*zs)**2, "-.", label='DI')
#Analytic expression for DII/III in the absence of loss  
ax.plot(zs,-np.sinh(0.5*zs)**2*np.cosh(0.5*zs)**2, "-.", label='DII/III')
#Analytic expression for root(product of EPR variances)-1 the absence of loss  
ax.plot(zs,-1+np.exp(-zs), "-.", label='EPR')
plt.xlabel(r'$\zeta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Determinant',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(loc=(0.2,0.2))
plt.title(r'Determinant for M={}, $\eta$={}, N={}, $\alpha$={}'.format(M,eta,N,amp))


ax.set_xlim([zs[0],zs[-1]])


# symmetric log scale for both positive and negative y values
plt.yscale("symlog")

plt.grid()
ax.yaxis.set_major_locator(MultipleLocator(1))

plt.savefig('TMSV_fig1.svg')

plt.show()



"FIG 2 TWO-MODE SQUEEZED STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(zs,0.95*np.ones(len(zs)),'k--',label='95%')

# Different linestyles to cycle through
lines = ["-","--"] 
linecycler = cycle(lines)

ys_12 = data_y_12
errs_12 = data_err_12
for i in range(len(combs_12)):  
    c1 = ax.plot(zs, confidence_level(ys_12[0][i],errs_12[0][i]), next(linecycler), label=str(combs_12[i]))

ys_13 = data_y_13
errs_13 = data_err_13
for i in range(len(combs_13)):   
    c1 = ax.plot(zs, confidence_level(ys_13[0][i],errs_13[0][i]), next(linecycler), label=str(combs_13[i]))
    
ys_14 = data_y_14
errs_14 = data_err_14
for i in range(len(combs_14)):   
    c1 = ax.plot(zs, confidence_level(ys_14[0][i],errs_14[0][i]), next(linecycler), label=str(combs_14[i]))

ys_15 = data_y_15
errs_15 = data_err_15
for i in range(len(combs_15)):   
    c1 = ax.plot(zs, confidence_level(ys_15[0][i],errs_15[0][i]), next(linecycler), label=str(combs_15[i]))


plt.xlabel(r'$\zeta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.xticks(fontsize=12)
fig.legend(loc=(0.4,0.2))
plt.title(r'M={}, $\eta$={}, N={}, $\alpha$={}'.format(M,eta,N,amp))
ax.set_xlim([zs[0],zs[-1]])
ax.set_ylim([0.5,1.01])
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

plt.grid()

plt.savefig('TMSV_fig2.svg')

plt.show()




# Parameters for fig 3 (dep on measurement number)
Ms = [10**i for i in np.linspace(1,3,10)] # increase from 40 in actual plots
amp=0.0
eta = 0.8
zeta = 1
N = 0

tms_state=TMSV_vec(N_fock,zeta,amp)


## This generates the data for value of determinant plus error bars. 
# as a function of M, with fixed alpha, t, kappa, N

# 2x2
data_ms_y_12=[]
data_ms_err_12=[]
for ls in combs_12:
    L = len(Ms)
    y_12=np.zeros((L))
    err_12=np.zeros((L))
    for i in range(L):
        y_12[i]=TD_det(N_fock,tms_state,2,ls,eta,N)
        err_12[i]=error(N_fock,tms_state,2,ls,eta,N,Ms[i])
    data_ms_y_12.append(y_12)
    data_ms_err_12.append(err_12)

# 3x3
data_ms_y_13=[]
data_ms_err_13=[]
for ls in combs_13:
    L = len(Ms)
    y_13=np.zeros((L))
    err_13=np.zeros((L))
    for i in range(L):
        y_13[i]=TD_det(N_fock,tms_state,2,ls,eta,N)
        err_13[i]=error(N_fock,tms_state,2,ls,eta,N,Ms[i])
    data_ms_y_13.append(y_13)
    data_ms_err_13.append(err_13)
    
# 4x4
data_ms_y_14=[]
data_ms_err_14=[]
for ls in combs_14:
    L = len(Ms)
    y_14=np.zeros((L))
    err_14=np.zeros((L))
    for i in range(L):
        y_14[i]=TD_det(N_fock,tms_state,2,ls,eta,N)
        err_14[i]=error(N_fock,tms_state,2,ls,eta,N,Ms[i])
    data_ms_y_14.append(y_14)
    data_ms_err_14.append(err_14)
    
# 5x5
data_ms_y_15=[]
data_ms_err_15=[]
for ls in combs_15:
    L = len(Ms)
    y_15=np.zeros((L))
    err_15=np.zeros((L))
    for i in range(L):
        y_15[i]=TD_det(N_fock,tms_state,2,ls,eta,N)
        err_15[i]=error(N_fock,tms_state,2,ls,eta,N,Ms[i])
    data_ms_y_15.append(y_15)
    data_ms_err_15.append(err_15)


"FIG 3 TWO-MODE SQUEEZED STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(Ms,0.95*np.ones(len(Ms)),'k--',label='95%')

# Different linestyles to cycle through
lines = ["-","--"]  
linecycler = cycle(lines)

ys_12 = data_ms_y_12
errs_12 = data_ms_err_12
for i in range(len(combs_12)):
    c1 = ax.plot(Ms,  confidence_level(ys_12[i],errs_12[i]), next(linecycler),label=str(combs_12[i]))

ys_13 = data_ms_y_13
errs_13 = data_ms_err_13
for i in range(len(combs_13)):
    c1 = ax.plot(Ms,  confidence_level(ys_13[i],errs_13[i]), next(linecycler),label=str(combs_13[i]))

ys_14 = data_ms_y_14
errs_14 = data_ms_err_14
for i in range(len(combs_14)):
    c1 = ax.plot(Ms,  confidence_level(ys_14[i],errs_14[i]), next(linecycler),label=str(combs_14[i]))
    
ys_15 = data_ms_y_15
errs_15 = data_ms_err_15
for i in range(len(combs_15)):
    c1 = ax.plot(Ms,  confidence_level(ys_15[i],errs_15[i]), next(linecycler),label=str(combs_15[i]))


plt.xlabel(r'$M_{tot}$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.xticks(fontsize=12)
fig.legend(loc=(0.7,0.2))
plt.title(r'$\zeta$={}, $\eta$={}, N={}, $\alpha$={}'.format(zeta,eta,N,amp))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xscale('log')

plt.grid(which='both')
ax.set_xlim([Ms[0],Ms[-1]])
plt.ylim(0.65,1.02)

plt.show()

plt.savefig('TMSV_fig3.svg')


# Parameters for fig 4 (dep on eta)
etas = np.linspace(1, 0, 10)
amp=0.0
M = 200
zeta = 1.0
N = 0

########################################################################################################
# To vary the values of zeta and amp between Fig 3 and 4
########################################################################################################
tms_state=TMSV_vec(N_fock,zeta,amp)


## This generates the data for value of determinant plus error bars. 
## as a function of kappa with fixed alpha, t, M, N

#2x2
data_ks_y_12=[]
data_ks_err_12=[]
for ls in combs_12:
    L = len(etas)
    y_12=np.zeros((L))
    err_12=np.zeros((L))
    for i in range(L):
        y_12[i]=TD_det(N_fock,tms_state,2,ls,etas[i],N)
        err_12[i]=error(N_fock,tms_state,2,ls,etas[i],N,M)
    data_ks_y_12.append(y_12)
    data_ks_err_12.append(err_12)
    
#3x3
data_ks_y_13=[]
data_ks_err_13=[]
for ls in combs_13:
    L = len(etas)
    y_13=np.zeros((L))
    err_13=np.zeros((L))
    for i in range(L):
        y_13[i]=TD_det(N_fock,tms_state,2,ls,etas[i],N)
        err_13[i]=error(N_fock,tms_state,2,ls,etas[i],N,M)
    data_ks_y_13.append(y_13)
    data_ks_err_13.append(err_13)
    
#4x4
data_ks_y_14=[]
data_ks_err_14=[]
for ls in combs_14:
    L = len(etas)
    y_14=np.zeros((L))
    err_14=np.zeros((L))
    for i in range(L):
        y_14[i]=TD_det(N_fock,tms_state,2,ls,etas[i],N)
        err_14[i]=error(N_fock,tms_state,2,ls,etas[i],N,M)
    data_ks_y_14.append(y_14)
    data_ks_err_14.append(err_14)
    
#5x5
data_ks_y_15=[]
data_ks_err_15=[]
for ls in combs_15:
    L = len(etas)
    y_15=np.zeros((L))
    err_15=np.zeros((L))
    for i in range(L):
        y_15[i]=TD_det(N_fock,tms_state,2,ls,etas[i],N)
        err_15[i]=error(N_fock,tms_state,2,ls,etas[i],N,M)
    data_ks_y_15.append(y_15)
    data_ks_err_15.append(err_15)


"FIG 4 TWO-MODE SQUEEZED STATE"

fig = plt.figure(figsize=(8,8))
fig, ax = plt.subplots()

plt.plot(etas,0.95*np.ones(len(etas)),'k--',label='95%')

# Different linestyles to cycle through
lines = ["-","--"] 
linecycler = cycle(lines)

ys_12 = data_ks_y_12
errs_12 = data_ks_err_12
for i in range(len(combs_12)): 
    c1 = ax.plot(etas,  confidence_level(ys_12[i],errs_12[i]), next(linecycler), label=str(combs_12[i]))
    
ys_13 = data_ks_y_13
errs_13 = data_ks_err_13
for i in range(len(combs_13)): 
    c1 = ax.plot(etas,  confidence_level(ys_13[i],errs_13[i]), next(linecycler), label=str(combs_13[i]))

ys_14 = data_ks_y_14
errs_14 = data_ks_err_14
for i in range(len(combs_14)):  
    c1 = ax.plot(etas,  confidence_level(ys_14[i],errs_14[i]), next(linecycler), label=str(combs_14[i]))

ys_15 = data_ks_y_15
errs_15 = data_ks_err_15
for i in range(len(combs_15)):
    c1 = ax.plot(etas,  confidence_level(ys_15[i],errs_15[i]), next(linecycler), label=str(combs_15[i]))


plt.xlabel(r'$\eta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Confidence',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(loc=(0.3,0.2))
plt.title(r'M={}, $\zeta$={}, N={}, $\alpha$={}'.format(M,zeta,N,amp))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xlim([etas[0],etas[-1]])
ax.set_ylim([0.40,1.01])


plt.grid(which="both")

plt.savefig('TMSV_fig4.svg')

plt.show()

#%%

print(datetime.now() - start)
# The mission was succesful and you can now find optimal NPT criteria :) 
speak("Mission complete. N P T optimised.")