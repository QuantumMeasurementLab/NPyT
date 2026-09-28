#%%
"""
Authors: Lydia A. Kanari-Naish, Amaya Calvo-Sánchez, and Arjun Gupta. Imperial College London
Last update: July 2026

EXAMPLE 1: Two-mode squeezed vacuum state
                                                       

Imports functions from NPyT 
Code for optimizing over all suitable NPT criteria for the example of a 
two-mode squeezed vacuum state
The code to generate the data is provided in .npz files for speed.

"""


"Import packages"
from itertools import cycle
from appscript import con
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
# from win32com.client import Dispatch
# speak = Dispatch("SAPI.SpVoice").Speak
# from datetime import datetime
# start = datetime.now()
import string


"Import the functions from the NPyT"
from NPyT import *

#%%

from tqdm import tqdm
from scipy.optimize import root
from scipy.optimize import least_squares

#%%

solidlinewidth = 3
dottedlinewidth = 6

#%%

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

#%%

# Parameters for fig 1 and 2 (dep on measurement sq parameter)
zs = np.linspace(0.01,2,50) 
amp=0.0
eta = 0.8#1.0 #0.8
N=0
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

#%%

M=500


#%%

# # COMMENTED OUT FOR SPEED

# ## This generates the data for value of determinant and relevant percentiles.
# # as a function of squeezing parameter z for the TMSV state.
# # fixed M, t, kappa, N

# #2x2
# MC_SAMPS = 1000000
# sh = len(zs)
# data_y_12=np.zeros(sh)
# statistical_power_12 = np.zeros(sh)
# statistical_power_12_int = np.zeros(sh)
# percentiles_12 = np.zeros((sh, 4)) 
# percentiles_12_int = np.zeros((sh, 4))


# # Flag: whether kurtosis check ever exceeded threshold for comb [2,4]
# kurt_exceeded_12 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing comb [2,4]'):
#     state = tms_states[i]
#     M_opt = optimize_Mij_gradient_descent(N_fock, state, 1, combs_12[0], eta, N, M, tol=1e-30)[1]
#     M_opt_int = integerisation(M, M_opt[0] + 1j * M_opt[1])
#     # compute detector values
#     data_y_12[i] = TD_det(N_fock, state, 2, combs_12[0], eta, N)
#     _, statistical_power_12[i], percentiles_12[i] = statistical_power_mc(N_fock, state, 2, combs_12[0], eta, N, M, M_opt, mc_samps=MC_SAMPS)
#     _, statistical_power_12_int[i], percentiles_12_int[i], kurt_flag = statistical_power_mc(N_fock, state, 2, combs_12[0], eta, N, M, M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_12[i] = bool(kurt_flag)


# #3x3
# data_y_13=np.zeros(sh)
# statistical_power_13 = np.zeros(sh)
# statistical_power_13_int = np.zeros(sh)
# percentiles_13 = np.zeros((sh, 4))
# percentiles_13_int = np.zeros((sh, 4))
# # Flag: whether kurtosis check ever exceeded threshold for each comb
# kurt_exceeded_13 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing comb {combs_13[1]}'):
#     state = tms_states[i]
#     M_opt = optimize_Mij_gradient_descent(N_fock, state, 1, combs_13[1], eta, N, M, tol=1e-30)[1]
#     M_opt_int = integerisation(M, M_opt[0] + 1j * M_opt[1])
#     data_y_13[i] = TD_det(N_fock, state, 2, combs_13[1], eta, N)
#     _, statistical_power_13[i], percentiles_13[i] = statistical_power_mc(N_fock, state, 2, combs_13[1], eta, N, M, M_opt, mc_samps=MC_SAMPS)
#     _, statistical_power_13_int[i], percentiles_13_int[i], kurt_flag = statistical_power_mc(N_fock, state, 1, combs_13[1], eta, N, M, M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_13[i] = bool(kurt_flag)


# #4x4
# data_y_14=np.zeros(sh)
# statistical_power_14 = np.zeros(sh)
# statistical_power_14_int = np.zeros(sh)
# percentiles_14 = np.zeros((sh, 4))
# percentiles_14_int = np.zeros((sh, 4))
# # Flag: whether kurtosis check ever exceeded threshold for each comb
# kurt_exceeded_14 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing comb {combs_14[1]}'):
#     state = tms_states[i]
#     M_opt = optimize_Mij_gradient_descent(N_fock, state, 1, combs_14[1], eta, N, M, tol=1e-30)[1]
#     M_opt_int = integerisation(M, M_opt[0] + 1j * M_opt[1])
#     data_y_14[i] = TD_det(N_fock, state, 2, combs_14[1], eta, N)
#     _, statistical_power_14[i], percentiles_14[i] = statistical_power_mc(N_fock, state, 2, combs_14[1], eta, N, M, M_opt, mc_samps=MC_SAMPS)
#     _, statistical_power_14_int[i], percentiles_14_int[i], kurt_flag = statistical_power_mc(N_fock, state, 1, combs_14[1], eta, N, M, M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_14[i] = bool(kurt_flag)


#%%
data_y_12 = np.load('data_files/TMSV_data.npz')['data_y_12']
data_y_13 = np.load('data_files/TMSV_data.npz')['data_y_13']
data_y_14 = np.load('data_files/TMSV_data.npz')['data_y_14']
percentiles_12_int = np.load('data_files/TMSV_data.npz')['percentiles_12_int']
percentiles_13_int = np.load('data_files/TMSV_data.npz')['percentiles_13_int']
percentiles_14_int = np.load('data_files/TMSV_data.npz')['percentiles_14_int']

#%%

"FIG 1 TWO-MODE SQUEEZED STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.hlines(0, zs[0], zs[-1], colors='black', zorder=0.5)

colors = ['#0072B2', '#FF9956', '#F0E442']

ax.fill_between(zs, percentiles_14_int[:,3], percentiles_14_int[:,1],alpha=0.3, facecolor='#FFDD00', lw=0.5, zorder=0.5)
ax.fill_between(zs, percentiles_13_int[:,3], percentiles_13_int[:,1],alpha=0.3, facecolor='#DD8452', lw=0.5, zorder=0.5)

ax.fill_between(zs, percentiles_12_int[:,3], percentiles_12_int[:,1],alpha=0.3, facecolor='#4C72B0', lw=0.5, zorder=0.5)

ax.plot(zs, data_y_12, label = r'D$_{I}$', color=colors[0], zorder=3, lw=solidlinewidth)
ax.plot(zs, data_y_13, label = r'D$_{II}$', color=colors[1], zorder=3, lw=solidlinewidth)
ax.plot(zs, data_y_14, label = r'D$_{III}$', color=colors[2], zorder=3, lw=solidlinewidth)


#Analytic expression for DI in the absence of loss  
#ax.plot(zs,-np.sinh(0.5*zs)**2, "-.", label='DI')
#Analytic expression for DII/III in the absence of loss  
#ax.plot(zs,-np.sinh(0.5*zs)**2*np.cosh(0.5*zs)**2, "-.", label='DII/III')
#Analytic expression for root(product of EPR variances)-1 the absence of loss  
#ax.plot(zs,-1+np.exp(-zs), "-.", label='EPR')

plt.xlabel(r'$\zeta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Determinant',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(loc=(0.3,0.25))
plt.title(r'$M_{{tot}}$ = {} $\eta$ = {} $\bar N$ = {}'.format(M,eta,N))


ax.set_xlim([zs[2],zs[-1]])


# symmetric log scale for both positive and negative y values
plt.yscale("symlog")

ax.set_box_aspect(0.65)

plt.grid(zorder=0)
ax.yaxis.set_major_locator(MultipleLocator(1))

plt.savefig('TMSV_fig1.pdf')
plt.savefig('TMSV_fig1.svg')

plt.show()

#%%

statistical_power_12 = np.load('data_files/TMSV_statistical_powers.npz')['statistical_power_12']
statistical_power_13 = np.load('data_files/TMSV_statistical_powers.npz')['statistical_power_13']
statistical_power_14 = np.load('data_files/TMSV_statistical_powers.npz')['statistical_power_14']
statistical_power_12_int = np.load('data_files/TMSV_statistical_powers.npz')['statistical_power_12_int']
statistical_power_13_int = np.load('data_files/TMSV_statistical_powers.npz')['statistical_power_13_int']
statistical_power_14_int = np.load('data_files/TMSV_statistical_powers.npz')['statistical_power_14_int']
kurt_exceeded_12 = np.load('data_files/TMSV_statistical_powers.npz')['kurt_exceeded_12']
kurt_exceeded_13 = np.load('data_files/TMSV_statistical_powers.npz')['kurt_exceeded_13']
kurt_exceeded_14 = np.load('data_files/TMSV_statistical_powers.npz')['kurt_exceeded_14']

#%%
"FIG 2 TWO-MODE SQUEEZED STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.plot(zs,0.95*np.ones(len(zs)),'k--',label=r'95%', lw=solidlinewidth)

plot_clt(kurt_exceeded_12, statistical_power_12_int, zs, 'solid', colors[0])
plot_clt(kurt_exceeded_13, statistical_power_13_int, zs, 'solid', colors[1])
plot_clt(kurt_exceeded_14, statistical_power_14_int, zs, 'solid', colors[2])



plt.xlabel(r'$\zeta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.xticks(fontsize=12)
#generate labels
labelhandles = [plt.Line2D([0], [0], color=colors[i], lw=solidlinewidth) for i in range(3)]
labelhandles.append(plt.Line2D([0], [0], color='k', lw=solidlinewidth, ls='--'))
fig.legend(handles=labelhandles, labels=['D$_{I}$', 'D$_{II}$', 'D$_{III}$', '95%'], loc=(0.4,0.25))
plt.title(r'$M_{{tot}}$ = {} $\eta$ = {} $\bar N$ = {}'.format(M,eta,N))
ax.set_xlim([zs[0],zs[-1]])
ax.set_ylim([0.5,1.02])
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_box_aspect(0.65)

plt.grid()

plt.savefig('TMSV_fig2.svg')
plt.savefig('TMSV_fig2.pdf')

plt.show()


#%%

# Parameters for fig 3 (dependence on measurement number)
Ms = np.round([10**i for i in np.linspace(1,3,50)]) 
# must always ensure integer inputs for M

amp=0.0
eta = 0.8
zeta = 1
N = 0

tms_state=TMSV_vec(N_fock,zeta,amp)

#%%

# #2x2
# MC_SAMPS = 1000000
# sh = len(Ms)
# statistical_power_ms_12 = np.zeros(sh)
# statistical_power_ms_12_int = np.zeros(sh)


# # Flag: whether kurtosis check ever exceeded threshold for comb [2,4]
# kurt_exceeded_ms_12 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing comb [2,4]'):
#     M_opt = optimize_Mij_gradient_descent(N_fock, tms_state, 1, combs_12[0], eta, N, Ms[i], tol=1e-30)[1]
#     M_opt_int = integerisation(Ms[i], M_opt[0] + 1j * M_opt[1])
#     _, statistical_power_ms_12[i], _ = statistical_power_mc(N_fock, tms_state, 2, combs_12[0], eta, N, Ms[i], M_opt, mc_samps=MC_SAMPS)
#     _, statistical_power_ms_12_int[i], _, kurt_flag = statistical_power_mc(N_fock, tms_state, 2, combs_12[0], eta, N, Ms[i], M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_ms_12[i] = bool(kurt_flag)


# #3x3
# statistical_power_ms_13 = np.zeros(sh)
# statistical_power_ms_13_int = np.zeros(sh)
# # Flag: whether kurtosis check ever exceeded threshold for each comb
# kurt_exceeded_ms_13 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing comb {combs_13[1]}'):
#     M_opt = optimize_Mij_gradient_descent(N_fock, tms_state, 1, combs_13[1], eta, N, Ms[i], tol=1e-30)[1]
#     M_opt_int = integerisation(Ms[i], M_opt[0] + 1j * M_opt[1])
#     _, statistical_power_ms_13[i], _ = statistical_power_mc(N_fock, tms_state, 2, combs_13[1], eta, N, Ms[i], M_opt, mc_samps=MC_SAMPS)
#     _, statistical_power_ms_13_int[i], _, kurt_flag = statistical_power_mc(N_fock, tms_state, 1, combs_13[1], eta, N, Ms[i], M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_ms_13[i] = bool(kurt_flag)


# #4x4
# statistical_power_ms_14 = np.zeros(sh)
# statistical_power_ms_14_int = np.zeros(sh)

# # Flag: whether kurtosis check ever exceeded threshold for each comb
# kurt_exceeded_ms_14 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing comb {combs_14[1]}'):
#     M_opt = optimize_Mij_gradient_descent(N_fock, tms_state, 1, combs_14[1], eta, N, Ms[i], tol=1e-30)[1]
#     M_opt_int = integerisation(Ms[i], M_opt[0] + 1j * M_opt[1])
#     _, statistical_power_ms_14[i], _ = statistical_power_mc(N_fock, tms_state, 2, combs_14[1], eta, N, Ms[i], M_opt, mc_samps=MC_SAMPS)
#     _, statistical_power_ms_14_int[i], _, kurt_flag = statistical_power_mc(N_fock, tms_state, 1, combs_14[1], eta, N, Ms[i], M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_ms_14[i] = bool(kurt_flag)


#%%
statistical_power_ms_12 = np.load('data_files/TMSV_statistical_powers_ms.npz')['statistical_power_ms_12']
statistical_power_ms_13 = np.load('data_files/TMSV_statistical_powers_ms.npz')['statistical_power_ms_13']
statistical_power_ms_14 = np.load('data_files/TMSV_statistical_powers_ms.npz')['statistical_power_ms_14']
statistical_power_ms_12_int = np.load('data_files/TMSV_statistical_powers_ms.npz')['statistical_power_ms_12_int']
statistical_power_ms_13_int = np.load('data_files/TMSV_statistical_powers_ms.npz')['statistical_power_ms_13_int']
statistical_power_ms_14_int = np.load('data_files/TMSV_statistical_powers_ms.npz')['statistical_power_ms_14_int']
kurt_exceeded_ms_12 = np.load('data_files/TMSV_statistical_powers_ms.npz')['kurt_exceeded_ms_12']
kurt_exceeded_ms_13 = np.load('data_files/TMSV_statistical_powers_ms.npz')['kurt_exceeded_ms_13']
kurt_exceeded_ms_14 = np.load('data_files/TMSV_statistical_powers_ms.npz')['kurt_exceeded_ms_14']
#%%

"FIG 3 TWO-MODE SQUEEZED STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.plot(Ms,0.95*np.ones(len(Ms)),'k--',label='95%', lw=solidlinewidth)

plot_clt(kurt_exceeded_ms_12, statistical_power_ms_12_int, Ms, 'solid', colors[0])
plot_clt(kurt_exceeded_ms_13, statistical_power_ms_13_int, Ms, 'solid', colors[1])
plot_clt(kurt_exceeded_ms_14, statistical_power_ms_14_int, Ms, 'solid', colors[2])


plt.xlabel(r'$M_{tot}$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.xticks(fontsize=12)
#generate labels
labelhandles = [plt.Line2D([0], [0], color=colors[i], lw=solidlinewidth) for i in range(3)]
labelhandles.append(plt.Line2D([0], [0], color='k', lw=solidlinewidth, ls='--'))
fig.legend(handles=labelhandles, labels=['D$_{I}$', 'D$_{II}$', 'D$_{III}$', '95%'], loc=(0.7,0.25))
plt.title(r'$\zeta$ = {} $\eta$ = {} $\bar N$ = {}'.format(zeta,eta,N))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xscale('log')

plt.grid(which='both')
ax.set_xlim([Ms[0],Ms[-1]])
plt.xlim(40,1000)
plt.ylim(0.5,1.02)

ax.set_box_aspect(0.65)

plt.savefig('TMSV_fig3.svg')
plt.savefig('TMSV_fig3.pdf')


plt.show()

#%%

# Parameters for fig 4 (dependence on eta)
etas = np.linspace(1, 0.001, 100)
amp=0.0
zeta = 1.0
N = 0

# tms_state=TMSV_vec(N_fock,zeta,amp)

#%%

# Code to generate figure 4 data 

# # Monte Carlo samples 
# MC_SAMPS = 1000000

# sh = len(etas)
# statistical_power_ks_12 = np.zeros(sh)
# statistical_power_ks_12_int = np.zeros(sh)
# statistical_power_ks_13 = np.zeros(sh)
# statistical_power_ks_13_int = np.zeros(sh)
# statistical_power_ks_14 = np.zeros(sh)
# statistical_power_ks_14_int = np.zeros(sh)

# kurt_exceeded_ks_12 = np.zeros(sh, dtype=bool)
# kurt_exceeded_ks_13 = np.zeros(sh, dtype=bool)
# kurt_exceeded_ks_14 = np.zeros(sh, dtype=bool)

# # Process each combination and eta, save intermediate results
# for i in tqdm(range(len(etas)), desc=f'Processing comb {combs_12[0]} (2x2)'):
#     eta_val = etas[i]
#     M_opt = optimize_Mij_gradient_descent(N_fock, tms_state, 1, combs_12[0], eta_val, N, M, tol=1e-30)[1]
#     M_opt_int = integerisation(M, M_opt[0] + 1j*M_opt[1])
#     _, pow,_ = statistical_power_mc(N_fock, tms_state, 1, combs_12[0], eta_val, N, M, M_opt, mc_samps=MC_SAMPS)
#     _, pow_int, _, kurt_flag = statistical_power_mc(N_fock, tms_state, 1, combs_12[0], eta_val, N, M, M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     statistical_power_ks_12[i] = pow
#     statistical_power_ks_12_int[i] = pow_int
#     kurt_exceeded_ks_12[i] = bool(kurt_flag)
#     # save progress so the run can be resumed if interrupted
#     # np.save(f'statistical_power_ks_12_comb{j}.npy', statistical_power_ks_12[j])

# for i in tqdm(range(len(etas)), desc=f'Processing comb {combs_13[1]} (3x3)'):
#     eta_val = etas[i]
#     M_opt = optimize_Mij_gradient_descent(N_fock, tms_state, 1, combs_13[1], eta_val, N, M, tol=1e-30)[1]
#     M_opt_int = integerisation(M, M_opt[0] + 1j*M_opt[1])
#     _, pow,_ = statistical_power_mc(N_fock, tms_state, 1, combs_13[1], eta_val, N, M, M_opt, mc_samps=MC_SAMPS)
#     _, pow_int, _, kurt_flag = statistical_power_mc(N_fock, tms_state, 1, combs_13[1], eta_val, N, M, M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     statistical_power_ks_13[i] = pow
#     statistical_power_ks_13_int[i] = pow_int
#     kurt_exceeded_ks_13[i] = bool(kurt_flag)
# # np.save(f'statistical_power_ks_13_comb{j}.npy', statistical_power_ks_13[j])

# for i in tqdm(range(len(etas)), desc=f'Processing comb {combs_14[1]} (4x4)'):
#     eta_val = etas[i]
#     M_opt = optimize_Mij_gradient_descent(N_fock, tms_state, 1, combs_14[1], eta_val, N, M, tol=1e-30)[1]
#     M_opt_int = integerisation(M, M_opt[0] + 1j*M_opt[1])
#     _, pow,_ = statistical_power_mc(N_fock, tms_state, 1, combs_14[1], eta_val, N, M, M_opt, mc_samps=MC_SAMPS)
#     _, pow_int, _, kurt_flag = statistical_power_mc(N_fock, tms_state, 1, combs_14[1], eta_val, N, M, M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     statistical_power_ks_14[i] = pow
#     statistical_power_ks_14_int[i] = pow_int
#     kurt_exceeded_ks_14[i] = bool(kurt_flag)
# # np.save(f'statistical_power_ks_14_comb{j}.npy', statistical_power_ks_14[j])


#%%

statistical_power_ks_12 = np.load('data_files/TMSV_statistical_powers_ks.npz')['statistical_power_ks_12']
statistical_power_ks_13 = np.load('data_files/TMSV_statistical_powers_ks.npz')['statistical_power_ks_13']
statistical_power_ks_14 = np.load('data_files/TMSV_statistical_powers_ks.npz')['statistical_power_ks_14']
statistical_power_ks_12_int = np.load('data_files/TMSV_statistical_powers_ks.npz')['statistical_power_ks_12_int']
statistical_power_ks_13_int = np.load('data_files/TMSV_statistical_powers_ks.npz')['statistical_power_ks_13_int']
statistical_power_ks_14_int = np.load('data_files/TMSV_statistical_powers_ks.npz')['statistical_power_ks_14_int']
kurt_exceeded_ks_12 = np.load('data_files/TMSV_statistical_powers_ks.npz')['kurt_exceeded_ks_12']
kurt_exceeded_ks_13 = np.load('data_files/TMSV_statistical_powers_ks.npz')['kurt_exceeded_ks_13']
kurt_exceeded_ks_14 = np.load('data_files/TMSV_statistical_powers_ks.npz')['kurt_exceeded_ks_14']
#%%

"FIG 4 TWO-MODE SQUEEZED STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.plot(etas,0.95*np.ones(len(etas)),'k--',label=r'95%', lw=solidlinewidth)


plot_clt(kurt_exceeded_ks_12, statistical_power_ks_12_int, etas, 'solid', colors[0])
plot_clt(kurt_exceeded_ks_13, statistical_power_ks_13_int, etas, 'solid', colors[1])
plot_clt(kurt_exceeded_ks_14, statistical_power_ks_14_int, etas, 'solid', colors[2])

plt.xlabel(r'$\eta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.yticks(fontsize=12)
# generate labels
labelhandles = [plt.Line2D([0], [0], color=colors[i], lw=solidlinewidth) for i in range(3)]
labelhandles.append(plt.Line2D([0], [0], color='k', lw=solidlinewidth, ls='--'))
fig.legend(handles=labelhandles, labels=['D$_{I}$', 'D$_{II}$', 'D$_{III}$', '95%'], loc=(0.3,0.25))

plt.title(r'$M_{{tot}}$ = {} $\zeta$ = {} $\bar N$ = {}'.format(M,int(zeta),N))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xlim([etas[0],etas[-1]])
ax.set_ylim([0.50,1.02])

ax.set_box_aspect(0.65)

plt.grid(which="both")

plt.savefig('TMSV_fig4.svg')
plt.savefig('TMSV_fig4.pdf')

plt.show()

#%%

# print(datetime.now() - start)
# # The mission was succesful and you can now find optimal NPT criteria :) 
# speak("Mission complete. N P T optimised.")

