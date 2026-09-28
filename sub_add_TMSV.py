#%%
"""
Authors: Lydia A. Kanari-Naish, Amaya Calvo-Sánchez, and Arjun Gupta. Imperial College London
Last update: July 2026

EXAMPLE 2: Photon-subtracted/added two-mode squeezed vacuum state


Imports functions from NPyT
Code for optimizing over all suitable NPT criteria for the example of a 
photon-subtracted/added two-mode squeezed vacuum state
The code to generate the data is provided in .npz files for speed.

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
#from win32com.client import Dispatch
#speak = Dispatch("SAPI.SpVoice").Speak
# from datetime import datetime
# start = datetime.now()


"Import the functions from the NPyT"
from NPyT import *

#%%

from tqdm import tqdm

#%%

solidlinewidth = 3
dottedlinewidth = 6

#%%

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


#%%

# Squeezing parameter
zs = np.round(np.linspace(0.0001,2,50), decimals=4)

sub_states = sub_vec(N_fock,zs,1,1)


eta = 0.8
N=0
M=1000

#%%

"Preliminary search of successful determinants"
vals_22,combs_22=my_state(2,2,N_fock,np.array([sub_states],dtype=object))
#Note EIV and EV are not unique as noted in the manuscript so let's removes these from combs_22
combs_22=np.array([[ 2,  4],[ 6, 11],[ 7, 12],[ 7, 14]])

#%%

# # COMMENTED OUT FOR SPEED

# #2x2
# MC_SAMPS = 1000000
# sh = len(zs)
# data_y_22=np.zeros((len(combs_22), sh))
# statistical_power_22 = np.zeros((len(combs_22), sh))
# statistical_power_22_int = np.zeros((len(combs_22), sh))
# percentiles_22 = np.zeros((len(combs_22), sh, 4)) 
# percentiles_22_int = np.zeros((len(combs_22), sh, 4))


# # Flag: whether kurtosis check ever exceeded threshold for comb [2,4]
# kurt_exceeded_22 = np.zeros((len(combs_22), sh), dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing 2x2 combs'):
#     state = sub_states[i]
#     for j in range(len(combs_22)):
#         ls = combs_22[j]
#         M_opt = optimize_Mij_gradient_descent(N_fock, state, 2, combs_22[j], eta, N, M, tol=1e-30)[1]
#         M_opt_int = integerisation(M, M_opt[0] + 1j * M_opt[1], 2, combs_22[j])
#         # compute detector values
#         data_y_22[j, i] = TD_det(N_fock, state, 2, combs_22[j], eta, N)
#         _, statistical_power_22[j, i], percentiles_22[j, i] = statistical_power_mc(N_fock, state, 2, combs_22[j], eta, N, M, M_opt, mc_samps=MC_SAMPS)
#         _, statistical_power_22_int[j, i], percentiles_22_int[j, i], kurt_flag = statistical_power_mc(N_fock, state, 2, combs_22[j], eta, N, M, M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#         kurt_exceeded_22[j, i] = bool(kurt_flag)

#%%

data_y_22 = np.load('data_files/sub_data_22.npz')['data_y_22']
percentiles_22 = np.load('data_files/sub_data_22.npz')['percentiles_22']
percentiles_22_int = np.load('data_files/sub_data_22.npz')['percentiles_22_int']
statistical_power_22 = np.load('data_files/sub_statistical_powers_22.npz')['statistical_power_22']
statistical_power_22_int = np.load('data_files/sub_statistical_powers_22.npz')['statistical_power_22_int']
kurt_exceeded_22 = np.load('data_files/sub_statistical_powers_22.npz')['kurt_exceeded_22']

#%%

"FIG 1 PHOTON SUBTRACTED TWO-MODE SQUEEZED STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.hlines(0, zs[0], zs[-1], colors='black', zorder=0.2)

linestyles_clt = [(0, (1,3))]*len(combs_22)
linestyles = ['-', (0, (6,1)), '-', 'dashdot']

colors = ['#D55F01','#069E73', '#E69F00','#0072B2']
alpha_black = [0.2,0,0.2,0.2]
labels = [r'D$_{I}$', r'E$_{I}$', r'E$_{II}$', r'E$_{III}$']
cs = [None]*len(combs_22)
for i in [1,2,0,3]:
    cs[i] =ax.fill_between(zs, percentiles_22_int[i,:,3], percentiles_22_int[i,:,1], alpha=0.3,facecolor=colors[::-1][i])

for i in range(len(combs_22)-1,-1,-1):
    ax.plot(zs, data_y_22[i], label=labels[i], linewidth=solidlinewidth, color=cs[i].get_facecolor()[0], alpha=1, ls=linestyles[i])

    
plt.xlabel(r'$\zeta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Determinant',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(loc=(0.2,0.22))
plt.title(r'$M_{{tot}}$ = {} $\eta$ = {} $\bar N$ = {}'.format(M,eta,N))

ax.set_xlim([zs[0],zs[-1]])
ax.set_box_aspect(0.65)

# symmetric log scale for both positive and negative y values
plt.yscale("symlog")

plt.grid(zorder=0)

plt.savefig('subTMSV_fig1.pdf')
plt.savefig('subTMSV_fig1.svg')

plt.show()

#%%

"FIG 2 PHOTON SUBTRACTED TWO-MODE SQUEEZED STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.plot(zs,0.95*np.ones(len(zs)),'k--',label='95%', lw=solidlinewidth)
#colors= ['C0','orange', 'green','C2']
#colors = ['#0072B2','#069E73', '#E69F00', '#D55F01']
colors2 = colors[::-1]

for i in range(len(combs_22)):
    plot_clt(kurt_exceeded_22[i], statistical_power_22_int[i], zs, linestyles[i], colors2[i])


plt.xlabel(r'$\zeta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.xticks(fontsize=12)
labelhandles = [plt.Line2D([0], [0], color=colors2[i], lw=solidlinewidth) for i in range(4)]
labelhandles.append(plt.Line2D([0], [0], color='k', lw=solidlinewidth, ls='--'))
fig.legend(handles=labelhandles, labels=['D$_{I}$', 'E$_{I}$', 'E$_{II}$', 'E$_{III}$', '95%'], loc=(0.5,0.25))
plt.title(r'$M_{{tot}}$ = {} $\eta$ = {} $\bar N$ = {}'.format(M,eta,N))
ax.set_xlim([zs[0],zs[-1]])
ax.set_ylim([0.5,1.02])
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))
ax.set_box_aspect(0.65)

plt.grid()

plt.savefig('subTMSV_fig2.pdf')
plt.savefig('subTMSV_fig2.svg')

plt.show()

#%%

# Data for how statistical power depends on measurement number
Ms = np.round([10**i for i in np.linspace(np.log10(30),4,50)])
#Ms = [10**i for i in np.linspace(2,4,100)]
eta = 0.8
zeta = 1
N = 0


sub_state=sub_vec(N_fock,zeta,1,1)

#%%

# #2x2
# MC_SAMPS = 1000000
# sh = len(Ms)
# data_y_ms_22=np.zeros((len(combs_22), sh))
# statistical_power_ms_22 = np.zeros((len(combs_22), sh))
# statistical_power_ms_22_int = np.zeros((len(combs_22), sh))
# percentiles_ms_22 = np.zeros((len(combs_22), sh, 4)) 
# percentiles_ms_22_int = np.zeros((len(combs_22), sh, 4))


# # Flag: whether kurtosis check ever exceeded threshold for comb [2,4]
# kurt_exceeded_ms_22 = np.zeros((len(combs_22), sh), dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing 2x2 combs'):
#     state = sub_state
#     for j in range(len(combs_22)):
#         ls = combs_22[j]
#         M_opt = optimize_Mij_gradient_descent(N_fock, state, 2, combs_22[j], eta, N, Ms[i], tol=1e-30)[1]
#         M_opt_int = integerisation(Ms[i], M_opt[0] + 1j * M_opt[1], 2, combs_22[j])
#         # compute detector values
#         data_y_ms_22[j, i] = TD_det(N_fock, state, 2, combs_22[j], eta, N)
#         _, statistical_power_ms_22[j, i], percentiles_ms_22[j, i] = statistical_power_mc(N_fock, state, 2, combs_22[j], eta, N, Ms[i], M_opt, mc_samps=MC_SAMPS, integer_measurements=True)
#         _, statistical_power_ms_22_int[j, i], percentiles_ms_22_int[j, i], kurt_flag = statistical_power_mc(N_fock, state, 2, combs_22[j], eta, N, Ms[i], M_opt_int, mc_samps=MC_SAMPS, integer_measurements=True, return_kurtosis=True)
#         kurt_exceeded_ms_22[j, i] = bool(kurt_flag)

statistical_power_ms_22 = np.load('data_files/sub_statistical_powers_ms_22.npz')['statistical_power_ms_22']
statistical_power_ms_22_int = np.load('data_files/sub_statistical_powers_ms_22.npz')['statistical_power_ms_22_int']
kurt_exceeded_ms_22 = np.load('data_files/sub_statistical_powers_ms_22.npz')['kurt_exceeded_ms_22']

#%%

"FIG 3 PHOTON SUBTRACTED TWO-MODE SQUEEZED STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.plot(Ms,0.95*np.ones(len(Ms)),'k--',label='95%', lw=solidlinewidth)

for i in range(len(combs_22)):
    plot_clt(kurt_exceeded_ms_22[i], statistical_power_ms_22_int[i], Ms, linestyles[i], colors2[i])

plt.xlabel(r'$M_{tot}$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.xticks(fontsize=12)
labelhandles = [plt.Line2D([0], [0], color=colors2[i], lw=solidlinewidth) for i in range(4)]
labelhandles.append(plt.Line2D([0], [0], color='k', lw=solidlinewidth, ls='--'))
fig.legend(handles=labelhandles, labels=['D$_{I}$', 'E$_{I}$', 'E$_{II}$', 'E$_{III}$', '95%'], loc=(0.7,0.25))
plt.title(r'$\zeta$ = {} $\eta$ = {} $\bar N$ = {}'.format(zeta,eta,N))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xscale('log')

plt.grid(which='both')
ax.set_xlim([Ms[0],Ms[-1]])
plt.ylim(0.55,1.02)

ax.set_box_aspect(0.65)

plt.savefig('subTMSV_fig3.pdf')
plt.savefig('subTMSV_fig3.svg')

plt.show()

#%%

# Data for how statistical power depends on loss
etas = np.linspace(1, 0.001, 40)
M = 1000
zeta = 1
N = 0

#%%

# #2x2
# MC_SAMPS = 1000000
# sh = len(etas)
# data_y_ks_22=np.zeros((len(combs_22), sh))
# statistical_power_ks_22 = np.zeros((len(combs_22), sh))
# statistical_power_ks_22_int = np.zeros((len(combs_22), sh))
# percentiles_ks_22 = np.zeros((len(combs_22), sh, 4)) 
# percentiles_ks_22_int = np.zeros((len(combs_22), sh, 4))

# # Flag: whether kurtosis check ever exceeded threshold for comb [2,4]
# kurt_exceeded_ks_22 = np.zeros((len(combs_22), sh), dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing 2x2 combs'):
#     state = sub_state
#     for j in range(len(combs_22)):
#         ls = combs_22[j]
#         M_opt = optimize_Mij_gradient_descent(N_fock, state, 2, combs_22[j], etas[i], N, M, tol=1e-30)[1]
#         M_opt_int = integerisation(M, M_opt[0] + 1j * M_opt[1], 2, combs_22[j])
#         # compute detector values
#         data_y_ks_22[j, i] = TD_det(N_fock, state, 2, combs_22[j], etas[i], N)
#         _, statistical_power_ks_22[j, i], percentiles_ks_22[j, i] = statistical_power_mc(N_fock, state, 2, combs_22[j], etas[i], N, M, M_opt, mc_samps=MC_SAMPS, integer_measurements=True)
#         _, statistical_power_ks_22_int[j, i], percentiles_ks_22_int[j, i], kurt_flag = statistical_power_mc(N_fock, state, 2, combs_22[j], etas[i], N, M, M_opt_int, mc_samps=MC_SAMPS, integer_measurements=True, return_kurtosis=True)
#         kurt_exceeded_ks_22[j, i] = bool(kurt_flag)

statistical_power_ks_22 = np.load('data_files/sub_statistical_powers_ks_22.npz')['statistical_power_ks_22']
statistical_power_ks_22_int = np.load('data_files/sub_statistical_powers_ks_22.npz')['statistical_power_ks_22_int']
kurt_exceeded_ks_22 = np.load('data_files/sub_statistical_powers_ks_22.npz')['kurt_exceeded_ks_22']

#%%

"FIG 4 PHOTON SUBTRACTED TWO-MODE SQUEEZED STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))
plt.plot(etas,0.95*np.ones(len(etas)),'k--',label='95%', lw=solidlinewidth)

for i in range(len(combs_22)):
    plot_clt(kurt_exceeded_ks_22[i], statistical_power_ks_22_int[i], etas, linestyles[i], colors2[i])


plt.xlabel(r'$\eta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.yticks(fontsize=12)
legend_elements = [plt.Line2D([0], [0], color=colors2[i], lw=solidlinewidth) for i in range(4)]
legend_elements.append(plt.Line2D([0], [0], color='k', lw=solidlinewidth, ls='--'))
fig.legend(handles=legend_elements, labels=['D$_{I}$', 'E$_{I}$', 'E$_{II}$', 'E$_{III}$', '95%'], loc=(0.3,0.2))
plt.title(r'$M_{{tot}}$ = {} $\zeta$ = {} $\bar N$ = {}'.format(M,zeta,N))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xlim([etas[0],etas[-1]])
ax.set_ylim([0.5,1.02])

plt.grid(which="both")
ax.set_box_aspect(0.65)

plt.savefig('subTMSV_fig4.pdf')
plt.savefig('subTMSV_fig4.svg')

plt.show()


#%%

# save statistical powers, percentiles, and kurtosis flags

# np.savez('data_files/sub_data_22.npz', data_y_22=data_y_22, percentiles_22=percentiles_22, percentiles_22_int=percentiles_22_int)
# np.savez('data_files/sub_statistical_powers_22.npz', statistical_power_22=statistical_power_22, statistical_power_22_int=statistical_power_22_int, kurt_exceeded_22=kurt_exceeded_22)
# np.savez('data_files/sub_statistical_powers_ms_22.npz', statistical_power_ms_22=statistical_power_ms_22, statistical_power_ms_22_int=statistical_power_ms_22_int, kurt_exceeded_ms_22=kurt_exceeded_ms_22)
# np.savez('data_files/sub_statistical_powers_ks_22.npz', statistical_power_ks_22=statistical_power_ks_2２, statistical_power_ks_２２_int=statistical_power_ks_２２_int, kurt_exceeded_ks_２２=kurt_exceeded_ks_２２)

#%%

# print(datetime.now() - start)
# # The mission was succesful and you can now find optimal NPT criteria :) 
# speak("Mission complete. N P T optimised.")