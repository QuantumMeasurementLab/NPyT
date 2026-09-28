#%%
"""
Authors: Lydia A. Kanari-Naish, Amaya Calvo-Sánchez, and Arjun Gupta. Imperial College London
Last update: July 2026

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
Code for optimizing over all suitable NPT criteria for the example of a 
two-mode SCS
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
# from win32com.client import Dispatch
# speak = Dispatch("SAPI.SpVoice").Speak
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

"Define the two-mode SCS"


def state_generator(fock_dims, alpha, beta, gamma, delta, phi):
    #this gives a general 2 mode cat state of the form |alpha>|beta>+exp(i\phi)|gamma>|delta>
    state=(tensor(coherent(fock_dims,alpha),coherent(fock_dims,beta))+np.exp(1j*phi)*tensor(coherent(fock_dims,gamma),coherent(fock_dims,delta))).unit()
    return state

#this creates a vectorised version
vec_states=np.frompyfunc(state_generator,6,1)

#%%

# create range of alphas to consider
alphas = np.linspace(0.0001,2,50) 

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


#%%

# # COMMENTED OUT FOR SPEED

# #2x2
# MC_SAMPS = 1000000
# sh = len(alphas)
# data_y_22=np.zeros((len(combs_22), sh))
# statistical_power_22 = np.zeros((len(combs_22), sh))
# statistical_power_22_int = np.zeros((len(combs_22), sh))
# percentiles_22 = np.zeros((len(combs_22), sh, 4)) 
# percentiles_22_int = np.zeros((len(combs_22), sh, 4))


# # Flag: whether kurtosis check ever exceeded threshold for comb [2,4]
# kurt_exceeded_22 = np.zeros((len(combs_22), sh), dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing 2x2 combs'):
#     state = cat_alphas[i]
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
data_y_22=np.load('data_files/TMSCS_data_22.npz')['data_y_22']
percentiles_22=np.load('data_files/TMSCS_data_22.npz')['percentiles_22']
percentiles_22_int=np.load('data_files/TMSCS_data_22.npz')['percentiles_22_int']
statistical_power_22=np.load('data_files/TMSCS_statistical_powers_22.npz')['statistical_power_22']
statistical_power_22_int=np.load('data_files/TMSCS_statistical_powers_22.npz')['statistical_power_22_int']
kurt_exceeded_22=np.load('data_files/TMSCS_statistical_powers_22.npz')['kurt_exceeded_22']

#%%

# # S3

# MC_SAMPS = 1000000
# sh = len(alphas)

# data_y_S3=np.zeros(sh)
# statistical_power_S3 = np.zeros(sh)
# statistical_power_S3_int = np.zeros(sh)
# percentiles_S3 = np.zeros((sh, 4))
# percentiles_S3_int = np.zeros((sh, 4))

# kurt_exceeded_S3 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc='Processing S3'):
#     state = cat_alphas[i]
#     M_opt_S3 = optimize_Mij_gradient_descent(N_fock, state, 2, [0, 4, 12], eta, N, M, tol=1e-30)[1]
#     M_opt_S3_int = integerisation(M, M_opt_S3[0] + 1j * M_opt_S3[1], 2, [0,4,12], identityincluded=True)
#     data_y_S3[i] = TD_det(N_fock, state, 2, [0, 4, 12], eta, N)
#     _, statistical_power_S3[i], percentiles_S3[i] = statistical_power_mc(N_fock, state, 2, [0, 4, 12], eta, N, M, M_opt_S3, mc_samps=MC_SAMPS)
#     _, statistical_power_S3_int[i], percentiles_S3_int[i], kurt_flag = statistical_power_mc(N_fock, state, 2, [0, 4, 12], eta, N, M, M_opt_S3_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_S3[i] = bool(kurt_flag)

# %%

data_y_S3=np.load('data_files/TMSCS_data_S3.npz')['data_y_S3']
percentiles_S3=np.load('data_files/TMSCS_data_S3.npz')['percentiles_S3']
percentiles_S3_int=np.load('data_files/TMSCS_data_S3.npz')['percentiles_S3_int']
statistical_power_S3=np.load('data_files/TMSCS_statistical_powers_S3.npz')['statistical_power_S3']
statistical_power_S3_int=np.load('data_files/TMSCS_statistical_powers_S3.npz')['statistical_power_S3_int']
kurt_exceeded_S3=np.load('data_files/TMSCS_statistical_powers_S3.npz')['kurt_exceeded_S3']

#%%

"FIG 1 TWO MODE SCHRODINGER CAT STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

labels_22 = ['F$_{I}$', 'F$_{II}$', 'F$_{III}$', 'F$_{IV}$']

colors_22 = ['darkblue', 'gold', 'lightseagreen', 'hotpink']
color_33 = 'mediumorchid'

colors_22 = ['#332288','#DDCC77','#44AA99', '#CC6577']
color_33 = '#882154'

ax.axhline(0,color='black',alpha=0.8, zorder=0.2)



cs = [None]*len(combs_22)
for i in [0,1,2,3]:
    cs[i] =ax.fill_between(alphas, percentiles_22_int[i,:,3], percentiles_22_int[i,:,1], alpha=0.5, facecolor=colors_22[i])
    ax.plot(alphas, data_y_22[i], label=labels_22[i], linewidth=solidlinewidth, color=cs[i].get_facecolor()[0], alpha=1)
    #ax.plot(alphas, data_y_22[i], linewidth=2, color='black', alpha=0.2)

c2 = ax.plot(alphas, data_y_S3,label=r'$S_{III}$', color=color_33, linewidth=solidlinewidth, ls='dashdot')
ax.fill_between(alphas,percentiles_S3_int[:,3],percentiles_S3_int[:,1],alpha=0.4, facecolor=c2[0].get_color())


plt.xlabel(r'$\alpha$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Determinant',fontsize=14)
plt.yticks(fontsize=12)
plt.ylim([-0.25, 0.3])

fig.legend(loc=(0.2,0.55))
plt.title(r'$M_{{tot}}$ = {} $\eta$ = {} $\bar N$ = {}'.format(int(M),eta,N))

ax.set_xlim([alphas[0],alphas[-1]])
ax.set_box_aspect(0.65)
plt.grid(zorder=0)

plt.savefig('TMSCS_fig1.pdf')
plt.savefig('TMSCS_fig1.svg')

plt.show()


#%%

"FIG 2 TWO MODE SCHRODINGER CAT STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.hlines(0.5, alphas[0], alphas[-1], colors='black', zorder=0.2, alpha=0.8)

plt.plot(alphas,0.95*np.ones(len(alphas)),'k--',label=r'95$\%$', lw=solidlinewidth)

for i in range(len(combs_22)):
    plot_clt(kurt_exceeded_22[i], statistical_power_22_int[i], alphas, 'solid', colors_22[i])

plot_clt(kurt_exceeded_S3, statistical_power_S3_int, alphas, 'dashdot', color_33)

plt.ylim([0.35, 1.02])

plt.xlabel(r'$\alpha$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.yticks(fontsize=12)
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1, decimals=None, symbol='%', is_latex=False))

legend_elements = [
    plt.Line2D([0], [0], color=colors_22[i], linestyle='-', label=labels_22[i], lw=solidlinewidth) for i in range(len(combs_22))
]
legend_elements.append(plt.Line2D([0], [0], color=color_33, linestyle='dashdot', label=r'$S_{III}$', lw=solidlinewidth))
legend_elements.append(plt.Line2D([0], [0], color='k', linestyle='--', label=r'95$\%$', lw=solidlinewidth))

fig.legend(handles=legend_elements, loc=(0.4,0.2))
plt.title(r'$M_{{tot}}$ = {} $\eta$ = {} $\bar N$ = {}'.format(int(M),eta,N))
ax.set_xlim([alphas[0],alphas[-1]])
ax.set_box_aspect(0.65)

plt.grid(zorder=0)

plt.savefig('TMSCS_fig2.pdf')
plt.savefig('TMSCS_fig2.svg')

plt.show()

#%%

# Data for how statistical power depends on measurement number

# Plot as a function of M, number of measurement
Ms = np.round([10**i for i in np.linspace(np.log10(30),4,50)])

# Keeping optical coupling fixed
eta = 0.95
# alpha is cat state parameter
alpha = 1
# phonon occupation number of bath
N = 0.01


cat_state=state_generator(N_fock,alpha,0,0,alpha,np.pi)

#%%


# # 2x2
# MC_SAMPS = 1000000
# sh = len(Ms)
# statistical_power_ms_22 = np.zeros((len(combs_22), sh))
# statistical_power_ms_22_int = np.zeros((len(combs_22), sh))


# # Flag: whether kurtosis check ever exceeded threshold for comb [2,4]
# kurt_exceeded_ms_22 = np.zeros((len(combs_22), sh), dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing 2x2 combs'):
#     state = cat_state
#     for j in range(len(combs_22)):
#         ls = combs_22[j]
#         M_opt = optimize_Mij_gradient_descent(N_fock, state, 2, combs_22[j], eta, N, Ms[i], tol=1e-30)[1]
#         M_opt_int = integerisation(Ms[i], M_opt[0] + 1j * M_opt[1], 2, combs_22[j])
#         # compute detector values
#         _, statistical_power_ms_22[j, i], _ = statistical_power_mc(N_fock, state, 2, combs_22[j], eta, N, Ms[i], M_opt, mc_samps=MC_SAMPS)
#         _, statistical_power_ms_22_int[j, i], _, kurt_flag = statistical_power_mc(N_fock, state, 2, combs_22[j], eta, N, Ms[i], M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#         kurt_exceeded_ms_22[j, i] = bool(kurt_flag)

statistical_power_ms_22 = np.load('data_files/TMSCS_statistical_powers_ms_22.npz')['statistical_power_ms_22']
statistical_power_ms_22_int = np.load('data_files/TMSCS_statistical_powers_ms_22.npz')['statistical_power_ms_22_int']
kurt_exceeded_ms_22 = np.load('data_files/TMSCS_statistical_powers_ms_22.npz')['kurt_exceeded_ms_22']

# statistical_power_ms_S3 = np.zeros(sh)
# statistical_power_ms_S3_int = np.zeros(sh)

# kurt_exceeded_ms_S3 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc='Processing S3'):
#     state = cat_state
#     M_opt_S3 = optimize_Mij_gradient_descent(N_fock, state, 2, [0, 4, 12], eta, N, Ms[i], tol=1e-30)[1]
#     M_opt_S3_int = integerisation(Ms[i], M_opt_S3[0] + 1j * M_opt_S3[1], 2, [0, 4, 12], identityincluded=True)
#     _, statistical_power_ms_S3[i], _ = statistical_power_mc(N_fock, state, 2, [0, 4, 12], eta, N, Ms[i], M_opt_S3, mc_samps=MC_SAMPS)
#     _, statistical_power_ms_S3_int[i], _, kurt_flag = statistical_power_mc(N_fock, state, 2, [0, 4, 12], eta, N, Ms[i], M_opt_S3_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_ms_S3[i] = bool(kurt_flag)

statistical_power_ms_S3 = np.load('data_files/TMSCS_statistical_powers_ms_S3.npz')['statistical_power_ms_S3']
statistical_power_ms_S3_int = np.load('data_files/TMSCS_statistical_powers_ms_S3.npz')['statistical_power_ms_S3_int']
kurt_exceeded_ms_S3 = np.load('data_files/TMSCS_statistical_powers_ms_S3.npz')['kurt_exceeded_ms_S3']

#%%

"FIG 3 TWO MODE SCHRODINGER CAT STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.plot(Ms,0.95*np.ones(len(Ms)),'k--',label=r'95$\%$', lw=solidlinewidth)

for i in range(len(combs_22)):
    plot_clt(kurt_exceeded_ms_22[i], statistical_power_ms_22_int[i], Ms, 'solid', colors_22[i])
plot_clt(kurt_exceeded_ms_S3, statistical_power_ms_S3_int, Ms, 'dashdot', color_33)

legend_elements = [
    plt.Line2D([0], [0], color=colors_22[i], linestyle='-', label=labels_22[i], lw=solidlinewidth) for i in range(len(combs_22))
]
legend_elements.append(plt.Line2D([0], [0], color=color_33, linestyle='dashdot', label=r'$S_{III}$', lw=solidlinewidth))
legend_elements.append(plt.Line2D([0], [0], color='k', linestyle='--', label=r'95$\%$', lw=solidlinewidth))

plt.xlabel(r'$M_{tot}$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(handles=legend_elements, loc=(0.7,0.3))
plt.title(r'$\alpha$ = {} $\eta$ = {} $\bar N$ = {}'.format(alpha,eta,N))

ax.set_xscale('log')
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))
ax.set_box_aspect(0.65) # 5/8 aspect ratio

plt.grid(which='both')
#ax.set_xlim([Ms[0],Ms[-1]])
ax.set_xlim(40, Ms[-1])
ax.set_ylim(0.5,1.02)

plt.savefig('TMSCS_fig3.pdf')
plt.savefig('TMSCS_fig3.svg')
plt.show()

#%%

# Data for how statistical power depends on eta
etas = np.round(np.linspace(1,0.001,50), decimals=3)
#etas = np.linspace(1,0.00,100)
M = 1000
alpha = 1
N = 0.01


cat_state=state_generator(N_fock,alpha,0,0,alpha,np.pi)

#%%

# #2x2
# MC_SAMPS = 1000000
# sh = len(etas)
# statistical_power_ks_22 = np.zeros((len(combs_22), sh))
# statistical_power_ks_22_int = np.zeros((len(combs_22), sh))


# # Flag: whether kurtosis check ever exceeded threshold for comb [2,4]
# kurt_exceeded_ks_22 = np.zeros((len(combs_22), sh), dtype=bool)
# for i in tqdm(range(sh), desc=f'Processing 2x2 combs'):
#     state = cat_state
#     for j in range(len(combs_22)):
#         ls = combs_22[j]
#         M_opt = optimize_Mij_gradient_descent(N_fock, state, 2, combs_22[j], etas[i], N, M, tol=1e-30)[1]
#         M_opt_int = integerisation(M, M_opt[0] + 1j * M_opt[1], 2, combs_22[j])
#         # compute detector values
#         _, statistical_power_ks_22[j, i], _ = statistical_power_mc(N_fock, state, 2, combs_22[j], etas[i], N, M, M_opt, mc_samps=MC_SAMPS)
#         _, statistical_power_ks_22_int[j, i], _, kurt_flag = statistical_power_mc(N_fock, state, 2, combs_22[j], etas[i], N, M, M_opt_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#         kurt_exceeded_ks_22[j, i] = bool(kurt_flag)

# statistical_power_ks_S3 = np.zeros(sh)
# statistical_power_ks_S3_int = np.zeros(sh)

# kurt_exceeded_ks_S3 = np.zeros(sh, dtype=bool)
# for i in tqdm(range(sh), desc='Processing S3'):
#     state = cat_state
#     M_opt_S3 = optimize_Mij_gradient_descent(N_fock, state, 2, [0, 4, 12], etas[i], N, M, tol=1e-30)[1]
#     M_opt_S3_int = integerisation(M, M_opt_S3[0] + 1j * M_opt_S3[1], 2, [0, 4, 12], identityincluded=True)
#     _, statistical_power_ks_S3[i], _ = statistical_power_mc(N_fock, state, 2, [0, 4, 12], etas[i], N, M, M_opt_S3, mc_samps=MC_SAMPS)
#     _, statistical_power_ks_S3_int[i], _, kurt_flag = statistical_power_mc(N_fock, state, 2, [0, 4, 12], etas[i], N, M, M_opt_S3_int, mc_samps=MC_SAMPS, return_kurtosis=True)
#     kurt_exceeded_ks_S3[i] = bool(kurt_flag)

statistical_power_ks_22 = np.load('data_files/TMSCS_statistical_powers_ks_22_N01.npz')['statistical_powers_ks_22']
statistical_power_ks_22_int = np.load('data_files/TMSCS_statistical_powers_ks_22_N01.npz')['statistical_powers_ks_22_int']
kurt_exceeded_ks_22 = np.load('data_files/TMSCS_statistical_powers_ks_22_N01.npz')['kurt_exceeded_ks_22']
statistical_power_ks_S3 = np.load('data_files/TMSCS_statistical_powers_ks_S3_N01.npz')['statistical_powers_ks_S3']
statistical_power_ks_S3_int = np.load('data_files/TMSCS_statistical_powers_ks_S3_N01.npz')['statistical_powers_ks_S3_int']
kurt_exceeded_ks_S3 = np.load('data_files/TMSCS_statistical_powers_ks_S3_N01.npz')['kurt_exceeded_ks_S3']

#%%

"FIG 4 TWO MODE SCHRODINGER CAT STATE"

fig, ax = plt.subplots(figsize=(4.6,2.875))

plt.plot(etas,0.95*np.ones(len(etas)),'k--',label=r'95$\%$', lw=solidlinewidth)
plt.hlines(0.5, etas[0], etas[-1], colors='black', zorder=0.2, alpha=0.8)

for i in range(len(combs_22)):
    plot_clt(kurt_exceeded_ks_22[i], statistical_power_ks_22_int[i], etas, 'solid', colors_22[i])
plot_clt(kurt_exceeded_ks_S3, statistical_power_ks_S3_int, etas, 'dashdot', color_33)

legend_elements = [
    plt.Line2D([0], [0], color=colors_22[i], linestyle='-', label=labels_22[i], lw=solidlinewidth) for i in range(len(combs_22))
]
legend_elements.append(plt.Line2D([0], [0], color=color_33, linestyle='dashdot', label=r'$S_{III}$', lw=solidlinewidth))
legend_elements.append(plt.Line2D([0], [0], color='k', linestyle='--', label=r'95$\%$', lw=solidlinewidth))

plt.xlabel(r'$\eta$',fontsize=14)
plt.xticks(fontsize=12)
plt.ylabel('Statistical power',fontsize=14)
plt.yticks(fontsize=12)
fig.legend(handles=legend_elements, loc=(0.25,0.2))
plt.title(r'$M_{{tot}}$ = {} $\alpha$ = {} $\bar N$ = {}'.format(int(M),alpha,N))
ax.yaxis.set_major_formatter(mtick.PercentFormatter(xmax=1.0, decimals=0, symbol='%', is_latex=False))

ax.set_xlim([etas[0],etas[-1]])
ax.set_ylim([0.0,1.05])
ax.set_box_aspect(0.65)

ax.set_yticks(np.arange(0,1.1,0.1))

plt.grid(which="both")

plt.savefig('TMSCS_fig4.pdf')
plt.savefig('TMSCS_fig4.svg')

plt.show()

#%%

# np.savez('data_files/TMSCS_data_22.npz', data_y_22=data_y_22, percentiles_22=percentiles_22, percentiles_22_int=percentiles_22_int)
# np.savez('data_files/TMSCS_data_S3.npz', data_y_S3=data_y_S3, percentiles_S3=percentiles_S3, percentiles_S3_int=percentiles_S3_int)
# np.savez('data_files/TMSCS_statistical_powers_22.npz', statistical_powers_22=statistical_power_22, statistical_powers_22_int=statistical_power_22_int, kurt_exceeded_22=kurt_exceeded_22)
# np.savez('data_files/TMSCS_statistical_powers_S3.npz', statistical_powers_S3=statistical_power_S3, statistical_powers_S3_int=statistical_power_S3_int, kurt_exceeded_S3=kurt_exceeded_S3)
# np.savez('data_files/TMSCS_statistical_powers_ms_22.npz', statistical_powers_ms_22=statistical_power_ms_22, statistical_powers_ms_22_int=statistical_power_ms_22_int, kurt_exceeded_ms_22=kurt_exceeded_ms_2２)
# np.savez('data_files/TMSCS_statistical_powers_ms_S3.npz', statistical_powers_ms_S3=statistical_power_ms_S3, statistical_powers_ms_S3_int=statistical_power_ms_S3_int, kurt_exceeded_ms_S3=kurt_exceeded_ms_S3)
# np.savez('data_files/TMSCS_statistical_powers_ks_22_N01.npz', statistical_powers_ks_22=statistical_power_ks_22, statistical_powers_ks_2₂_int=statistical_power_ks_₂_int, kurt_exceeded_ks_₂=kurt_exceeded_ks_₂)
# np.savez('data_files/TMSCS_statistical_powers_ks_S3_N01.npz', statistical_powers_ks_S3=statistical_power_ks_S3, statistical_powers_ks_S3_int=statistical_power_ks_S3_int, kurt_exceeded_ks_S3=kurt_exceeded_ks_S3)

#%%

# print(datetime.now() - start)
# # The mission was succesful and you can now find optimal NPT criteria :) 
# speak("Mission complete. N P T optimised.")

#%%

