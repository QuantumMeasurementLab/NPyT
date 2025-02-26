"""
(N)PT test optimization (Py)thon (T)oolbox
"""
#___________________________________________________________________________________ 
#
#                 _   ______       ______
#                / | / / __ \__  _/_  __/
#               /  |/ / /_/ / / / // /   
#              / /|  / ____/ /_/ // /    
#             /_/ |_/_/    \__, //_/     
#                         /____/         
# ___________________________________________________________________________________
# ___________________________________________________________________________________                         
 #  _  _ ___ _____   _          _              _   _       _         _   _          
 # | \| | _ \_   _| | |_ ___ __| |_   ___ _ __| |_(_)_ __ (_)_____ _| |_(_)___ _ _  
 # | .` |  _/ | |   |  _/ -_|_-<  _| / _ \ '_ \  _| | '  \| |_ / _` |  _| / _ \ ' \ 
 # |_|\_|_|  _|_|_   \__\___/__/\__|_\___/ .__/\__|_|_|_|_|_/__\__,_|\__|_\___/_||_|
 # | _ \_  _| |_| |_  ___ _ _   |_   _|__|_|__| | |__  _____ __                     
 # |  _/ || |  _| ' \/ _ \ ' \    | |/ _ \/ _ \ | '_ \/ _ \ \ /                     
 # |_|  \_, |\__|_||_\___/_||_|   |_|\___/\___/_|_.__/\___/_\_\                     
 #      |__/                                                                                                             
# ___________________________________________________________________________________                           
                                                                 
"""
Author: Lydia A. Kanari-Naish, Imperial College London
Last update: September 2023

For the research paper:
"Optimizing confidence in negative-partial-transpose-based 
 entanglement criteria"

Code for optimizing over all suitable NPT criteria for a given state,
with examples of two-mode squeezed vacuum (TMSV) state, the photon 
subtracted/ added TMSV, and the two mode Schrodinger cat state given in
seperate .py files
"""

#!/usr/bin/env python
# coding: utf-8


from functools import reduce
from itertools import combinations
import numpy as np
from qutip import *
from scipy.special import comb as comb
import math
from scipy.stats import norm
from scipy import optimize


"Fock space dimension"
#turn up/down as required
N_fock=20;


#%%
"Functions in the toolbox"

# In order to generate the index label of pqrs, nmkl we use the ordering outlined in Section II A. 
# We have chosen to group according to length=number of indices and total=p+q+r+s.
def get_tuples(length, total):
    if length == 1:
        yield (total,)
        return
    for i in range(total + 1):
        for t in get_tuples(length - 1, total - i):
            yield t+(i,)


#Note the order can only be even
#so order=1->order=2 in the paper,
# order=2-> order =4 in the paper, and so on
def full_matrix(order):
    k=np.array(reduce(lambda x,y:x+y, [list(get_tuples(4, i)) for i in range(order+1)]))
    pqrs=k[:,np.array([1,0,2,3])]*np.array(['a','b','c','d'], object)
    nmkl=k[:,np.array([0,1,3,2])]*np.array(['a','b','c','d'],object)
    pq=pqrs[:,:2]
    rs=pqrs[:,2:]
    nm=nmkl[:,:2]
    kl=nmkl[:,2:]
    amat=np.sum(pq,axis=-1)[...,None]+np.sum(nm,axis=-1)[None,...]
    bmat=np.sum(kl,axis=-1)[None,...]+np.sum(rs,axis=-1)[...,None]
#     pqrs=np.sum(k[:,np.array([1,0,2,3])]*np.array(['a','b','c','d'], object),axis=-1)
#     nmkl=np.sum(k[:,np.array([0,1,3,2])]*np.array(['a','b','c','d'],object),axis=-1)
#     print(k[:,np.array([1,0,2,3])]*np.array(['a','b','c','d'], object))
#     return pqrs[...,None]+nmkl[None,...]
    out=amat+bmat
    out[0,0]='I'
    return out



def sub_matrix(order, rows):
    # Rows selects which rows and columns we keep. The rows and columns are chosen in a pairwise fashion.
    # e.g. rows=[0,1,2] selects 0th, 1st, 2nd columns/rows
    # Note the difference (-1) between python indexing and the indexing used in the paper.
    # e.g. DI is sub_matrix(1,[2,4]), which is d=2, n=2, rows=(3,5) in the paper
                   
    take=np.array(rows)
    matrix=full_matrix(order)
    
#     limit=matrix.shape[0]
#     if rows[-1]>limit:
#         print('Rows are outside range of matrix, last row index must be less than or equal to {}'.format(limit))
    
    return matrix[take,:][:,take]




# To convert strings to qutip operators.
def string_to_op(s,fock_dims,):
    op_dict = {
      'a': tensor(create(fock_dims),qeye(fock_dims)),
      'b': tensor(destroy(fock_dims),qeye(fock_dims)),
      'c': tensor(qeye(fock_dims),create(fock_dims)),
      'd': tensor(qeye(fock_dims),destroy(fock_dims)),
      'I' : tensor(qeye(fock_dims),qeye(fock_dims))
    }
    return reduce(lambda x,y:x*y, [op_dict[i] for i in list(s)])




# Vectorizes function string_to_op.
string_convert=np.frompyfunc(string_to_op,2,1)




# Converts submatrix of strings to submatrix of qutip operators.
def sub_matrix_op(order,rows,fock_dims):
    k=sub_matrix(order,rows)
#     k[0,0]='I'
    return string_convert(k,fock_dims)




# Calculates expectation value of each entry in matrix.
# fock_dims is fock dimensions, must match fock dimensions of state
# state can be an array of states

def matrix_det(order, rows, fock_dims, state):
    #submatrix of operators
    matrix_op=sub_matrix_op(order,rows,fock_dims)
    
    #expect function from Qutip into numpy
    expectfn=np.frompyfunc(expect,2,1)
    state=np.array(state,object)
    nd=state.ndim
    p=expectfn(np.expand_dims(matrix_op,[i for i in range(nd)]), state[...,None,None]).astype(complex)

    return np.real(np.linalg.det(p))



def generate_index_combinations(order, submat_size):
    k=np.array(reduce(lambda x,y:x+y, [list(get_tuples(4, i)) for i in range(order+1)]),dtype=object)
    limit=k.shape[0]
    all_subs=np.array(list(combinations(range(limit),submat_size)))
    return all_subs



# A test of the submatrix function

sub_matrix(2,[12,14])



#to get initial max value of n for order 4 dims 2x2



# A preliminary search, which is done:
# (i) on the pure state to identify any determinants that are negative in any region of parameter space,
# (ii) for a certain matrix size and max order.

def my_state(order,submat_size,fock_dims,state):
    subs=generate_index_combinations(order, submat_size)
    dets=[]
    combs=[]

    #states=vec_states(fock_dims,state)

    for i,c in enumerate(subs):
        det=matrix_det(order, c, fock_dims, state)
        if (det<-1e-10).any():
            dets.append(det)
            combs.append(c)
    
    return np.array(dets),np.array(combs)




# This function gives the factor that arises from normally ordering two operators,
# e.g. a a^dag reordered to give a^\dag a for arbitrary powers of each operator
# cf Eq. 16 from E. Shchukin and W. Vogel PRL 95, 230502 (2005).
def G(n,m,k):
    num = math.factorial(n)*math.factorial(m)
    den = math.factorial(k)*(math.factorial(n-k))*(math.factorial(m-k))
    return num/den



# kronecker delta
def d(i,j):
    return 0 if i!=j else 1




def qutip_ops(fock_dims):
    
    a = tensor(create(fock_dims),qeye(fock_dims))
    b = tensor(destroy(fock_dims),qeye(fock_dims))
    c = tensor(qeye(fock_dims),create(fock_dims))
    d = tensor(qeye(fock_dims),destroy(fock_dims))

    return a, b, c, d



#a, b, c, d, = qutip_ops(20)
a, b, c, d, = qutip_ops(N_fock)




# Finds the expectation value of operator <a^dag n a^m b^dag k b^l (t)>
# 0 < eta < 1
def H(dims, state, n, m, k, l, eta , N):
    #kappaT = -np.log(np.sqrt(eta))
    a, b, c, d, = qutip_ops(dims)
    ans = 0
    if eta==0:
        for p in range(0, min(n,m)+1):
            for r in range(0, min(k,l)+1):
                res = comb(n,p)*comb(m,p)*comb(k,r)*comb(l,r)
                res*=(0**(n+m+k+l-2*p-2*r))
                res*= math.factorial(p)*math.factorial(r)
                res*=((N*(1-0))**p)
                res*=((N*(1-0))**r)
                    #qu_exp = expect(power(dims,a,n)*power(dims,b,m)*power(dims,c,k)*power(dims,d,l), state)
                qu_exp = expect((a**(n-p))*(b**(m-p))*(c**(k-r))*(d**(l-r)), state)
                res*=qu_exp
                ans += res
    else: 
        kappaT = -np.log(np.sqrt(eta))
        a, b, c, d, = qutip_ops(dims)
        ans = 0
        for p in range(0, min(n,m)+1):
            for r in range(0, min(k,l)+1):
                res = comb(n,p)*comb(m,p)*comb(k,r)*comb(l,r)
                res*=(np.exp(-kappaT)**(n+m+k+l-2*p-2*r))
                res*= math.factorial(p)*math.factorial(r)
                res*=((N*(1-np.exp(-2*kappaT)))**p)
                res*=((N*(1-np.exp(-2*kappaT)))**r)
                    #qu_exp = expect(power(dims,a,n)*power(dims,b,m)*power(dims,c,k)*power(dims,d,l), state)
                qu_exp = expect((a**(n-p))*(b**(m-p))*(c**(k-r))*(d**(l-r)), state)
                res*=qu_exp
                ans += res
                    
    return ans




def string_to_ind(string):
    out = np.zeros(8,dtype=int)
    if string=='I':
        return out
    curr_character = string[0]
    ind_dict = {'a':0,'b':1,'c':4,'d':5}
    pairs={'b':'a','d':'c'}
    i=0
    while i<len(string):
        char=string[i]
        if char!=curr_character and char not in pairs:
            ind_dict[curr_character]+=2
        if curr_character in pairs:
            ind_dict[pairs[curr_character]]+=2
            pairs.pop(curr_character)
        out[ind_dict[char]]+=1
        curr_character = char
        i+=1
    return out
#[string_to_ind(i) for i in case_s]




# Time dependent determinant
def TD_det(dims, state, order, combs, eta, N):
    mat_str = sub_matrix(order,combs)
    sh = len(combs)
    mat_ind = np.empty((sh,sh),dtype=complex)
    for i in range(sh):
        for j in range(sh):
            n,m,k,l,p,q,r,s = tuple(string_to_ind(mat_str[i][j]))
            mom = 0.
            for f in range(min(m,k)+1):
                for g in range(min(q,r)+1):
                    mom+=G(m,k,f)*G(q,r,g)*H(dims,state,n+k-f,m-f+l,p+r-g,q-g+s,eta,N)
#                     mom+=G(m,k,f)*G(q,r,g)*H_vec(dims,state,n+k-f,m-f+l,p+r-g,q-g+s,kappa,t,N)
            mat_ind[i][j]=mom
    return np.real(np.linalg.det(mat_ind))



# Time dependent matrix
def TD_mat(dims, state, order, combs, eta, N):
    mat_str = sub_matrix(order,combs)
    sh = len(combs)
    mat_ind = np.empty((sh,sh),dtype=complex)
    for i in range(sh):
        for j in range(sh):
            n,m,k,l,p,q,r,s = tuple(string_to_ind(mat_str[i][j]))
            mom = 0.
            for f in range(min(m,k)+1):
                for g in range(min(q,r)+1):
                    mom+=G(m,k,f)*G(q,r,g)*H(dims,state,n+k-f,m-f+l,p+r-g,q-g+s,eta,N)
            mat_ind[i][j]=mom
    return mat_ind




# matrix adjugate
def adj(dims, state, order, combs, eta, N):
    A = TD_mat(dims, state, order, combs, eta, N)
    
    dim1 = A.shape[-1]
    dim2 = A.shape[-2]
    out=np.zeros_like(A)
    
    for i in range(dim1):
        for j in range(dim2):
            slice1=np.concatenate((np.arange(i),np.arange(i+1,dim1)))
            slice2=np.concatenate((np.arange(j),np.arange(j+1,dim2)))
            submat=A[...,slice1,:][...,:,slice2]
            out[...,i,j] = ((-1)**(i+j)) * np.linalg.det(submat)
    return out.T



# Variance of an operator in the form a^dag n a^m a^dag k a^l b^dag p b^q b^dag r b^s
def full_var(dims, state, order, combs, eta, N):
    sh = len(combs)
    X_real = np.zeros((sh,sh),dtype=complex)
    X_im = np.zeros((sh,sh),dtype=complex)
    mat_str = sub_matrix(order,combs)
    
    # i, j is index of A matrix
    for i in range(sh):
        for j in range(i,sh):
            n,m,k,l,p,q,r,s = tuple(string_to_ind(mat_str[i][j]))

            #this handles the decomposition of moments into 2 Hermitian operators
            # A operator, real operator
            exp_AA = 0.
            for f in range(min(m,k)+1):
                for g in range(min(q,r)+1):
                    for u in range(min(m,k)+1):
                        for v in range(min(q,r)+1):
                            # actual expectation values
                            for x in range(min(m+l-f,n+k-u)+1):
                                for y in range(min(q+s-g,p+r-v)+1):
                                    res = G(m+l-f,n+k-u,x)*G(q+s-g,p+r-v,y)
                                    res *= H(dims, state,(n+k-f+n+k-u-x),(m+l-f+m+l-u-x),(p+r-g+p+r-v-y),(q+s-g+q+s-v-y),eta,N)
                                    exp_AA += res
#                                     print(i,j,(n+k-f+n+k-u-x)(m+l-f+m+l-u-x),(p+r-g+p+r-v-y),(q+s-g+q+s-v-y),res)
                            for x in range(min(m+l-f,m+l-u)+1):
                                for y in range(min(q+s-g,q+s-v)+1):
                                    res = G(m+l-f,m+l-u,x)*G(q+s-g,q+s-v,y)
                                    res *= H(dims, state,(n+k-f+m+l-u-x),(m+l-f+n+k-u-x),(p+r-g+q+s-v-y),(q+s-g+p+r-v-y),eta,N)
                                    exp_AA += res
#                                     print((n+k-f+m+l-u-x),(m+l-f+n+k-u-x),(p+r-g+q+s-v-y),(q+s-g+p+r-v-y),res)
                            for x in range(min(n+k-f,n+k-u)+1):
                                for y in range(min(p+r-g,p+r-v)+1):
                                    res = G(n+k-f,n+k-u,x)*G(p+r-g,p+r-v,y)
                                    res *= H(dims, state,(m+l-f+n+k-u-x),(n+k-f+m+l-u-x),(q+s-g+p+r-v-y),(p+r-g+q+s-v-y),eta,N)
                                    exp_AA += res
#                                     print((m+l-f+n+k-u-x),(n+k-f+m+l-u-x),(q+s-g+p+r-v-y),(p+r-g+q+s-v-y),res)
                            for x in range(min(n+k-f,m+l-u)+1):
                                for y in range(min(p+r-g,q+s-v)+1):
                                    res = G(n+k-f,m+l-u,x)*G(p+r-g,q+s-v,y)
                                    res *= H(dims, state,(m+l-f+m+l-u-x),(n+k-f+n+k-u-x),(q+s-g+q+s-v-y),(p+r-g+p+r-v-y),eta,N)
                                    exp_AA += res
#                                     print((m+l-f+m+l-u-x),(n+k-f+n+k-u-x),(q+s-g+q+s-v-y),(p+r-g+p+r-v-y),res)
                            exp_AA*= G(m,k,f)*G(q,r,g)*G(m,k,u)*G(q,r,v)
            exp_AA*=0.25
            
            exp_A = 0.
            for f in range(min(m,k)+1):
                for g in range(min(q,r)+1):
                    res = H(dims, state,(n+k-f),(m+l-f),(p+r-g),(q+s-g),eta,N)+H(dims,state,(m+l-f),(n+k-f),(q+s-g),(r+p-g),eta,N)
                    res *= 0.5*G(m,k,f)*G(q,r,g)
                    exp_A+= res
            var_A = exp_AA - exp_A**2
            X_real[i][j]=var_A
            
            if i!=j:
            
                exp_BB = 0.
                for f in range(min(m,k)+1):
                    for g in range(min(q,r)+1):
                        for u in range(min(m,k)+1):
                            for v in range(min(q,r)+1):
                                # actual expectation values
                                for x in range(min(m+l-f,n+k-u)+1):
                                    for y in range(min(q+s-g,p+r-v)+1):
                                        res = G(m+l-f,n+k-u,x)*G(q+s-g,p+r-v,y)
                                        res *= H(dims, state,(n+k-f+n+k-u-x),(m+l-f+m+l-u-x),(p+r-g+p+r-v-y),(q+s-g+q+s-v-y),eta,N)
                                        exp_BB += res
                                for x in range(min(m+l-f,m+l-u)+1):
                                    for y in range(min(q+s-g,q+s-v)+1):
                                        res = G(m+l-f,m+l-u,x)*G(q+s-g,q+s-v,y)
                                        res *= H(dims, state,(n+k-f+m+l-u-x),(m+l-f+n+k-u-x),(p+r-g+q+s-v-y),(q+s-g+p+r-v-y),eta,N)
                                        exp_BB -= res
                                for x in range(min(n+k-f,n+k-u)+1):
                                    for y in range(min(p+r-g,p+r-v)+1):
                                        res = G(n+k-f,n+k-u,x)*G(p+r-g,p+r-v,y)
                                        res *= H(dims, state,(m+l-f+n+k-u-x),(n+k-f+m+l-u-x),(q+s-g+p+r-v-y),(p+r-g+q+s-v-y),eta,N)
                                        exp_BB -= res
                                for x in range(min(n+k-f,m+l-u)+1):
                                    for y in range(min(p+r-g,q+s-v)+1):
                                        res = G(n+k-f,m+l-u,x)*G(p+r-g,q+s-v,y)
                                        res *= H(dims, state,(m+l-f+m+l-u-x),(n+k-f+n+k-u-x),(q+s-g+q+s-v-y),(p+r-g+p+r-v-y),eta,N)
                                        exp_BB += res
                                exp_BB*= G(m,k,f)*G(q,r,g)*G(m,k,u)*G(q,r,v)
                exp_BB*=-0.25

                exp_B = 0.
                for f in range(min(m,k)+1):
                    for g in range(min(q,r)+1):
                        res = H(dims, state,(n+k-f),(m+l-f),(p+r-g),(q+s-g),eta,N)-H(dims,state,(m+l-f),(n+k-f),(q+s-g),(r+p-g),eta,N)
                        res *= -1j*0.5*G(m,k,f)*G(q,r,g)
                        exp_B+= res

                var_B =exp_BB - exp_B**2
                X_im[i][j]=var_B


    return [X_real,X_im]



# This is the capital Gamma in the paper. The overall error is Gamma/\sqrt{M}
def Gamma(dims, state, order, combs, eta, N):
    sh = len(combs)
    Adj = adj(dims, state, order, combs, eta, N)
    re_adj = np.real(Adj)
    im_adj = np.imag(Adj)
    X = full_var(dims, state, order, combs, eta, N)
    if Adj is None:
        return np.nan
    
    else:
        gamma = 0
        for i in range(sh):
            gamma+=np.abs(re_adj[i][i])*np.sqrt(X[0][i][i])

        for i in range(sh):
            for j in range(i+1,sh):
                gamma+= 2*(np.abs(re_adj[i][j])*np.sqrt(X[0][i][j])+np.abs(im_adj[i][j])*np.sqrt(X[1][i][j]))

        return np.real(gamma)



def sample_matrix(dims, state, order, combs, eta, N, n_samps, M):
    sh = len(combs)
    Adj = adj(dims, state, order, combs, eta, N)
    re_adj = np.real(Adj)
    im_adj = np.imag(Adj)
    X_r, X_i = full_var(dims, state, order, combs, eta, N)
    true_mat = TD_mat(dims, state, order, combs, eta, N)
    rng = np.random.default_rng()
    out = np.zeros((n_samps, sh, sh))
    gamma = Gamma(dims, state, order, combs, eta, N)
    smu = gamma/M
    for i in range(sh):
        m = np.abs(np.real(Adj[i][i]))*np.sqrt(abs(np.real(X_r[i][i])))/smu
        if m==0:
            m=1
        out[:,i,i]=rng.normal(np.real(true_mat[i][i]),np.sqrt(abs(np.real(X_r[i][i])/m)),size=n_samps)
    for i in range(sh):
        for j in range(i+1,sh):
            m_r = 2*np.abs(np.real(Adj[i][j]))*np.sqrt(abs(np.real(X_r[i][j])))/smu
            m_i = 2*np.abs(np.imag(Adj[i][j]))*np.sqrt(abs(np.real(X_i[i][j])))/smu
            if m_i==0:
                m_i=1
            if m_r==0:
                m_r=1
            out[:,i,j]=rng.normal(np.real(true_mat[i][j]),np.sqrt(abs(np.real(X_r[i][j])/m_r)),size=n_samps)+1j*rng.normal(np.imag(true_mat[i][j]),np.sqrt(abs(np.real(X_i[i][j])/m_i)),size=n_samps)
            out[:,j,i]=np.conj(out[:,i,j])
    return out

def sample_determinant(dims, state, order, combs, eta, N, n_samps=10000, M=1000):
    matrix_samples = sample_matrix(dims, state, order, combs, eta, N, n_samps, M)
#     print(matrix_samples)
    return np.linalg.det(matrix_samples)




def error(dims, state, order, combs, eta, N, M):
    if Gamma == np.nan:
        return np.nan
    else:
        return Gamma(dims, state, order, combs, eta, N)/np.sqrt(M)



def confidence_level(y, err):
    rat=y/err
    rat[np.isnan(rat)]=0 
    return 1-norm.cdf(rat, loc=0, scale=1)


