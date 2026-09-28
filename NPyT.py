"""
NPT test optimization Python Toolbox
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
Authors: Lydia A. Kanari-Naish, Amaya Calvo-Sánchez, and Arjun Gupta. Imperial College London
Last update: July 2026

For the research paper:
"Optimizing the statistical power of negative-partial-transpose-based 
 entanglement tests"

Code for optimizing over all suitable NPT criteria for a given state,
with examples of two-mode squeezed vacuum (TMSV) state, the photon 
subtracted/ added TMSV, and the two mode Schrodinger cat state given in
separate .py files
"""

#!/usr/bin/env python
# coding: utf-8


from functools import reduce
from itertools import combinations
from itertools import cycle
import numpy as np
from qutip import *
from scipy.special import comb as comb
import math
from scipy.stats import norm
from scipy import optimize
import sparse
from typing import Tuple, Generator
import string


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
# order=2-> order=4 in the paper, and so on
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



def is_positive_definite(operatorstr):
    return operatorstr.count('a')==operatorstr.count('b') and operatorstr.count('c')==operatorstr.count('d')

def positivedefinite_count(order, rows):
    submatrix = sub_matrix(order, rows)
    return np.vectorize(is_positive_definite)(submatrix)
    


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


# second-order adjugate: determinants from removing 2 arbitrary rows and 2 arbitrary columns
# returns a 4D array. Indices [i1,i2,j1,j2] correspond to removing rows i1,i2 and columns j1,j2
def adj2(dims, state, order, combs, eta, N):
    A = TD_mat(dims, state, order, combs, eta, N)
    
    sh = A.shape[-1]  # Assuming square matrix
    out = np.zeros((sh, sh, sh, sh), dtype=complex)
    
    # Iterate over all pairs of rows to remove (all combinations, not just i1 < i2)
    for i1 in range(sh):
        for i2 in range(sh):
            if i1 == i2:  # Skip if rows are the same
                continue
            # Iterate over all pairs of columns to remove (all combinations, not just j1 < j2)
            for j1 in range(sh):
                for j2 in range(sh):
                    if j1 == j2:  # Skip if columns are the same
                        continue
                    # Create index arrays excluding i1, i2 and j1, j2
                    row_indices = [i for i in range(sh) if i != i1 and i != i2]
                    col_indices = [j for j in range(sh) if j != j1 and j != j2]
                    
                    # Extract the (sh-2) x (sh-2) submatrix
                    submat = A[np.ix_(row_indices, col_indices)]
                    
                    # Compute determinant
                    out[i1, i2, j1, j2] = np.linalg.det(submat)
    
    return out


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
                            exp_AA_term = 0.
                            for x in range(min(m+l-f,n+k-u)+1):
                                for y in range(min(q+s-g,p+r-v)+1):
                                    res = G(m+l-f,n+k-u,x)*G(q+s-g,p+r-v,y)
                                    res *= H(dims, state,(n+k-f+n+k-u-x),(m+l-f+m+l-u-x),(p+r-g+p+r-v-y),(q+s-g+q+s-v-y),eta,N)
                                    exp_AA_term += res
#                                     print(i,j,(n+k-f+n+k-u-x)(m+l-f+m+l-u-x),(p+r-g+p+r-v-y),(q+s-g+q+s-v-y),res)
                            for x in range(min(m+l-f,m+l-u)+1):
                                for y in range(min(q+s-g,q+s-v)+1):
                                    res = G(m+l-f,m+l-u,x)*G(q+s-g,q+s-v,y)
                                    res *= H(dims, state,(n+k-f+m+l-u-x),(m+l-f+n+k-u-x),(p+r-g+q+s-v-y),(q+s-g+p+r-v-y),eta,N)
                                    exp_AA_term += res
#                                     print((n+k-f+m+l-u-x),(m+l-f+n+k-u-x),(p+r-g+q+s-v-y),(q+s-g+p+r-v-y),res)
                            for x in range(min(n+k-f,n+k-u)+1):
                                for y in range(min(p+r-g,p+r-v)+1):
                                    res = G(n+k-f,n+k-u,x)*G(p+r-g,p+r-v,y)
                                    res *= H(dims, state,(m+l-f+n+k-u-x),(n+k-f+m+l-u-x),(q+s-g+p+r-v-y),(p+r-g+q+s-v-y),eta,N)
                                    exp_AA_term += res
#                                     print((m+l-f+n+k-u-x),(n+k-f+m+l-u-x),(q+s-g+p+r-v-y),(p+r-g+q+s-v-y),res)
                            for x in range(min(n+k-f,m+l-u)+1):
                                for y in range(min(p+r-g,q+s-v)+1):
                                    res = G(n+k-f,m+l-u,x)*G(p+r-g,q+s-v,y)
                                    res *= H(dims, state,(m+l-f+m+l-u-x),(n+k-f+n+k-u-x),(q+s-g+q+s-v-y),(p+r-g+p+r-v-y),eta,N)
                                    exp_AA_term += res
#                                     print((m+l-f+m+l-u-x),(n+k-f+n+k-u-x),(q+s-g+q+s-v-y),(p+r-g+p+r-v-y),res)
                            exp_AA += exp_AA_term * G(m,k,f)*G(q,r,g)*G(m,k,u)*G(q,r,v)
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
                                exp_BB_term = 0.
                                # actual expectation values
                                for x in range(min(m+l-f,n+k-u)+1):
                                    for y in range(min(q+s-g,p+r-v)+1):
                                        res = G(m+l-f,n+k-u,x)*G(q+s-g,p+r-v,y)
                                        res *= H(dims, state,(n+k-f+n+k-u-x),(m+l-f+m+l-u-x),(p+r-g+p+r-v-y),(q+s-g+q+s-v-y),eta,N)
                                        exp_BB_term += res
                                for x in range(min(m+l-f,m+l-u)+1):
                                    for y in range(min(q+s-g,q+s-v)+1):
                                        res = G(m+l-f,m+l-u,x)*G(q+s-g,q+s-v,y)
                                        res *= H(dims, state,(n+k-f+m+l-u-x),(m+l-f+n+k-u-x),(p+r-g+q+s-v-y),(q+s-g+p+r-v-y),eta,N)
                                        exp_BB_term -= res
                                for x in range(min(n+k-f,n+k-u)+1):
                                    for y in range(min(p+r-g,p+r-v)+1):
                                        res = G(n+k-f,n+k-u,x)*G(p+r-g,p+r-v,y)
                                        res *= H(dims, state,(m+l-f+n+k-u-x),(n+k-f+m+l-u-x),(q+s-g+p+r-v-y),(p+r-g+q+s-v-y),eta,N)
                                        exp_BB_term -= res
                                for x in range(min(n+k-f,m+l-u)+1):
                                    for y in range(min(p+r-g,q+s-v)+1):
                                        res = G(n+k-f,m+l-u,x)*G(p+r-g,q+s-v,y)
                                        res *= H(dims, state,(m+l-f+m+l-u-x),(n+k-f+n+k-u-x),(q+s-g+q+s-v-y),(p+r-g+p+r-v-y),eta,N)
                                        exp_BB_term += res
                                exp_BB += exp_BB_term * G(m,k,f)*G(q,r,g)*G(m,k,u)*G(q,r,v)
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


# This function takes the variance and divides by the number of measurements Mij,p to give the weighted variance
def weighted_var(Mij,var,combs):
    # Mij shape should be
    # (2, len(combs), len(combs)) where Mij[0] is the weighting for the real part of the variance and Mij[1] is the weighting for the imaginary part of the variance
    if Mij.shape != (2, len(combs), len(combs)):
        raise ValueError("Mij must have shape (2, len(combs), len(combs))")
    # Make a copy to avoid modifying the input array
    var_copy = np.array(var, copy=True)
    var_copy[0] = (var_copy[0] + var_copy[0].T) 
    var_copy[1] = (var_copy[1] + var_copy[1].T)
    var_copy[0][np.diag_indices_from(var_copy[0])] /= 2
    weighted_var = np.zeros_like(var_copy)
    weighted_var[0] = var_copy[0] / Mij[0]
    # For Mij[1], handle the case where diagonal entries are 0 (constrained, no imaginary measurements on diagonal)
    weighted_var[1] = var_copy[1] / Mij[1]
    weighted_var[1][np.diag_indices(len(combs))] =0
    if 0 in combs:
        weighted_var[0][0][0] = 0
    return weighted_var


# This function calculates the error in the determinant to second order
def delta_detA2_secondorder(Mij,var, Adj1, Adj2, combs, first_order=False):
    # Adj1 = adj(dims, state, order, combs, eta, N)
    # Adj2 = adj2(dims, state, order, combs, eta, N)
    # if Mijbool:
    #     X = weighted_var(combs, Mij,var)
    # else:
    #     X = full_var(dims, state, order, combs, eta, N)
    X = weighted_var(Mij,var,combs)

    re_adj1 = np.real(Adj1)
    im_adj1 = np.imag(Adj1)
    re_adj2 = np.real(Adj2)
    im_adj2 = np.imag(Adj2)

    delta = 0

    sh = len(combs)

    # first order

    for i1 in range(sh):
        delta += re_adj1[i1][i1]**2 * abs(np.real(X[0][i1][i1]))
    for i1 in range(sh):
        for j1 in range(i1+1, sh):
            delta += 2*re_adj1[i1][j1]**2*abs(np.real(X[0][i1][j1]))
            delta += 2*im_adj1[i1][j1]**2*abs(np.real(X[1][i1][j1]))

    # second order

    if first_order: 
        return delta

    for i1 in range(sh):
        for j1 in range(sh):
            for i2 in range(sh):
                for j2 in range(sh):
                    delta += 1/2 * re_adj2[i1][i2][j1][j2]**2 * (abs(np.real(X[0][i1][j1])) * abs(np.real(X[0][i2][j2])) + abs(np.real(X[1][i1][j1])) * abs(np.real(X[1][i2][j2])))
                    delta += 1/2 * im_adj2[i1][i2][j1][j2]**2 * (abs(np.real(X[0][i1][j1])) * abs(np.real(X[1][i2][j2])) + abs(np.real(X[1][i1][j1])) * abs(np.real(X[0][i2][j2])))
    
    return delta


# Optimize the measurement allocation using constrained gradient descent to minimize the error in the determinant
# to second order
def optimize_Mij_gradient_descent(dims, state, order, combs, eta, N, M_total, 
                                   initial_Mij=None, method='SLSQP', 
                                   verbose=False, max_iterations=1000, tol=1e-30):
    
    sh = len(combs)
    var = full_var(dims, state, order, combs, eta, N)
    Adj1 = adj(dims, state, order, combs, eta, N)
    Adj2 = adj2(dims, state, order, combs, eta, N)
    
    # Create initial guess if not provided, then symmetrize it by construction.
    if initial_Mij is None:
        initial_Mij = np.zeros((2, sh, sh), dtype=float)
        total_unique = sh * (sh + 1) // 2 + sh * (sh - 1) // 2
        uniform_val = M_total / total_unique
        initial_Mij[0] = uniform_val
        initial_Mij[1] = uniform_val
        initial_Mij[1][np.diag_indices(sh)] = 0

        if 0 in combs:
            total_unique -= 1
            uniform_val = M_total / total_unique
            initial_Mij[0] = uniform_val
            initial_Mij[1] = uniform_val
            initial_Mij[1][np.diag_indices(sh)] = 0
            initial_Mij[0][0][0] = 0
    else:
        initial_Mij = np.array(initial_Mij, dtype=float, copy=True)
        if initial_Mij.shape != (2, sh, sh):
            raise ValueError("initial_Mij must have shape (2, len(combs), len(combs))")

    # Optimize only unique entries: upper triangle of Mij[0] and upper triangle of Mij[1]
    # excluding the diagonal of Mij[1]. The full matrices are reconstructed symmetrically.
    free_indices = []
    for i in range(sh):
        for j in range(i, sh):
            free_indices.append((0, i, j))
    for i in range(sh):
        for j in range(i + 1, sh):
            free_indices.append((1, i, j))

    if 0 in combs:
        free_indices = free_indices[1:]

    # Flatten only unique free variables for optimization
    x0 = np.array([initial_Mij[p, i, j] for p, i, j in free_indices], dtype=float)

    def unflatten_Mij_from_free(x):
        """Convert free variables back to a symmetric full (2, sh, sh) array."""
        Mij = np.zeros((2, sh, sh), dtype=float)
        for idx, (p, i, j) in enumerate(free_indices):
            Mij[p, i, j] = x[idx]
            if i != j:
                Mij[p, j, i] = x[idx]
        Mij[1][np.diag_indices(sh)] = 0.0
        return Mij
    
    def objective(x):
        """Objective function: minimize delta_detA2_secondorder."""
        Mij = unflatten_Mij_from_free(x)
        try:
            delta = delta_detA2_secondorder(Mij=Mij, var=var, Adj1=Adj1, Adj2=Adj2, combs=combs, first_order=False)
            result_val = np.real(delta)
            if np.isnan(result_val) or np.isinf(result_val):
                if verbose:
                    print(f"Warning: NaN/inf in objective")
                return 1e10
            return result_val
        except Exception as e:
            if verbose:
                print(f"Objective evaluation failed: {e}")
            return 1e10
    
    # Set bounds: all unique free variables positive
    bounds = [(1e-6, M_total*10) for _ in free_indices]
    
    # Constraint: sum of the full symmetric matrix entries equals M_total.
    def budget_constraint(x):
        """Total measurement budget constraint."""
        return np.sum(unflatten_Mij_from_free(x)) - M_total
    
    constraints = {'type': 'eq', 'fun': budget_constraint}
    
    # Run optimization
    result = optimize.minimize(
        objective, x0, 
        method=method, 
        bounds=bounds,
        tol=tol,
        constraints=constraints,
        options={'verbose': 1 if verbose else 0, 'maxiter': max_iterations}
    )
    
    optimal_Mij = unflatten_Mij_from_free(result.x)
    
    if verbose:
        print(f"\nOptimization result:")
        print(f"  Success: {result.success}")
        print(f"  Final delta: {result.fun:.6e}")
        print(f"  Iterations: {result.nit}")
        print(f"  Message: {result.message}")
        print(f"  Total budget used: {np.sum(optimal_Mij):.2f}")
        print(f"  Mij[0] sum: {np.sum(optimal_Mij[0]):.2f}")
        print(f"  Mij[1] off-diagonal sum: {np.sum(optimal_Mij[1]) - np.sum(optimal_Mij[1, np.diag_indices(sh)]):.2f}")
    
    return result, optimal_Mij


# The optimal measurement allocation is generally non-integer.
# This function finds a nearby integer distribution of measurements that 
# preserves the total number of measurements Mtot and approximates the best integer solution.
def integerisation(M, M_ij, order, combs, print_m=False, identityincluded=False):
    # smu = gamma / M

    # check that M is integer (or a float with zero fractional part)
    if not M.is_integer():
        raise ValueError("M must be an integer or a float with zero fractional part.")

    n = M_ij[0].shape[0]

    # First set unrounded M_ijs for the real and imaginary parts
    raw_re = np.zeros((n, n))
    raw_im = np.zeros((n, n))

    for i in range(n):
        raw_re[i][i] = np.real(M_ij[i][i])

    for i in range(n):
        for j in range(i + 1, n):
            raw_re[i][j] = np.real(M_ij[i][j]) 
            raw_im[i][j] = np.imag(M_ij[i][j]) 
            raw_re[j][i] = raw_re[i][j]
            raw_im[j][i] = raw_im[i][j]

    # Round all measurement values to nearest integer
    M_re = np.round(raw_re)
    M_im = np.round(raw_im)

    # enforce minimum of 1 for entries for which rounding gives 0
    # do not enforce for the imaginary part when moment is positive definite
    posdef = positivedefinite_count(order, combs)

    # M_re = np.where(M_re < 1, 1, M_re)
    # M_im = np.where(M_im < 1, 1, M_im)
    M_re = np.maximum(M_re, 1)
    M_im = np.maximum(M_im, 1)
    M_im[posdef] = 0
    if identityincluded:
        M_re[0][0] = 0

    if print_m:
        print('Raw real part:')
        print(raw_re)
        print('Raw imaginary part:')
        print(raw_im)
        print('Rounded + subbed real part:')
        print(M_re)
        print('Rounded imaginary part:')
        print(M_im)

    # Weight array for excess / deficit reallocation
    # Diagonal weights are real only, off-diagonal contribute both re and im
    weight_re = np.zeros((n, n))
    weight_im = np.zeros((n, n))

    for i in range(n):
        weight_re[i][i] = np.real(M_ij[i][i])

    for i in range(n):
        for j in range(i + 1, n):
            weight_re[i][j] = weight_re[j][i] = np.real(M_ij[i][j])
            weight_im[i][j] = weight_im[j][i] = np.imag(M_ij[i][j])

    # Redistribution guided by the weights
    current_sum = int(np.sum(M_re) + np.sum(M_im))
    residual = M - current_sum

    # Stack real and imaginary weights into a single pool
    # First n*n entries = real part, next n*n = imaginary part
    flat_weight_re = weight_re.flatten()
    flat_weight_im = weight_im.flatten()
    all_weights = np.concatenate([flat_weight_re, flat_weight_im])
    n2 = n * n

    if residual > 0:
        # Add 1 to the highest-weight entries first
        indices = np.argsort(-all_weights)  # descending

        # for k in range(int(residual)):
        #     idx = indices[k]
        #     if idx < n2:
        #         M_re[np.unravel_index(idx, (n, n))] += 1
        #     else:
        #         M_im[np.unravel_index(idx - n2, (n, n))] += 1
        
        removed = 0 
        counter = 0
        for idx in cycle(indices):
            counter += 1
            if removed == residual:
                break
            if counter > abs(residual) * 2 * n2**2:  # safety break to prevent infinite loop
                print("Warning: Integerisation failed (residual > 0). Double check input parameters.")
                break
            if idx < n2:
                M_re[np.unravel_index(idx, (n, n))] += 1
            elif idx >= n2 and M_im[np.unravel_index(idx - n2, (n, n))] < 1:
                continue
            else:
                M_im[np.unravel_index(idx - n2, (n, n))] += 1
            removed += 1
        

    elif residual < 0:
        # Remove 1 from the lowest-weight entries first, before minimum is enforced
        # Ensure we don't remove from entries that are equal to 1
        indices = np.argsort(all_weights)   # ascending
        removed = 0
        counter = 0
        for idx in cycle(indices):
            counter += 1
            if removed == abs(residual):
                break
            if counter > abs(residual) * 2 * n2**2:  # safety break to prevent infinite loop
                print("Warning: Could not remove enough entries without violating minimum constraints.")
                break
            # if entry is already at minimum, skip it
            # if idx < n2 and M_re[np.unravel_index(idx, (n, n))] <= 1:
            #     continue
            # if idx >= n2 and M_im[np.unravel_index(idx - n2, (n, n))] <= 1:
            #     continue
            if idx < n2 and M_re[np.unravel_index(idx, (n, n))] <= 1:
                continue
            if idx >= n2 and M_im[np.unravel_index(idx - n2, (n, n))] <= 1:
                continue
            if idx < n2:
                M_re[np.unravel_index(idx, (n, n))] -= 1
            else:
                M_im[np.unravel_index(idx - n2, (n, n))] -= 1
            removed += 1

    return np.array([M_re, M_im])


# For a given M_ij, list of matrix indices, and state, this function creates a sample 
# of matrices by sampling from each of the matrix element distributions.
# The function checks the kurtosis of the distributions to ensure that the number 
# of measurements is sufficient for the Gaussian approximation to be valid. It stores 
# a warning/flag if the excess kurtosis is above a specified threshold.
def sample_matrix(dims, state, order, combs, eta, N, n_samps, M, M_ij, print_m=False, equal_distrib_test=False, kurtosis_threshold=4, return_kurtosis=False):
    sh = len(combs)
    Adj = adj(dims, state, order, combs, eta, N)
    re_adj = np.real(Adj)
    im_adj = np.imag(Adj)
    X_r, X_i = full_var(dims, state, order, combs, eta, N)
    true_mat = TD_mat(dims, state, order, combs, eta, N)
    rng = np.random.default_rng()
    out = np.zeros((n_samps, sh, sh), dtype=np.complex128)
    # gamma = Gamma(dims, state, order, combs, eta, N)
    # smu = gamma/M
    
    # determine distribution of measurement numbers
    identityincluded = False
    if 0 in combs:
        identityincluded = True

    minima_matrix = np.ceil(abs(full_kurtosis(dims, state, order, combs, eta, N) / kurtosis_threshold))
    # to avoid kurtosis of diagonal entries blowing up for low squeezing
    # set minimum of 1, calculate measurements from weight matrix, sample from gamma function

    posdef = positivedefinite_count(order, combs)
    notposdef = ~posdef

    minima_matrix[0][posdef] = 1
    minima_matrix[minima_matrix < 1] = 1
    minima_matrix[1][posdef] = 0

    minima_matrix[np.isnan(minima_matrix)] = 0

    if print_m:
        print('Minima matrix:')
        print(minima_matrix)

    kurtosistest = minima_matrix / M_ij
    kurtosistest[np.isnan(kurtosistest)]=0

    if identityincluded:
        # safe guard, shape-dependent: only set when shapes match
        try:
            kurtosistest[0][0][0] = 0
        except Exception:
            pass

    kurt_flag = False
    if np.any(kurtosistest > kurtosis_threshold):
        kurt_flag = True
        print(str(combs), "Optimal distribution does not satisfy CLT")
        print('Overshot elements: ',  np.argwhere(kurtosistest>kurtosis_threshold))
        print('Values: ', kurtosistest[kurtosistest > kurtosis_threshold])

    if sum(minima_matrix.flatten()) > M:
        kurt_flag = True
        #raise ValueError("Minimum measurement numbers required for CLT exceed total M_tot. Consider reducing kurtosis threshold or increasing M_tot.")
        print(str(combs), "Minimum measurement numbers required for CLT exceed total M_tot. Consider reducing kurtosis threshold or increasing M_tot.")


    if equal_distrib_test:
        equaldistrib = np.ceil(np.sum(np.real(M_ij) + np.imag(M_ij))/(2 * sh**2 - sh))
        M_ij = np.full_like(M_ij, equaldistrib + 1j*equaldistrib)
        M_ij[np.diag_indices(sh)] = equaldistrib


    M_ij = M_ij[0] + 1j*M_ij[1]

    if print_m:
        print('M_ij:')
        print(M_ij)

    posdef = positivedefinite_count(order, combs)
    notposdef = ~posdef

    mval = np.real(M_ij[posdef])
    out[:,posdef] = rng.gamma(shape = mval * np.real(true_mat[posdef])**2 / abs(np.real(X_r[posdef])), scale=abs(np.real(X_r[posdef]))/(mval * np.real(true_mat[posdef])), size=(n_samps, len(mval)))

    # off diagonal elements which are not Hermitian,
    # these stand for gaussian upper and gaussian lower
    gupper = np.triu(notposdef, k=1)
    glower = np.tril(notposdef, k=-1)

    mval_r = np.real(M_ij[gupper])
    mval_i = np.imag(M_ij[gupper])
    out[:,gupper] = rng.normal(np.real(true_mat[gupper]),np.sqrt(abs(np.real(X_r[gupper])/mval_r)),size=(n_samps, len(mval_r))) + 1j*rng.normal(np.imag(true_mat[gupper]),np.sqrt(abs(np.real(X_i[gupper])/mval_i)),size=(n_samps, len(mval_i)))
    out[:,glower] = rng.normal(np.real(true_mat[glower]),np.sqrt(abs(np.real(X_r[gupper])/mval_r)),size=(n_samps, len(mval_r))) + 1j*rng.normal(np.imag(true_mat[glower]), np.sqrt(abs(np.real(X_i[gupper])/mval_i)),size=(n_samps, len(mval_i)))

    if identityincluded:
        out[:,0,0] = 1.0
    if return_kurtosis:
        return out, kurt_flag
    return out


# Calculate the determinant of the sampled matrices
def sample_determinant(dims, state, order, combs, eta, N, M_ij, n_samps=10000, M=1000, equal_distrib_test=False, return_kurtosis=False):
    
    result = sample_matrix(dims, state, order, combs, eta, N, n_samps, M, M_ij, equal_distrib_test=equal_distrib_test, return_kurtosis=return_kurtosis)
    if return_kurtosis:
        matrix_samples, kurt_flag = result
        dets = np.linalg.det(matrix_samples)
        return np.real(dets), kurt_flag
    else:
        matrix_samples = result
        dets = np.linalg.det(matrix_samples)
        return np.real(dets)
    


def statistical_power_mc(dims, state, order, combs, eta, N, M, M_ij, mc_samps=10000, equal_distrib_test=False, return_kurtosis=False):
    # monte carlo estimate of the distribution of determinant values, from which we can calculate the statistical power
    if return_kurtosis:
        det_samples, kurt_flag = sample_determinant(dims, state, order, combs, eta, N, M_ij, n_samps=mc_samps, M=M, equal_distrib_test=equal_distrib_test, return_kurtosis=True)
    else:
        det_samples = sample_determinant(dims, state, order, combs, eta, N, M_ij, n_samps=mc_samps, M=M, equal_distrib_test=equal_distrib_test)
    # return det samples and number of samples below zero
    percentiles = np.array([np.percentile(det_samples, 68.27 + 0.5 * (100 - 68.27)), np.percentile(det_samples, 97.5), np.percentile(det_samples, (100-68.27) * 0.5), np.percentile(det_samples, 2.5)])
    if return_kurtosis:
        return det_samples, len(det_samples[det_samples<0])/len(det_samples), percentiles, kurt_flag
    return det_samples, len(det_samples[det_samples<0])/len(det_samples), percentiles


def full_kurtosis(dims, state, order, combs, eta, N):

    sh = len(combs)
    mat_str = sub_matrix(order, combs)

    K_real = np.full((sh, sh), np.nan, dtype=float)
    K_im = np.full((sh, sh), np.nan, dtype=float)

    def to_4index_terms(n, m, k, l, p, q, r, s):
        terms = []
        for f in range(min(m, k) + 1):
            for g in range(min(q, r) + 1):
                coeff = G(m, k, f) * G(q, r, g)
                terms.append([coeff, (
                n + k - f,
                m - f + l,
                p + r - g,
                q - g + s
                )])
        return terms

    def multiply_strings(terms, N2, M2, K2, L2):
        new_terms = []
        for coeff, (N1, M1, K1, L1) in terms:
            for f in range(min(M1, N2) + 1):
                for g in range(min(L1, K2) + 1):
                    nc = coeff * G(M1, N2, f) * G(L1, K2, g)
                    new_terms.append([nc, (
                    N1 + N2 - f,
                    M1 - f + M2,
                    K1 + K2 - g,
                    L1 - g + L2
                    )])
        return new_terms

    def expect_4index_terms(terms):
        result = 0.
        for coeff, (nn, mm, kk, ll) in terms:
            result += coeff * H(dims, state, nn, mm, kk, ll, eta, N)
        return result

    def multiply_by_op(current_terms, op_terms):
        new_terms = []
        for c2, (N2, M2, K2, L2) in op_terms:
            partial = multiply_strings(current_terms, N2, M2, K2, L2)
            for cp, idx in partial:
                new_terms.append([c2 * cp, idx])
        return new_terms

    def compute_moment(O_terms, Od_terms, power, scale, s):
        if power == 0:
            return 1.0

        total_terms = []

        for t in range(power + 1):
            for o_positions in combinations(range(power), t):
                o_pos_set = set(o_positions)
                current = [[1.0, (0, 0, 0, 0)]]
                for pos in range(power):
                    if pos in o_pos_set:
                        current = multiply_by_op(current, O_terms)
                    else:
                        current = multiply_by_op(current, Od_terms)
                sign = s ** (power - t)
                for c, idx in current:
                    total_terms.append([sign * c, idx])

        result = expect_4index_terms(total_terms)
        return result / (scale ** power)

    def excess_kurtosis(mu1, mu2, mu3, mu4):
        var = mu2 - mu1**2
        if np.abs(var) < 1e-8:
            return np.nan
        c4 = mu4 - 4*mu3*mu1 + 6*mu2*mu1**2 - 3*mu1**4
        return float(np.real(c4 / var**2)) - 3.0
    
    # Main loop over upper triangle only
    for i in range(sh):
        for j in range(i, sh):
            n, m, k, l, p, q, r, s = tuple(string_to_ind(mat_str[i][j]))
            nd, md, kd, ld, pd, qd, rd, sd = tuple(string_to_ind(mat_str[j][i]))

            O_terms = to_4index_terms(n, m, k, l, p, q, r, s )
            Od_terms = to_4index_terms(nd, md, kd, ld, pd, qd, rd, sd)

            if i == j:
            # Diagonal: O is Hermitian, compute <O^r> directly
                mu1 = np.real(compute_moment(O_terms, Od_terms, 1, 1, +1))
                mu2 = np.real(compute_moment(O_terms, Od_terms, 2, 1, +1))
                mu3 = np.real(compute_moment(O_terms, Od_terms, 3, 1, +1))
                mu4 = np.real(compute_moment(O_terms, Od_terms, 4, 1, +1))

                K_real[i][i] = excess_kurtosis(mu1, mu2, mu3, mu4)
                # Diagonal imaginary part is zero (O is Hermitian, no imaginary part)
                K_im[i][i] = 0.0

            else:
                # Off-diagonal: A = (O + O†)/2
                mu1_A = np.real(compute_moment(O_terms, Od_terms, 1, 2, +1))
                mu2_A = np.real(compute_moment(O_terms, Od_terms, 2, 2, +1))
                mu3_A = np.real(compute_moment(O_terms, Od_terms, 3, 2, +1))
                mu4_A = np.real(compute_moment(O_terms, Od_terms, 4, 2, +1))

                K_real[i][j] = excess_kurtosis(mu1_A, mu2_A, mu3_A, mu4_A)
                # Lower triangle mirrors upper triangle
                K_real[j][i] = K_real[i][j]

                # C = (O - O†)/(2i)
                mu1_C = np.real(compute_moment(O_terms, Od_terms, 1, 2j, -1))
                mu2_C = np.real(compute_moment(O_terms, Od_terms, 2, 2j, -1))
                mu3_C = np.real(compute_moment(O_terms, Od_terms, 3, 2j, -1))
                mu4_C = np.real(compute_moment(O_terms, Od_terms, 4, 2j, -1))

                K_im[i][j] = excess_kurtosis(mu1_C, mu2_C, mu3_C, mu4_C)
                # Lower triangle mirrors upper triangle
                K_im[j][i] = K_im[i][j]

    return np.array([K_real, K_im])


# Function to assist plotting
# Solid lines when the kurtosis is below the threshold,
# dotted lines when it is above

def plot_clt(kurt_exceeded, statistical_power, xvals, linestyle, color, dottedlinewidth=6, solidlinewidth=3, splitmarkersize=15):
        import matplotlib.pyplot as plt
        transition = np.where(kurt_exceeded)[0]
        splits = np.where(np.diff(transition) != 1)[0] + 1
        runs = np.split(transition, splits)

        if len(transition) == 0:
            plt.plot(xvals, statistical_power, linestyle=linestyle, color=color, lw=solidlinewidth)
            return

        start = 0
        for run in runs:
            # solid
            if start < run[0]:
                plt.plot(xvals[start:run[0]], statistical_power[start:run[0]], linestyle=linestyle, color=color, lw=solidlinewidth)

            # dotted
            i0 = max(run[0]-1, 0)
            i1 = min(run[-1]+1, len(xvals)-1)

            plt.plot(xvals[i0:i1+1], statistical_power[i0:i1+1], linestyle=(0, (1, 3)), color=color, lw=dottedlinewidth)

            if i0 != 0:
                plt.plot(xvals[i0], statistical_power[i0], '.', color=color, markersize=splitmarkersize)
            if i1 != len(xvals)-1:
                plt.plot(xvals[i1], statistical_power[i1], '.', color=color, markersize=splitmarkersize)

            start = run[-1]

        if start < len(xvals)-1:
            plt.plot(xvals[start:], statistical_power[start:], linestyle=linestyle, color=color, lw=solidlinewidth)

        return


# The functions below can be useful to confirm that the measurement allocation is optimal
# (used in development)

def sjt_permutations(n: int) -> Generator[Tuple[int], None, None]:
    """An implementation of the Steinhaus–Johnson–Trotter
    algorithm with Even's speedup for generating
    permutations in order of alternating parity."""

    # Each element of the permutation is
    # initialized with a direction.
    # All elements face left initially.
    yield tuple(p := list(range(n))), (parity := 1)
    d = [0] + [-1] * (n - 1)

    # While any element is still directional,
    while any(d):
        # get the greatest element with nonzero direction,
        # its index, and the index of the element it faces.
        c, i, j = max((p[i], i, i + j) for i, j in enumerate(d) if j)

        # All undirected elements greater than c have
        # their directions set to face toward c.
        d = [d[k] or (p[k] > c) * ((i > k) * 2 - 1) for k in range(n)]

        # Swap the directions and elements at indices i and j.
        # If the swap puts element c at the start or end
        # of the permutation, or if the next element in the same direction
        # is greater than c, zero out c's direction instead.
        d[i], d[j] = d[j], 0 if j in {0, n - 1} or p[j + d[i]] > c else d[i]
        p[i], p[j] = p[j], p[i]

        # Permutation parity alternates.
        yield tuple(p), (parity := -parity)


def levi_civita(n: int) -> sparse.COO:
    """Generates the sparse Levi-Civita symbol in n-dimensions."""
    p, v = zip(*sjt_permutations(n))
    return sparse.COO(tuple(zip(*p)), v, shape=(n,) * n)

def sanity_check(dims, state, order, ls, eta, N, optimal_Mij_ls):
    sh = len(ls)
    mat = TD_mat(dims, state, order, ls, eta, N)
    levi_civita0 = np.array(levi_civita(sh).todense())
    adjoint = adj(dims,state,order,ls,eta,N)
    lambdas0 = np.zeros((sh,sh), dtype=np.complex128)
    lambdas1 = np.zeros((sh,sh), dtype=np.complex128)
    mat_tensor = 1
    for i in range(sh-2):
        mat_tensor = np.tensordot(mat_tensor, mat, 0)
    indices = list(string.ascii_lowercase)[:2*sh]
    levi_civita_tensor = np.tensordot(levi_civita0, levi_civita0, 0)
    var = full_var(dims, state, order, ls, eta, N)
    weightvar = weighted_var(optimal_Mij_ls,var,ls)
    all_indices = ''.join(indices)
    cont_indices = ''.join(indices[2:-2])
    cont_indices_order = ''
    for i in range(len(cont_indices)//2):
        cont_indices_order +=cont_indices[i]
        cont_indices_order +=cont_indices[i+len(cont_indices)//2]
    # print(all_indices+','+cont_indices_order+'->'+all_indices[:2]+all_indices[-2:])
    # print(np.shape(levi_civita_tensor))
    # print(np.shape(mat_tensor))
    second_der_real=(np.real(np.einsum(all_indices+','+cont_indices_order+'->'+all_indices[:2]+all_indices[-2:],levi_civita_tensor,mat_tensor)))**2
    second_der_imag=(np.imag(np.einsum(all_indices+','+cont_indices_order+'->'+all_indices[:2]+all_indices[-2:],levi_civita_tensor,mat_tensor)))**2
    second_order_realq0 = np.einsum('ikjl,kl->ij',second_der_real,weightvar[0])
    second_order_realq1 = np.einsum('ikjl,kl->ij',second_der_imag,weightvar[1])
    second_order_imagq0 = np.einsum('ikjl,kl->ij',second_der_imag,weightvar[0])
    second_order_imagq1 = np.einsum('ikjl,kl->ij',second_der_real,weightvar[1])
    for i in range(sh):
        lambdas0[i][i] += 1/math.factorial(sh-2)**2*second_order_realq0[i][i]*weightvar[0][i][i]/optimal_Mij_ls[0][i][i]
        lambdas0[i][i] += 1/math.factorial(sh-2)**2*second_order_realq1[i][i]*weightvar[0][i][i]/optimal_Mij_ls[0][i][i]
        lambdas0[i][i] += np.real(adjoint[i][i])**2*weightvar[0][i][i]/optimal_Mij_ls[0][i][i] 

    for i in range(sh):
        for j in range(i+1,sh):
            lambdas0[i][j] += 1/math.factorial(sh-2)**2*second_order_realq0[i][j]*weightvar[0][i][j]/optimal_Mij_ls[0][i][j]
            lambdas0[i][j] += 1/math.factorial(sh-2)**2*second_order_realq1[i][j]*weightvar[0][i][j]/optimal_Mij_ls[0][i][j]
            lambdas0[i][j] += np.real(adjoint[i][j])**2*weightvar[0][i][j]/optimal_Mij_ls[0][i][j]   
            lambdas1[i][j] += 1/math.factorial(sh-2)**2*second_order_imagq0[i][j]*weightvar[1][i][j]/optimal_Mij_ls[1][i][j]
            lambdas1[i][j] += 1/math.factorial(sh-2)**2*second_order_imagq1[i][j]*weightvar[1][i][j]/optimal_Mij_ls[1][i][j]
            lambdas1[i][j] += np.imag(adjoint[i][j])**2*weightvar[1][i][j]/optimal_Mij_ls[1][i][j]

            lambdas0[j][i] = lambdas0[i][j]
            lambdas1[j][i] = lambdas1[i][j]
    return np.array([lambdas0, lambdas1])

