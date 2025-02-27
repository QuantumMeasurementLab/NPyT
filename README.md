# NPyT: The (N)PT test optimization (Py)thon (T)oolbox 

NPyT is a toolbox in Python for selecting the optimal NPT criterion by determining the level of confidence of criteria within the Shchukin and Vogel hierarchy [[Phys. Rev. Lett. 95, 230502 (2005)]](https://doi.org/10.1103/PhysRevLett.95.230502).

NPyT is developed by Lydia A. Kanari-Naish building on project discussions with Jack Clarke, Sofia Qvarfort, and Michael R. Vanner.


## Step 1

The function `my_state` performs a search over all submatrices of a given dimension $`d`$ and order $`n`$. This search assumes no coupling to the environment and no sampling errors, as only determinants that are negative in the absence of such environmental imperfections can be negative when such effects are included. In this way, the function `my_state` can used to perform a preliminary search for candidate NPT criteria up to a given dimension and order.

The `my_state` function outputs the values and rows/columns that identify the submatrix from which the determinant is calculated.
From the order and rows/columns of the submatrix, the function `sub_matrix` may be used to reconstruct the matrix in terms of annihilation/creation operators of subsystems A and B.

NPyT uses the following conventions: 
(i) `a`, `b`, `c`, and `d` correspond to the operators $`\hat{a}^\dagger`$, $`\hat{a}`$, $`\hat{b}^\dagger`$, and $`\hat{b}`$, respectively. 
(ii) The order can only be even so order is equal to $`n/2`$.
(iii) Python array indexing, which starts at 0 as opposed to 1 (as in the manuscript).
For example,`sub_matrix(1,[2,4])` outputs the submatrix `array([['ab', 'ac'],['bd', 'cd']], dtype=object)`, which is the submatrix that produces the determinant $`D_\mathrm{I}`$ parameterized by $`d=2`$, $`n=2`$, and rows/columns=(3,5), i.e. 

$$\begin{vmatrix}
  \langle{\hat{a}^\dagger \hat{a}}\rangle & \langle{\hat{a}^\dagger \hat{b}^{\dagger}}\rangle \\
  \langle{\hat{a}\hat{b}}\rangle & \langle{\hat{b}^\dagger \hat{b}}\rangle.
\end{vmatrix}$$


## Step 2

Following this preliminary search, the effects of environmental interactions and sampling errors are calculated on the subset of successful determinants identified from step 1.
The function `TD_det` calculates the determinants in the presence of environmental interactions.
The function `error` calculates the error on each determinant, in the presence of environmental interactions and sampling errors, according to the optimal allocation of measurements derived in the manuscript.

## Step 3

The confidence levels of each determinant may then be calculated. For a given parameter set, the determinant with the highest confidence level is identified as the optimal NPT criterion. 
The confidence levels may be calculated using the function `confidence_level`, which takes two inputs: the determinant and the error on the determinant.

# Installation

Download the *NPyT.py* file to your project folder and import:

```python
from NPyT import *
```

# Examples

Examples of how to calculate optimal NPT criteria using NPyT are given for the TMSV states, the photon subtracted/added TMSV state, and the two-mode Schrodinger cat state in the files *TMSV.py*, *sub_add_TMSV.py*, and *TMSCS.py*.


# Citation

If you use NPyT in your work, please cite the accompanying paper:
Lydia A. Kanari-Naish, Jack Clarke, Sofia Qvarfort, and Michael R. Vanner, "Optimizing confidence in negative-partial-transpose-based entanglement criteria", (link to appear shortly)
