# NPyT: The NPT test optimization Python Toolbox 

NPyT is a toolbox in Python for selecting the optimal NPT criterion by determining the statistical power of tests for criteria within the Shchukin and Vogel hierarchy [[Phys. Rev. Lett. 95, 230502 (2005)]](https://doi.org/10.1103/PhysRevLett.95.230502).

NPyT is developed by Lydia A. Kanari-Naish, Amaya Calvo-Sánchez, and Arjun Gupta, building on project discussions with Jack Clarke, Sofia Qvarfort, and Michael R. Vanner.


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
The function `statistical_power_mc` calculates the distribution of each determinant through a Monte Carlo method, accounting for the presence of environmental interactions and sampling errors, and following the optimal allocation of measurements described in the manuscript.

## Step 3

The statistical power may be calculated using the function `statistical_power_mc`, which takes information about the state and the optimized measurement allocation as inputs, as well as a sample number to build up the determinant distribution. For a given parameter set, the determinant with the highest statistical power is identified as the optimal NPT criterion. 

# Installation

Download the *NPyT.py* file to your project folder and import:

```python
from NPyT import *
```

# Examples

Examples of how to calculate optimal NPT criteria using NPyT are given for the TMSV state, the photon subtracted/added TMSV state, and the two-mode Schrodinger cat state in the files *TMSV.py*, *sub_add_TMSV.py*, and *TMSCS.py*.

For each state, data files are provided containing statistical powers with varying parameters such as total number of measurements and optical efficiency. The code used to generate the files is provided but commented out for speed.


# Citation

If you use NPyT in your work, please cite the accompanying paper:
Lydia A. Kanari-Naish, Amaya Calvo-Sánchez, Jack Clarke, Arjun Gupta, Sofia Qvarfort, and Michael R. Vanner, "Optimizing the statistical power of negative-partial-transpose-based entanglement tests", 
[[arXiv:2502.19624, (2025)]](https://doi.org/10.48550/arXiv.2502.19624).
