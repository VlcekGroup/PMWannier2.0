# PMWannier2.0
The second version (V2.0) of Pipek Mezey wannier codes with local, fragment, and sequential wannierization routines for large systems where molecules are place in an environment.

# Manual for PMwannier package
This code performs localization on canonical orbitals employing the Pipek-Mezey localization scheme. The objective functional is maximized using the steepest ascent algorithm. It is primarily designed for the localization of occupied states but also available for virtual states if a specific virtual space is defined. The main feature of this code is to efficiently obtain a set of regionally localized states on the atoms (molecules) of interest in a large system (usually with over 4000 electrons). The PM functional is modified to evaluate the selected atoms only. The orbital space can be either transformed all at once or sequentially exhausted even if there is no knowledge of the single-particle states a priori.\
\
The theoretical fundaments of local, fragment, and sequential wannierization can be found in (1) G. Weng, V. Vlcek  J. Chem. Phys. 155, 054104 (2021); https://doi.org/10.1063/5.0058410, and (2) G. Weng, V. Vlcek, et al., arXiv:2204.00252; https://doi.org/10.48550/arXiv.2204.00252. This program is mostly written by the authors, but it has two publicly available subroutines, primarily including the calculation of the exponential of a matrix provided by Cleve Moler and Charles Van Loan (SIAM Review, Volume 45, Number 1, March 2003, pages 3-49.), and the Gamma function provided by John Burkardt (Handbook of Mathematical Functions, National Bureau of Standards, 1964 and The Mathematica Book, Fourth Edition, Cambridge University Press, 1999)

The PMWannier1.0 work was supported by the NSF CAREER award through Grant No. DMR-1945098. The PMWannier2.0 covers all the functions in the old version. And this work was supported by the Office of Science of the U.S. Department of Energy under Contract No. DE-AC02-05CH11231 using NERSC award BES-ERCAP0020089.

## Installation
In the source code folder "src", simply hit "make" to install the program wannier.x, where the openmpi library is needed.

## Wannier.x
### Command to run a job
~/path_to_executable_file/wannier.x input_binfile_name >& output_file_name &

### Files needed to run a job
(1) cnt.ini : a text file that provides the coordinates of atoms.\
(2) wavefunction bin file : a binary file that provides canonical orbitals for start-over or transformed orbitals for restart. This file can have different names.\
(3) input : a text file that specifies the parameters for the calculation, see Input Specifications. The name of this file must be "input".\
(4) group_atom : a text file that contains a description line, the number of regionally localized states, the number of atoms of interest, and the atom list.
### Output files after running a job
(1) Wannier.bin : a binary file that contains all the unitary-transformed states in the occupied/unoccupied states.\
(2) WannierXXX.bin : a binary file that contains all the unitary-transformed states in the occupied/unoccupied states for the fragment XXX if fragment wannierization is used.\
(3) wannier_trun.bin : a truncated version of the wannier.bin, which contains the regionally localized states or the core space localized states only.\
(4) wannierXXX_trun.bin : a truncated version of the wannierXXX.bin, which contains the regionally localized states or the core space localized states one the fragment XXX only.\
(5) U_tmp.txt : a text file that stores the temporary unitary matrix at the current inner-loop step for restart.\
(6) stepXXX.bin : a binary file that stores the transformed states at the outer-loop step XXX if deterministic sequential wannierization is used. If the iteration exceeds 1000 steps, the 1001st-step file will rewrite the 1st-step file.\
(6) mcstepXXX.bin : a binary file that stores the transformed states at the outer-loop step XXX if stochastic sequential wannierization is used. If the iteration exceeds 1000 steps, the 1001st-step file will rewrite the 1st-step file.\
(7) OUTPUT (or other name) : by default the output information is written to screen. File with the name "OUTPUT" contains some of the input information, the objective function value at each iteration, the time spent at each step, and the delta_t value at each step, in both the inner- and the outer-loop.

### Input Specification
#### Basic Parameters
restart\
(Logical. This controls whether your job starts from canonical orbitals (f) or restarts from transformed orbitals (t). Default: f)\

max_iter\
(Integer. This specifies the maximum iteration steps in the steepest ascent process (inner-loop). Default: 500)\

delta_t\
(Real. This specifies the initial step of the step ascent process. Default: 5.0)\

delta_t_correction\
(Real. This specifies the correction factor for delta_t if the change of the objective function is negative. The delta_t is corrected by taking the new delta_t = delta_t / delta_t_correction. Default: 1.1)\

num_of_threads\
(Integer. This specifies the number of threads for Openmp processes. Default: 1)\

convergence_threshold_inner\
(Real. This specifies the convergence threshold for the steepest ascent process (the inner-loop). Default: 1E-7)\

nhit_convergence_inner\
(Integer. This specifies the number of times of consecutively hitting the convergence criterion in the inner-loop before exiting. Requirement of multiple hittings guarantees a smooth convergence behavior. Default: 3)\

atom_space\
(Character. This specifies whether the localization is performed on all atoms, selected atoms, or a giant "atom". There are three options: full, local, and fragment, where full means all atoms, local means selected atoms, and fragment means giant "atom". Default: full)\

step_prnt_info\
(Integer. This specifies the interval to print the output information during the steepest ascent process. Default: 50)\

space\
(Character. This specifies in which space the localization if performed. There are three options: occ, unocc, and random. Randoms means one can define the space by just providing the orbitals, either occupied or unoccupied. Default: occ\
If you choose "random", in the next three lines provid the following information:\
1st line: Integer. This specifies the number of states that define your random space. No default value\
2nd line: Logical. This specifies whether your state indices are continuous in your bin file. No default value\
3rd line: Integers. This specifies the state indices to read from your bin file. If continuous, provide the 1st and the last indices. Otherwise, provide the full list of indices. )


### All-atom localization
#### If you specify the atom_space as "full", no further parameter is needed. You are done here.

### Selected-atom localization
#### If you specify the atom_space as "local", specify the following parameters.
local_wannierization\
(1st line: Integer. This specifies the number of regionally localized states on the selected atoms. Default: -1)\
(2nd line: Integer. This specifies the number of selected atoms. Default: 1)\
(3rd line: Integer. This specifies the atom list in the cnt.ini file. No default value)\

#### If you specify the atom_space as "fragment", specify the following parameter.
fragment_wannierization\
(Integer. This specifies the number of fragments for the localization. Note that sequential localization applies to one fragment only in the current version. Default: 1)

#### You can choose different ways to transform the orbitals on the selected atoms by specifying "exhaust_space".
exhaust_space\
(Character. This specifies how the orbital space is transformed. There are three options: simple, sequential, and stochastic. Default: simple)

#### If you choose "simple" in "exhaust_space", specify the following parameters.
space_size\
(Character. This specifies the size of the occupied or the unoccupied space. There are two options: full, and sub. Full occupied space takes all the occupied states while full unoccupied states takes all the unoccupied states in the bin file. Default: full)

#### If you choose "sub" in "space_size", specify the following parameters.
define_subspace\
(1st line: Character. This specifies how to determine the size of the subspace. There are three options: loc_cut, nstates, energy. Default: loc_cut)\
(2nd line: Depends on the 1st line. If "loc_cut" is choose, provide a "Real" number (Default: 1E-3) to specify the locality cutoff when selecting the states. If "nstates" is chosen, provide an "Integer" number (No default value) to specify the number of states to take based on the locality order. If "energy"  is chosen, provide a "Real" number (No default value) to specify the energy window for the unoccpied states only.)

#### If you choose "sequential" or "stochastic" in "exhaust_space", specify the following parameters.
convergence_threshold_outer\
(Real. This specifies the convergence criterion for the outer-loop. Default: 5E-7)\

step_prnt_tmp\
(Integer. This specifies the interval to output the temporary transformed bin file for restart. Default: 100)\

#### If you choose "sequential" in "exhaust_space", specify the following parameters.
nstates_seqwan\
(Integer. This specifies the number of states in the full space that is being sequentially exhausted. The code reads these states from the 1st states to the specified states in the bin file. Default: -1)\

nstates_core\
(Integer. This specifies the size of the core space. It should be greater than or equal to the number of regionally localized states. Default: -1)\

nstates_r\
(Integer. This specifies the step of exhausting the rest space. Default: 1)\

istep_0\
(Integer: This specifies the initial step in the outer-loop. For restart, specify it from where the previous job ends. Default: 1)\

total_step\
5000\
(Integer. This specifies the maximum iteration steps in the outer-loop. Default: 1000)\

reorder_rest\
(Logical. This specifies whether to do reordering of the rest space or not when one complete accessing is done. Default: t)\

#### If you choose "sequential" in "exhaust_space", specify the following parameters.
nstates_stowan\
(Integer. This specifies the number of states in the full space that is being sequentially exhausted. The code reads these states from the 1st states to the specified states in the bin file. Default: -1)\

nstates_det\
(Integer. This specifies the size of the deterministic space. It should be greater than or equal to the number of regionally localized states. Default: -1)\

nstates_sto\
(Integer. This specifies the number of stochastic states that sample the rest space. Default: 1)\

i_mc_step0\
1\
(Integer: This specifies the initial step in the outer-loop. For restart, specify it from where the previous job ends. Default: 1)\

mc_step\
(Integer. This specifies the maximum iteration steps in the outer-loop. Default: 1000)\

nhit_cong_outer\
(Integer. This specifies the number of times of consecutively hitting the convergence criterion in the outer-loop before exiting. Default: 5)\

## Typical Input Examples for wannier.x
### One benzene molecule
restart\
f\
max_iter\
500\
delta_t\
5.0\
delta_t_correction\
1.1\
num_of_threads\
10\
convergence_threshold_inner\
1E-7\
nhit_convergence_inner\
3\
atom_space\
full\
step_prnt_info\
50\
space\
occ

### Nine benezene molecules
restart\
f\
max_iter\
2000\
delta_t\
5.0\
delta_t_correction\
1.1\
num_of_threads\
10\
convergence_threshold_inner\
1E-7\
nhit_convergence_inner\
3\
atom_space\
local\
local_wannierization\
15\
12\
48 35 9 22 61 73 27 40 1 14 54 66\
fragment_wannierization\
1\
step_findwf\
100\
step_prnt_info\
100\
exhaust_space\
sequential\
space\
occ\
space_size\
full\
define_subspace\
loc_cut\
1E-3\
convergence_threshold_outer\
5E-7\
step_prnt_tmp\
1000\
nstates_seqwan\
135\
nstates_core\
15\
nstates_r\
30\
istep_0\
1\
total_step\
5000\
reorder_rest\
t\
nstates_stowan\
-1\
nstates_det\
-1\
nstates_sto\
1\
i_mc_step0\
1\
mc_step\
10000\
nhit_cong_outer\
5


