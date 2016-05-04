# NEqMD-ToolKit

This toolkit is a set of stand-along command-line programs written for assisting theoretical/computational chemistry research where Non-Equilibrium MD simulation is needed. This toolkit is written primarily to work in the situation where:
1. [GROMACS](http://www.gromacs.org/) is used for MD propagation;
2. [Gaussian 09](http://www.gaussian.com) is used for electronic structure calculation with DFT;
3. [GaussView](http://www.gaussian.com/g_prod/gv5.htm) is used for molecule visualization;
4. [CNDO program](https://nanohub.org/resources/CNDO/) is used for electronic structure calculation with INDO/S semi-empirical method;
5. Block-diagonalization or general Mulliken-Hush (GMH) method is used for computing electronic coupling (H<sub>DA</sub>).

Detailed description can be found below for each tool.

# Usage

This toolkit is developed in Scientific Linux environment and included in each directory of this repo are the source code along with Makefile used to compile the source. [GNU compilers](https://gcc.gnu.org/) are used in the development and for the python scripts, python 2.7.x installation through [ANACONDA from Continuum](https://www.continuum.io/downloads) is recommended and NumPy is required. For most of the tools, the compiled binary excutable has a name started with `g_` to indicate the consistency with GROMACS. Most of them can display help information or usage information by typing `[excutable name] -h` in terminal and one should alway refer to that for detailed usage info. If the `-h` option is not available, then a simple `[excutable name]` call without command line argument could usually display the help information.

## util.lib
This directory contains the source code of the auxiliary library which is used accross the entire repository. All the C code are written by [Zheng Ma](https://github.com/zanemarkson/). The Fortran code which are written for  matrix manipulations, [BLAS][blas], [LAPACK][lapack] and [EXPOKIT][expokit] are truncated and included in the library. Part of the Fortran code is written by [Dr. Peng Zhang](https://scholars.duke.edu/display/per8149892) for facilitating the chemistry-related matrix operations.

[blas]:       http://www.netlib.org/blas/
[lapack]:     http://www.netlib.org/lapack/
[expokit]:    http://www.maths.uq.edu.au/expokit/

## <a id='cndo2bdiag'></a>cndo2bdiag : `g_cndo2bdiag_d`
This tool is used to compute H<sub>DA</sub> using Block-Diagonalization method. It has `-h` option. `-L` takes one argument specifying the name of a .CI file which is used to indicate the donor and acceptor orbital. Degenerate or near-degenerate orbital situation is allowed and needs to be sepecifed in `[degenerate]` entry of .CI file. The main input file, following `-l` tag, is the CNDO output file. The .ndx input file format is extened from the GROMACS indexing .ndx file format (EIF). It is not only compatible with original GROMACS .ndx format, but can also parse atom sequency specified like 
```bash
[ donor ]
3 - 10 19 20 - 25
```

## <a id='cndo2gmh'></a>cndo2gmh: `g_cndo2gmh_d`
This tools is used to compute H<sub>DA</sub> using CI-based GMH method. Usage is very similar to [cndo2bdiag](#cndo2bdiag).

## cndo3gmh: `g_cndo3gmh_d`
Similar to [cndo2gmh](#cndo2gmh), but in this case the GMH calculation is orbital-based, meaning no CI coefficient calculation is needed and corresponding .ORB file should specify the donor and acceptor _ORBITAL_ instead of states.

## dat2hsd: `g_dat2hsd_d`
To convert .dat file which is the input file format for CNDO program into .hsd file format which is used by [DFTB+](http://www.dftb.org).

## dat2inp: `g_dat2inp_d`
To convert .dat file format into Gaussian 09 input (could also with .gjf or .com file extension) file format.

## dft2bdiag: `g_dft2bdiag_d`
To use DFT Hamiltonian matrix obtained from G09 calculation. To obtain orthongonized Hamiltonian, please refer to `rwf2hao.sh` script included in [dft3gmh](#dft3gmh) directory.

## <a id='dft3gmh'></a>dft3gmh: `g_dft3gmh_mo.dev.py`
To perform GMH calculation using DFT Hamiltonian obtained from G09 calculations. This script is written in python and severy other auxiliary python/bash scripts are included here. In order to use these scripts, G09 read-write-files (RWFs) need to be explicitly saved. These scripts are name in a semantic way so that one can easily understand the usage of these scripts.

## fakeG09mo: `g_fakeG09mo_d`
The `g_*2bdiag_d` tools when `-g yes` option is specified, will generate block-diagonalized MO. By using `g_fakeG09mo_d`, one can feed this type of MO coefficients into GaussView to visualize the orbitals.

## fakefreq: `g_fakefreq_d`
One can use [`g_gmxfreq_d`](#gmxfreq) to calculate the vibrational frequencies and normal modes. By using this tool, one can feed self-calculated frequencies and modes into GaussView to visualize.

## fchkhess2gmx: `g_fchkhess2gmx_d`
To convert Hessian matrix printed in Guassian .fchk file into GROMACS Hessian file format. Note that the unit conversion is done too.

## g09hess2gmx: `g_g09hess2gmx_d`
To convert Hessian matrix printed in Guassian .rwf (RWF number 584) file into GROMACS Hessian file format. Note that the unit conversion is done too.

## <a id='gmxfreq'></a>gmxfreq: `g_gmxfreq_d`
To calculate vibrational modes and associated frequency based on GROMACS-format Hession matrix.

## gro2dat: `g_gro2dat_d`
To convert GROMACS .gro format file into CNDO .dat format file.

## gro2inp: `g_gro2inp_d`
To convert GROMACS .gro format file into Gaussian 09 input file format.

## inp2dat: `g_inp2dat_d`
To convert G09 input file format into CNDO .dat format file.

## modestretch: `g_neqicgen_wigner`
This tool is at the center of this toolkit since it is responsible for :
- Sampling initial geometry and velocity of molecule;
- and stretching molecular geometry along chosen normal modes.

For the initial sampling, Wigner-like probability distribution function is used to determine the position and velocity while the vibrational levels of non-excited modes are determined based on quantized Boltzmann distribution. Detailed usage info can be displayed by calling `g_neqicgen_wigner -h`.

## neqalign: `g_neqalign_vdwradii`
This tool could be used to substitube the solute in one configuration with another solute geometry and remove all the solvent molecules that have van der Waal clash with new solute molecule. The new solute molecule will be rotated and "_aligned_" according to the original solute molecule orientation.

## nmassgn: `g_nmassgn_d`
This tool takes the output files from [`g_gmxfreq_d`](#gmxfreq) and a [EIF](#cndo2bdiag) file to determine the localization of calculated vibrational modes.

## <a id='reorgdat'></a>reorgdat: `g_reorgdat_d`
This tool re-organizes the numbering of each atom in a CNDO .dat file according to the ordering provided by user. This is particularly useful in Block-Diagonalization calculations.

## reorginp: `g_reorginp_d`
Similar to [reorgdat](#reorgdat) but performs on G09 input file.

## rmsolshell: `g_rmsolshell_d`
It takes GROMACS geometry file (.gro file) to remove a shell of solvent molecules. The solvent, solute and the thickness of solvent shell are defined by user.

## rotateMol: `g_rotateMol_d`
Rotates the entire configuration based on user's choice. This tool will align the (donor group COM)-to-(acceptor group COM) vector to +z axis and the (donor group COM)-to-(L shape atom) vector to xOz plane.


# Credits

All the C code are written by [Zheng Ma](https://github.com/zanemarkson/) and part of the Fortran code is contributed by [Dr. Peng Zhang](https://scholars.duke.edu/display/per8149892) as stated above. 



