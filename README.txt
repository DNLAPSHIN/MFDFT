The initial code was provided by Peter Monson (University of Massachusetts, Amherst).
The provided code is modified by Dmitry Lapshin (University of Edinburgh) 2021.

These instructions include description of the files used in MFDFT simulations and provide a step-by-step procedure
to run simulations.
____________________

(1)
MFDFT simulations condected in the paper consist of 2 steps that correspond to two respective folders:
1. 'IL_distribution_NVT'
2. 'N2 adsorption_muVT'
____________________


(2)
Folder 'IL_distribution_NVT' contains files to run simulations in the canonical ensemble.
The key initial condition here is that the desnity of fluid is constant.

(2.1)
Folder 'Files_to_run_simulations'

(2.1.1)
Before running the code it should be compiled using

module add intel
ifort -O2 -o mft.exe mft_nvt.f

After this, file 'mft.exe' will be generated. Next, run the simulation using

./mft.exe


(2.1.2)
File 'mft_nvt.dat' contains initial conditions and simulation parameters.
File 'mft_nvt.f' contains the MFT code.

(2.1.3)
File 'mft_psd.dat' is the PSD in z-direction. This PSD is divided into two parts (2 lines) relative to the z-axis.
File 'Pore.py' is a python script used to generate the structure of the pore  and written to 'mft_psd.dat'.
File 'Pore.py' is used to generate a single mesopore with varying widths of intrawall pores.

(2.2)
Folder 'Example_of_obtained_results'

(2.2.1)
'isotherm_wpore' and 'isotherm_mpore' contains the following columns: tstar, cp, exp((cp+3.0d0)/tstar), rhoav (rhoav1), gp
tstar - Temperature
cp - chemical potential
exp((cp+3.0d0)/tstar) - relative activity
rhoav - average density of fluid for the whole pore
rhoav1 - average density of fluid for the middle of the pore
gp - grand potential

(2.2.2)
'rho001.dat' and 'rho314.dat' contain information about the pore and the adsorbed fluid (in our case ionic liquid).
These files are used further to generate a new pore containing IL.

Use RasMol to visualize the pore and the adsorbed fluid.
____________________


(3)
Folder 'N2_adsorption_muVT' contains files to run simulations in the grand canonical ensemble.

(3.1)
Folder 'Files_to_run_simulations'

(3.1.1)
Change the name of the file 'rho314.dat' to '314.dat' and run a python script 'PSD.py' using

./python 314

(3.1.2)
File 'PSD.py' is a python script which uses 'rho314.dat' to generate a new pore containing the ionic liquid
adsorbed from the previous simulation in the canonical ensemble.

The script generates file 'mft_psd.dat' used for N2 adsorption simulations.
The script also generates file '314_mft_psd_count' which counts the number of pores with different widths,
in other words this is pore size distribution for the pore containing the ionic liquid. The file contains:
1 column - width of the pore
2 column - the number of pores with the corresponding width.
Last line - the number of empty lattice sites i.e. sites not filled with the ionic liquid. If the pore is
completely filled with the ionic liquid this number is 0 (zero).

(3.1.3)
Files 'mft2.dat' and 'mft2.f' are simular to those discussed above but changed for simulation in the grand
canonical ensemble.