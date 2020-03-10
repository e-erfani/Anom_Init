# Anom_Init

Motivation:

Anomaly Initialization is a method to improve the seasonal/decadal prediction of a global climate model, especially the ocean oscillators (that are important in long-term predictions). 

In this experiment, we use a common climatology (CMIP6) and add an anomaly called "optimal perturbation" that would excite the ocean oscillator in Atlantic Ocean. We use a global climate model called CESM2.


Methods:

Creating modified CESM2 initial files for optimal perturbation project.

Copy the current directory for yourself, and go to that directory.
Make sure you installed Python and the required libraries (numpy, matplotlib, netCDF4).

Part I:
Creating modified CESM2 initial files.

Two python codes will be used. You don't need to change any address or directory. To apply different factors (1, 5, 10, 20, ...), just modify the variable "factor" at the beginning of the python codes.

1- TEMP_anom_init.py:
This code multiplies the optimal perturbation by the factor, and adds that to CMIP6 climatology
inputs: tntl_interp.mat (optimal perturbation), TEMP_mean_ens_cmip6_interp_1JAN.mat (CMIP6 climatology)
outputs: <factor>_1990_TEMP_CUR.pkl, <factor>1990_TEMP_CUR.pkl (<factor> can be 1, 5, ...), and one map .eps file (for verification)

run this code:
ipython TEMP_anom_init.py


Create a copy of POP initial file from the original one.
Be careful not to overwrite!
If there is already a POP initial file in this folder, make sure to copy it in another directory before envoking the below command line:

cp b.e20.BHIST.f09_g17.20thC.297_01.pop.r.1990-01-01-00000_ORIGINAL.nc b.e20.BHIST.f09_g17.20thC.297_01.pop.r.1990-01-01-00000.nc

2- BHIST_YEAR_01_01_output_anom_init.py
This code replaces the original TEMP_CUR and TEMP_OLD in the CESM2 POP initial file with the mosidied ones (peturbation + climatology) produced in the previous step.
inputs: b.e20.BHIST.f09_g17.20thC.297_01.pop.r.1990-01-01-00000_ORIGINAL.nc, <factor>_1990_TEMP_CUR.pkl, <factor>_1990_TEMP_CUR.pkl
outputs: b.e20.BHIST.f09_g17.20thC.297_01.pop.r.1990-01-01-00000_ORIGINAL.nc, and two map .eps files (for verification)

run this code:
ipython BHIST_YEAR_01_01_output_anom_init.py


