###########################################
## Created by Ehsan Erfani
##
## This python code reads the Florian Optimal perturbation, 
## multiply that by the factor, add the result to CMIP6 climatology,
## and create a plot for verification
###########################################


from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import sys, os, glob, math, pickle, pprint

######################
# Define parameters
# Read lat, lon, tarea 
######################
factor = 10  # multiply the Florian optimal perturbation by this factor
lev = 60     # vertical level
target_year = 1990 # initial year

fle_grid = Dataset('file_test_f09_g16.nc') #just a sample file to read lat, lon, TAREA
lat = fle_grid.variables['TLAT'][:]
lon = fle_grid.variables['TLONG'][:]
TAREA = fle_grid.variables['TAREA'][:]

###########################################
#######  Perturbation from Florian  #######
####### based on Germe et al. (2017) ######
###########################################
import h5py

# read the perturbation in .mat format. This file was prepared by Matlab. 
# The perturbation is already interpolated vertically and horizontally to be consistent with CESM ocean products.
pert_file = 'tntl_interp.mat'
file = h5py.File(pert_file, 'r')
pert_list = list(file['tntl_interp']) 
pert = np.asarray(pert_list, dtype=np.float32)
pert = pert * factor  # multiply the perturbation by factor 

# Calculate the global mean of perturbation to make sure it makes sense
II = np.where(pert[0,:,:] == 0.0)
JJ = np.where(np.isnan(pert[0,:,:]) == 1)
TAREA2 = TAREA
TAREA2[JJ] = np.nan
TAREA2[II] = np.nan
pert_glb = np.nansum(np.nansum(pert[0,:,:] * TAREA2)) / np.nansum(np.nansum(TAREA2))
print(pert_glb)

# Make sure there is no NaN in the perturbation, because such NaN values will ultimately make the CESM stop running.
KK = np.where(np.isnan(pert) == 1)
pert[KK] = 0.0

###################
##### CMIP6 #######
###################

# Read the CMIP6 climatology. 
# The .mat CMIP6 file was previously saved based on the average of CMIP6 11 ensemble members.
filepath = 'TEMP_mean_ens_cmip6_interp_1JAN.mat'
file = h5py.File(filepath, 'r')
clim_list = list(file['TEMP_mean_ens_interp_1JAN'])
clim = np.asarray(clim_list, dtype=np.float32)

# Calculate the global mean of climatology to make sure it makes sense
II = np.where(clim[0,:,:] == 0.0)
JJ = np.where(np.isnan(clim[0,:,:]) == 1)
TAREA3 = TAREA
TAREA3[JJ] = np.nan
TAREA3[II] = np.nan
clim_glb = np.nansum(np.nansum(clim[0,:,:] * TAREA3)) / np.nansum(np.nansum(TAREA3))
print(clim_glb)

########################################
##### perturbation + climatology #######
########################################
pert_clim = np.empty((lev,lat.shape[0],lat.shape[1])) # define the output
pert_clim[:] = 0.0
pert_clim[:] = pert[:] + clim[:]

# Calculate the global mean of climatology to make sure it makes sense
KK = np.where(pert_clim[0,:,:] == 0.0)
TAREA4 = TAREA
TAREA4[KK] = np.nan
pert_clim_glb = np.nansum(np.nansum(pert_clim[0,:,:] * TAREA4)) / np.nansum(np.nansum(TAREA4))

######### save the output ( pert + clim) in a file ########
# The factor, initial year, and variable name are mentioned in the file name.
# Both CUR and OLD variables will be used in the CESM ocean initial file.
out_name = "{}_1990_TEMP_CUR.pkl".format(str(factor)) 
output = open(out_name, 'wb')
pickle.dump(pert_clim[:], output)
output.close()

out_name2 = "{}_1990_TEMP_OLD.pkl".format(str(factor))
output2 = open(out_name2, 'wb')
pickle.dump(pert_clim[:], output2)
output2.close()

######################################################
# The rest is just making a plot for verification:
######################################################

###### First, modify POP variables and lat and lon to be compatible for plotting
lon = np.where(np.greater_equal(lon,min(lon[:,0])),lon-360,lon) # make longitudes monotonically increasing.
# stack grids side-by-side (in longitiudinal direction), so any range of longitudes may be plotted on a world map.
lon = np.concatenate((lon,lon+360),1)
lat = np.concatenate((lat,lat),1)
lon = lon-360.
pert_clim = np.concatenate((pert_clim,pert_clim),2)
pert      = np.concatenate((pert,pert),2)

##############################
## plot the map of pert + clim
plt.figure(figsize=(13,7))
map = Basemap(projection='cyl',llcrnrlat=-75,urcrnrlat=75,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawcoastlines(linewidth=.5)
parallels = np.arange(-90,91,30.) 
meridians = np.arange(-180,180,60.)
matplotlib.rcParams['font.size'] = 12
map.drawparallels(parallels,labels=[1,0,0,0],linewidth=0.0,fontsize=12)
map.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.0,fontsize=12)
x,y = map(lon,lat)
clevs = np.arange(-5.,33.,1)
csf = map.contourf(x,y,pert_clim[0,:,:],clevs,extend='both',cmap='PuRd',vmin=-5., vmax=33.)
cb = map.colorbar(csf,"right", extend='neither',size="3%", pad="1%",ticks=[-5,0,5,10,15,20,25,30,35])
title_name = "SST, Perturbation (Florian, 1 Jan {}) + climatology (CMIP6 1 Jan 1960-2010),  mean = ".format(str(target_year))
plt.title(title_name + str(pert_clim_glb))
#plt.show()
file_name = "TEMP_1_pert_CMIP6_1990_{}.eps".format(str(factor))
plt.savefig(file_name , format='eps', dpi=1000)

