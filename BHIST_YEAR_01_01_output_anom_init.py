########################################################################
## Created by Ehsan Erfani
##
## This python code replaces the original TEMP_CUR and TEMP_OLD in the CESM
## POP initial file with the modified ones (perturbation + climatology),
## and creates some plots for verification
#######################################################################

from netCDF4 import Dataset
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.basemap import Basemap
import sys, os, glob, math, pickle, pprint

################
target_year = 1990 # initial year
factor = 10        # perturbation is multiplied by this factor

# Load temperature variables (perturbation + climatology)
name = "{}_1990_TEMP_CUR.pkl".format(str(factor))
pkl_file = open(name, 'rb')
temp_cur = pickle.load(pkl_file)
pkl_file.close()
name = "{}_1990_TEMP_OLD.pkl".format(str(factor))
pkl_file = open(name, 'rb')
temp_old = pickle.load(pkl_file)
pkl_file.close()

# replace NaN values with zero. NaN values will make the CESM stop.
I1 = np.where(np.isnan(temp_cur))
temp_cur[I1] = 0.0
I2 = np.where(np.isnan(temp_old))
temp_old[I2] = 0.0

################
# Create a copy of POP initial file from the original one.
# Be careful not to overwrite!
# If there is already a POP initial file in this folder, make sure to copy it in another directory before running the below line.
cp b.e20.BHIST.f09_g17.20thC.297_01.pop.r.1990-01-01-00000_ORIGINAL.nc b.e20.BHIST.f09_g17.20thC.297_01.pop.r.1990-01-01-00000.nc

#### Replace variables in CESM BHIST POP restart file with new variables
fle_name = 'b.e20.BHIST.f09_g17.20thC.297_01.pop.r.{}-01-01-00000.nc'.format(str(target_year))
fle_grid = Dataset(fle_name,'r+')
fle_grid.variables['TEMP_CUR'][:] = temp_cur[:]
fle_grid.variables['TEMP_OLD'][:] = temp_old[:]
fle_grid.close()

#########################################################################
#### The rest is verficiation to make sure the replacement is correct.

#### read lat, lon, and tarea from a sample file
filee = Dataset('file_test_f09_g16.nc')
lat = filee.variables['TLAT'][:]
lon = filee.variables['TLONG'][:]
TAREA = filee.variables['TAREA'][:]

#### read new CESM POP initial file
fle_grid = Dataset(fle_name)
TEMP_CUR = fle_grid.variables['TEMP_CUR'][:]
# calculate global mean for verification
II = np.where(TEMP_CUR[0,:,:] == 0.0)
TAREA2 = TAREA
TAREA2[II] = np.nan
var_glb = np.nansum(np.nansum(TEMP_CUR[0,:,:] * TAREA2)) / np.nansum(np.nansum(TAREA2))
print(var_glb)

##### Read the original CESM BHIST POP initial file
fle_name_orig = 'b.e20.BHIST.f09_g17.20thC.297_01.pop.r.{}-01-01-00000_ORIGINAL.nc'.format(str(target_year))
fle_orig = Dataset(fle_name_orig)
TEMP_CUR_orig = fle_orig.variables['TEMP_CUR'][:]
# calculate global mean of difference
JJ = np.where(TEMP_CUR_orig[0,:,:] == 0.0)
TAREA3 = TAREA
TAREA3[JJ] = np.nan
diff_glb = np.nansum(np.nansum((TEMP_CUR_orig[0,:,:] - TEMP_CUR[0,:,:]) * TAREA3)) / np.nansum(np.nansum(TAREA3))
print(diff_glb)

##############################################
######### Making a couple of plots for verification

### preparing variables for plotting
lon = np.where(np.greater_equal(lon,min(lon[:,0])),lon-360,lon) # make longitudes monotonically increasing.
# stack grids side-by-side (in longitiudinal direction), so any range of longitudes may be plotted on a world map.
lon = np.concatenate((lon,lon+360),1)
lat = np.concatenate((lat,lat),1)
lon = lon-360.
TEMP_CUR = np.concatenate((TEMP_CUR,TEMP_CUR),2)
TEMP_CUR_orig = np.concatenate((TEMP_CUR_orig,TEMP_CUR_orig),2)

#####################
plt.figure(figsize=(13,7))
map = Basemap(projection='cyl',llcrnrlat=-75,urcrnrlat=75,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawcoastlines(linewidth=.5)
parallels = np.arange(-90,91,30.) 
meridians = np.arange(-180,180,60.)
matplotlib.rcParams['font.size'] = 12
map.drawparallels(parallels,labels=[1,0,0,0],linewidth=0.0,fontsize=12)
map.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.0,fontsize=12)
clevs = np.arange(-5,32,2)
x,y = map(lon,lat)
csf = map.contourf(x,y,TEMP_CUR[0,:,:],clevs,extend='both',cmap='PuRd',vmin=-5, vmax=32)
cb = map.colorbar(csf,"right", extend='both',size="3%", pad="1%",ticks=[-5, 0, 5, 10, 15, 20, 25, 30])
TITLE = 'SST, perturbation (Florian 1 Jan {}) + climatology (CMIP6 1 Jan 1960-2010),  mean = '.format(str(target_year))
plt.title(TITLE + str(var_glb))
FILE_NAME = 'TEMP_CUR_1_FILE_pert_clim_1990_{}.eps'.format(str(factor))
plt.savefig(FILE_NAME, format='eps', dpi=1000)
plt.close()

################## 
plt.figure(figsize=(13,7))
map = Basemap(projection='cyl',llcrnrlat=-75,urcrnrlat=75,\
            llcrnrlon=-180,urcrnrlon=180,resolution='c')
map.drawcoastlines(linewidth=.5)
parallels = np.arange(-90,91,30.) 
meridians = np.arange(-180,180,60.)
matplotlib.rcParams['font.size'] = 12
map.drawparallels(parallels,labels=[1,0,0,0],linewidth=0.0,fontsize=12)
map.drawmeridians(meridians,labels=[0,0,0,1],linewidth=0.0,fontsize=12)
clevs = np.arange(-4,4,0.5)
x,y = map(lon,lat)
csf = map.contourf(x,y,(TEMP_CUR_orig[0,:,:] - TEMP_CUR[0,:,:]),clevs,extend='both',cmap='RdBu_r',vmin=-4, vmax=4)
cb = map.colorbar(csf,"right", extend='both',size="3%", pad="1%",ticks=[-4, -3, -2, -1, 0, 1, 2, 3, 4])
TITLE = 'SST, (ctrl orig) - (pert+clim), 1 Jan {}, mean = '.format(str(target_year))
plt.title(TITLE + str(diff_glb))
FILE_NAME = 'diff_origctrl_pertclim_TEMP_CUR_1_FILE_1990_{}.eps'.format(str(factor))
plt.savefig(FILE_NAME, format='eps', dpi=1000)
plt.close()

