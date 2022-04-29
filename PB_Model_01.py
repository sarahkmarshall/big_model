# -*- coding: utf-8 -*-
"""
Created on Tue Jul  6 09:54:26 2021

@author: 00098687
"""

import numpy as np
import pylab as plt
import flopy
import flopy.utils.binaryfile as bf
from matplotlib.colors import LinearSegmentedColormap
import pandas as pd
import pickle
import os
import time


#------------------------------------------------------------------------------
# Options for this script

firstrun = 'no' #'yes' 'no'

#------------------------------------------------------------------------------
# Define classes and functions

class ASCII:
    
    def __init__(self, fname):
        f = open(fname,'r')
        l = f.readline().strip().split()
        self.ncols = int(l[1])
        l = f.readline().strip().split()
        self.nrows = int(l[1])
        l = f.readline().strip().split()
        self.xllc = float(l[1])
        l = f.readline().strip().split()
        self.yllc = float(l[1])
        l = f.readline().strip().split()
        self.cellsize = float(l[1])
        l = f.readline().strip().split()
        self.NODATA = float(l[1])
        f.close()
        self.data = np.loadtxt(fname,skiprows = 6)
        #self.data = np.ma.masked_where(self.data == self.NODATA, self.data)


#------------------------------------------------------------------------------
# Set up directories
        
os.getcwd()
this_model_fldr = r'C:\Users\mar886\WaterTableProject\Bens_Oasis\Pilbara'
os.chdir(this_model_fldr)
  
figureDirectory = os.path.join(this_model_fldr, "Figures")
dataDirectory = os.path.join(this_model_fldr, "data")

#------------------------------------------------------------------------------

# Data set up

     
h_dat = pd.read_csv('Head_data.txt', sep = '\t') # This is the data from Rio for the first head occurrence
rc_df_path = os.path.join(r'C:\Users\mar886\WaterTableProject\Bens_Oasis\Pilbara\data', 'rc_df.csv')

rc_df_already_saved = 'yes'

if rc_df_already_saved == 'yes':
    rc_df = pd.read_csv(rc_df_path)
else:
    
    rc_df = h_dat.copy()
    rc_df.head()
    rc_df.columns
    rc_df = rc_df.drop(labels=['Sample Date'], axis=1)
    rc_df = rc_df.rename(columns={'(Calc) LEVEL WATER (COMB) (mAHD)': "rswl",
                                  'Sample Point': "well_name"})
    rc_df.head()

    #Taken from code below
    '''
    rc_df["layer"] = jvals
    rc_df["column"] = nvals
    rc_df["row"] = mvals
    '''

    rc_df.to_csv(rc_df_path, encoding='utf-8')






#------------------------------------------------------------------------------
# Setting up the geological model

#####################Start importing the leapfrog surfaces####################
fname = r'./uwa_072021\FOR - MM Contacts.asc'
A = ASCII(fname)
plt.figure()
#plt.imshow(np.ma.masked_where(A.data => 20,A.data))

fname = r'./uwa_072021\MM - WF Contacts.asc'
B = ASCII(fname)

fname = r'./uwa_072021\WF - MCS Contacts.asc'
C = ASCII(fname)

fname = r'./uwa_072021\MCS - BRK Contacts.asc'
D = ASCII(fname)

fname = r'./uwa_072021\BRK - WW Contacts.asc'
E = ASCII(fname)



###################### End leapfrof surface importation#######################

# Values for colourbars
a = np.ma.masked_where(A.data == A.NODATA,A.data)
minv = np.nanmin(a)
maxa = np.nanmax(a)


# PLOT FIGURE FOR LEAPFROG SURFACES

fig = plt.figure(figsize=(10,15))

ax1 = plt.subplot(511)
im1 = plt.imshow(np.ma.masked_where(A.data == A.NODATA,A.data),
                 vmin=minv, vmax=maxa, cmap="gist_earth")
ax1.axes.get_xaxis().set_visible(False)
ax1.text(0.5, -5, "(a) FOR - MM contacts")

#cb1 = plt.colorbar(im1, orientation="vertical", format='%1.0f')

# - - - - - - - - - - - 
ax2 = plt.subplot(512)
im2 = plt.imshow(np.ma.masked_where(B.data == B.NODATA,B.data),
                 vmin=minv, vmax=maxa, cmap="gist_earth")
ax2.axes.get_xaxis().set_visible(False)
ax2.text(0.5, -5, "(b) MM - WF contacts")

#cb2 = plt.colorbar(im2, orientation="vertical", format='%1.0f')


# - - - - - - - - - - - 
ax3 = plt.subplot(513)
im3 = plt.imshow(np.ma.masked_where(C.data == C.NODATA,C.data),
                 vmin=minv, vmax=maxa, cmap="gist_earth")
ax3.axes.get_xaxis().set_visible(False)
ax3.text(0.5, -5, "(c) WF - MCS Contacts")

#cb3 = plt.colorbar(im3, orientation="vertical", format='%1.0f')


# - - - - - - - - - - - 
ax4 = plt.subplot(514)
im4 = plt.imshow(np.ma.masked_where(D.data == D.NODATA,D.data),
                 vmin=minv, vmax=maxa, cmap="gist_earth")
ax4.axes.get_xaxis().set_visible(False)
ax4.text(0.5, -5, "(d) MCS - BRK Contacts")

#cb4 = plt.colorbar(im4, orientation="vertical", format='%1.0f')


# - - - - - - - - - - - 
ax5 = plt.subplot(515)
im5 = plt.imshow(np.ma.masked_where(E.data == E.NODATA,E.data),
                 vmin=minv, vmax=maxa, cmap="gist_earth")
ax5.axes.get_xaxis().set_visible(False)
ax5.text(0.5, -5, "(e) BRK - WW Contacts")

#cb5 = plt.colorbar(im5, orientation="vertical", format='%1.0f')


# - - - - - - - - - - - 
plt.subplots_adjust(bottom=0.21)


cbaxes1 = fig.add_axes([0.133, 0.1, 0.32, 0.03]) 
cb1 = plt.colorbar(im1, cax = cbaxes1, orientation="horizontal", format='%1.0f')

plt.savefig(os.path.join("Figures", "leapfrog_surfaces"), dpi=300)  


#------------------------------------------------------------------------------
# Model grid set up

# Set up the grid and map the geology based on the leapfrog surfaces

z = np.arange(-220,580,20)
z = np.append(z,np.arange(580,700,10))

# Set up an empty grid
geology = np.zeros((len(z), A.nrows, A.ncols))

i=0
j=0
for i in range(A.ncols):
    for j in range(A.nrows):
        if A.data[j,i] > A.NODATA:
            n = np.argmin(abs(z - A.data[j,i]))
            geology[n:,j,i] = 1.
        else:
            geology[:,j,i] = 1.
        if B.data[j,i] > B.NODATA:
            n = np.argmin(abs(z - B.data[j,i]))
            geology[n:,j,i] = 2.
        elif A.data[j,i] == A.NODATA:
            geology[:,j,i] = 2.
        if C.data[j,i] > C.NODATA:
            n = np.argmin(abs(z - C.data[j,i]))
            geology[n:,j,i] = 3.
        elif A.data[j,i] == A.NODATA and  B.data[j,i] == B.NODATA:
            geology[:,j,i] = 3.
        if D.data[j,i] > D.NODATA:
            n = np.argmin(abs(z - D.data[j,i]))
            geology[n:,j,i] = 4.
        elif A.data[j,i] == A.NODATA and  B.data[j,i] == B.NODATA and  C.data[j,i] == C.NODATA:
            geology[:,j,i] = 4.
        if E.data[j,i] > E.NODATA:
            n = np.argmin(abs(z - E.data[j,i]))
            geology[n:,j,i] = 5.
        '''elif A.data[j,i] == A.NODATA and  B.data[j,i] == B.NODATA and  C.data[j,i] == C.NODATA and D.data[j,i] == D.NODATA:
            geology[:,j,i] = 5.'''
        if A.data[j,i] != A.NODATA and  B.data[j,i] == B.NODATA and  C.data[j,i] == C.NODATA and D.data[j,i] == D.NODATA and E.data[j,i] == E.NODATA:
            geology[:,j,i] = 0.

            
####there is one quirk I eed to fix so this just brute forces it
geology[:,100:,:235][geology[:,100:,:235]==4]= 0.

################Now we have the surfaces mapped to the grid##################
# Plot the geology for all layers     


aq_names = ["0", "1", "2", "3", "4", "5"]

aq_hsu   = [0,1,2,3,4,5]

aq_colours = ["slategrey", # 0
          (1.,204/255.,0.), # 1
          (160/255.,90/255.,44/255.), #2
          "cornflowerblue", # 3
          (1.,1/5.,0.), # 4
          (0.,90/265.,44/255.)] # 5

aq_cm = LinearSegmentedColormap.from_list('aq_list', aq_colours, N=len(aq_hsu))

lay=0
'''
for lay in range(len(z)):
    plt.figure()
    img = plt.imshow(geology[lay], cmap=aq_cm, vmin=-0.5, vmax=5.5)     
    cbr = plt.colorbar(img, ticks=[0, 1, 2, 3, 4, 5], orientation='vertical',shrink=0.5)
    cbr.ax.set_yticklabels(['No data', 'A','B','C', 'D', 'E', ])  # horizontal colorbar
    
    
    plt.savefig(os.path.join("Figures", "Geology", "geol_lay_%i"%lay), dpi=300)  
    plt.close()
    print("Plotted layer: %i" %lay)
   ''' 
#-------------------------------------------------------------------------------
    
#simultamneously map the higher res DEM to determine land surface, and add the detritals.
LSE = ASCII(r'./uwa_072021\Smooth_DEM_2_v2.asc')

x1 = np.arange(LSE.xllc,LSE.xllc+(LSE.ncols+1)*LSE.cellsize,LSE.cellsize)
y1 = np.arange(LSE.yllc,LSE.yllc+(LSE.nrows+1)*LSE.cellsize,LSE.cellsize)[::-1]
x2 = np.arange(A.xllc,A.xllc+(A.ncols+1)*A.cellsize,A.cellsize)
y2 = np.arange(A.yllc,A.yllc+(A.nrows+1)*A.cellsize,A.cellsize)[::-1]


LSE2 = np.zeros_like(A.data)
det_bot = np.zeros_like(A.data)

new_det = np.loadtxt('new_detritals.dat')[::-1,:] #This is my updated detritals
surface_geo = np.zeros_like(det_bot) #Here I will just track the surface geo for a figure...
for i in range (A.nrows):
    for j in range(A.ncols):
        m = np.argmin(abs(x2[j]-x1))
        n = np.argmin(abs(y2[i]-y1))
        LSE2[i,j] = LSE.data[n,m]
        #det_bot[i,j] = LSE2[i,j] - det_thick.data[nn,mm]
        if new_det[i,j] > 0.1:
            det_bot[i,j] = LSE2[i,j] - new_det[i,j]
            n = np.argmin(abs(z - det_bot[i,j]))
            geology[n:,i,j] = 6.
        n = np.argmin(abs(z - LSE2[i,j]))
        if z[n] < LSE2[i,j]:
            n+= 1
        geology[n:,i,j] = -1.
        surface_geo[i,j] = geology[n-1,i,j]
        if A.data[i,j] != A.NODATA:
            zz = np.min([A.data[i,j],z[-1],LSE2[i,j]])
            zdum = min([zz-20,550.])
            n = np.argmin(abs(zdum - z))
            geology[:n,i,j] = -1
        else:
            #find the top of the basement
            dum = -1
            for k in range(len(z)):
                if geology[k,i,j] == 0:
                    dum = k
            if dum > 0:
                zz = np.min([z[k],LSE2[i,j]])
                zdum = min([zz-20,550.])
                n = np.argmin(abs(zdum - z))
                geology[:n,i,j] = -1                    
                    
                
################ Finish land surface a detrital implementation################        


###### The next section adds "activated" dykes to the model. Only layers below detritals ##################          
################################
Dykes = True #controls 2 primary Dykes
Dyke3 = False #controls Dyke 3
Dyke4 = False  #controls Dyke 4
Dyke5 = False  #controls Dyke 5
################################
if Dykes:
    Dloc = np.loadtxt('Dykes.dat')
    for k in range(np.shape(geology)[0]):
        for i in range(A.ncols):
            for j in range(A.nrows):
                if geology[k,j,i] >=0 and geology[k,j,i] < 6:
                    if Dloc[j,i] > 0 and Dloc[j,i]<=2:
                        geology[k,j,i] = 7.
if Dyke3:
    Dloc = np.loadtxt('Dykes.dat')
    for k in range(np.shape(geology)[0]):
        for i in range(A.ncols):
            for j in range(A.nrows):
                if geology[k,j,i] >=0 and geology[k,j,i] < 6:
                    if Dloc[j,i] ==3:
                        geology[k,j,i] = 7.

if Dyke4:
    Dloc = np.loadtxt('Dykes.dat')
    for k in range(np.shape(geology)[0]):
        for i in range(A.ncols):
            for j in range(A.nrows):
                if geology[k,j,i] >=0 and geology[k,j,i] < 6:
                    if Dloc[j,i] ==4:
                        geology[k,j,i] = 7.
                        
if Dyke5:
    Dloc = np.loadtxt('Dykes.dat')
    for k in range(np.shape(geology)[0]):
        for i in range(A.ncols):
            for j in range(A.nrows):
                if geology[k,j,i] >=0 and geology[k,j,i] < 6:
                    if Dloc[j,i] ==5:
                        geology[k,j,i] = 7.


##################Added detritals#########################################
# Plot the geology for all layers 
                        
# How many unique geology entries? Check
geo_vals = np.unique(geology)
print(geo_vals)

aq_names = ["Inactive", "Basement", "MM", "Wittenoom", "McRae shale", "Brockman", 
            "WW", "Detritals", "Dyke"]

aq_hsu   = geo_vals

aq_colours = ["white", # -1
              "k", # 0
              (1.,204/255.,0.), # 1
              "silver", #2
              "cornflowerblue", # 3
              (1.,1/5.,0.), # 4
              (0.,90/265.,44/255.), # 5
              "lightgoldenrodyellow", # Detritals
              "teal"] # Dykes 

aq_cm = LinearSegmentedColormap.from_list('aq_list', aq_colours, N=len(aq_hsu))

# PLOTTING - don't want to do this every time
'''
lay=0
for lay in range(len(z)):
    plt.figure(figsize=(20,10))
    img = plt.imshow(geology[lay], cmap=aq_cm, vmin=-1.5, vmax=7.5)     
    cbr = plt.colorbar(img, ticks=aq_hsu, orientation='vertical',shrink=0.5)
    cbr.ax.set_yticklabels(aq_names)  # horizontal colorbar
    
    
    plt.savefig(os.path.join("Figures", "Geology_with_dykes", "geol_lay_%i"%lay), dpi=300)  
    plt.close()
    print("Plotted layer: %i" %lay)
    
    
# Plot some cross sections through the model domain
row=0
for row in range(0,A.nrows,4):
    plt.figure(figsize=(20,10))
    img = plt.imshow(geology[:,row,:], cmap=aq_cm, vmin=-1.5, vmax=7.5)     
    cbr = plt.colorbar(img, ticks=aq_hsu, orientation='vertical',shrink=0.5)
    cbr.ax.set_yticklabels(aq_names)  # horizontal colorbar
    
    
    plt.savefig(os.path.join("Figures", "Geology_with_dykes", "geol_row_%i"%row), dpi=300)  
    plt.close()
    print("Plotted row: %i" %row)    
   ''' 
# -----------------------------------------------------------------------------

geology=geology[::-1,:,:] # just to re-orientat geology mod for Modflow

       
#Ibound array, 1 for active, zero for inactive, -1 constant head
IBOUND = np.ones_like(geology,dtype = int)
IBOUND[geology == -1.] = 0
#IBOUND[geology == 0.] = 0
Extent = np.loadtxt(r'Extent.dat') ##### This is just the active model domain.
#print(dkjhkfk)
#print(jhsdf) 
for i in range(np.shape(geology)[0]):
    IBOUND[i,:,:][Extent == 0] = 0
    
#Determine head locations for the estern boundary. Could make this general head?
hd = np.zeros_like(Extent)
for j in range(110,174,1):
    for i in range(390,np.shape(Extent)[1]-1):
        if Extent[j,i] == 1 and Extent[j,i+1] == 0:
            hd[j,i] = 1
        

# Save the geology file
outputfp = os.path.join(dataDirectory, "geology.p")
with open(outputfp, "wb") as fp:
    pickle.dump(geology, fp)
    
#Layer property data
kxx = np.ones_like(geology)*1e-4

Kh_detritals = 10.
Kh_WW = 0.1
Kh_Brockman = 0.001
Kh_Mcrae = 0.001
Kh_Wittenoom = 5.0
Kh_MM = 0.1
Kh_Basement = 0.5
Kh_Dyke = 1e-4
 
kKh = np.ones_like(geology)*1e-4
kKh[geology == 7.] = Kh_Dyke
kKh[geology == 6.] = Kh_detritals
kKh[geology == 5.] = Kh_WW
kKh[geology == 4.] = Kh_Brockman
kKh[geology == 3.] = Kh_Mcrae
kKh[geology == 2.] = Kh_Wittenoom
kKh[geology == 1.] = Kh_MM
kKh[geology == 0.] = Kh_Basement

# Save the kKh 
outputfp = os.path.join(dataDirectory, "kKh.p")
with open(outputfp, "wb") as fp:
    pickle.dump(kKh, fp)


# Sy
Sy_detritals = 0.1
Sy_WW = 0.001
Sy_Brockman = 0.001
Sy_Mcrae = 0.001
Sy_Wittenoom = 0.003
Sy_MM = 0.001
Sy_Basement = 0.0005

ksy = np.ones_like(geology)*1e-4
ksy[geology == 6.] = Sy_detritals
ksy[geology == 5.] = Sy_WW
ksy[geology == 4.] = Sy_Brockman
ksy[geology == 3.] = Sy_Mcrae
ksy[geology == 2.] = Sy_Wittenoom
ksy[geology == 1.] = Sy_MM
ksy[geology == 0.] = Sy_Basement

# Save the ksy
outputfp = os.path.join(dataDirectory, "ksy.p")
with open(outputfp, "wb") as fp:
    pickle.dump(ksy, fp)

SS_detritals = 1e-4
SS_WW = 1e-5
SS_Brockman = 1e-5
SS_Mcrae = 2e-4
SS_Wittenoom = 3.1e-4
SS_MM = 1e-5
SS_Basement = 2e-4

kss = np.ones_like(geology)*1e-4
kss[geology == 6.] = SS_detritals
kss[geology == 5.] = SS_WW
kss[geology == 4.] = SS_Brockman
kss[geology == 3.] = SS_Mcrae
kss[geology == 2.] = SS_Wittenoom
kss[geology == 1.] = SS_MM
kss[geology == 0.] = SS_Basement

# Save the kss
outputfp = os.path.join(dataDirectory, "kss.p")
with open(outputfp, "wb") as fp:
    pickle.dump(kss, fp)

vka_detritals = Kh_detritals/10.
vka_WW = Kh_WW/10.
vka_Brockman = Kh_Brockman/10
vka_Mcrae = Kh_Mcrae/10.
vka_Wittenoom = Kh_Wittenoom/10.
vka_MM = Kh_MM/10
vka_Basement = Kh_Basement/10.

kvk = np.ones_like(geology)*1e-4
kvk[geology == 7.] = Kh_Dyke
kvk[geology == 6.] = vka_detritals
kvk[geology == 5.] = vka_WW
kvk[geology == 4.] = vka_Brockman
kvk[geology == 3.] = vka_Mcrae
kvk[geology == 2.] = vka_Wittenoom
kvk[geology == 1.] = vka_MM
kvk[geology == 0.] = vka_Basement


# Save the kKh 
outputfp = os.path.join(dataDirectory, "kKh.p")
with open(outputfp, "wb") as fp:
    pickle.dump(kKh, fp)

 
#laytyp = 7

plt.figure()
plt.imshow(np.log10(kxx[:,90,:]),origin = 'lower') 

plt.figure()
plt.imshow(np.flipud(geology[:,90,:]),origin = 'lower') 

#Define Stress period
#print(hsablahbcs)

#Define Modflow object    
mf = flopy.modflow.Modflow('BO', exe_name = 'MODFLOW-NWT_64.exe', version = 'mfnwt') 

dum = np.shape(geology)
nlay = dum[0]
nrow = dum[1]
ncol = dum[2]
delr = A.cellsize
delc = delr
botm = z-5
botm = botm[::-1]
top = z[-1]+5
nper = 1
perlen = 365.25 * np.ones((nper),dtype=np.float32)
perlen[0] = 0.0
steady = 0*np.ones((nper),dtype=np.bool)
steady[0] = True
nstp=12*np.ones((nper),dtype=np.int32)
nstp[0] = 1


#Define Modflow spatial discretisation
dis = flopy.modflow.ModflowDis(mf, nlay, nrow, ncol, delr = delr, delc= delc,
                           top= top , botm= botm, nper=nper, perlen=perlen, nstp=nstp,
                               steady=steady) 
ibound = IBOUND
'''ibound[:,:,0] = -1
ibound[:,:,-1] = -1'''
for i in range(nlay):
    ibound[i,:,:][hd == 1] = -1
ibound[geology == -1.] = 0
CH = np.loadtxt('CH.dat')
ibound[11,:,:][CH == 1] = -1

strt = np.ones_like(geology,dtype = float)
for i in range(nlay):
    strt[i,:,:] = LSE2
    strt[i,:,:][hd == 1] = 640


###############################################################################
    # FIRST RUN DON'T USE THIS
if firstrun == 'yes': # 'no'
    pass
else:
    #this just recycles the last head as the start.... If a model crashes, this won't be re-written
    hdobj = bf.HeadFile('BO.hds', precision='single')
    strt = hdobj.get_data(kstpkper=(0, 0))
    strt[11,:,:][CH == 1] = 570.

###############################################################################

#strt[:,:,0] = 2000.
#strt[:,:,-1] = 2010.

#Define BAS
bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=strt) 

#Building LPF for Modflow
#lpf = flopy.modflow.ModflowLpf(mf, hk= kKh,laywet=1, ihdwet=1, vka=kvk, ss=kss, sy=ksy, laytyp=1, ipakcb=53)
upw = flopy.modflow.ModflowUpw(mf, hk= kKh, vka=kvk, ss=kss, sy=ksy, laytyp=1, ipakcb=53)
'''pcg = flopy.modflow.ModflowPcg(mf, mxiter=100, iter1=50, 
                               hclose=1e-05, rclose=1e-03,damp=0.1)'''
#(model, laytyp=0, layavg=0, chani=1.0, layvka=0, laywet=0, ipakcb=None, hdry=-1e+30, iphdry=0, hk=1.0, hani=1.0, vka=1.0, ss=1e-05, sy=0.15, vkcb=0.0, noparcheck=False, extension='upw', unitnumber=None, filenames=None
nwt= flopy.modflow.ModflowNwt(mf,options='MODERATE', headtol = 1e-4, maxiterout=1000,  iprnwt=1)

##############Ephemeral streams
Det_rech = 0.1/365 #recharge rate to the detritals
Det_cond = 2000. #drain conductance of the detritals
wit_rech = 0.001/365 #recharge rate to the wittenoom
wit_cond = 1000. #drain conductance of the wittenoom
bed_rech = 0.001/365 #recharge to the bderock
bed_cond = 1000. #conductance of the bedrock
##############################################
# We will make drains and recharge rater based on this...
RECH = np.zeros((nrow,ncol))
a = np.loadtxt('Detrital_recharge.dat')
RECH[a==1.] = Det_rech
a = np.loadtxt('Wittenoom_recharge.dat')
RECH[a==1.] = wit_rech
a = np.loadtxt('Bedrock_recharge.dat')
RECH[a==1.] = bed_rech


# save the recharge file
outputfp = os.path.join(dataDirectory, "recharge.p")
with open(outputfp, "wb") as fp:
    pickle.dump(RECH, fp)

#### RECHARGE

rch = flopy.modflow.ModflowRch(mf, ipakcb=53,rech=RECH) 


deets = []
deets_rivr = []
a = np.loadtxt('Detrital_recharge.dat')
for j in range(nrow):
    for i in range(ncol):
        if a[j,i] == 1.:
            k = np.argmin(abs(botm-LSE2[j,i]))
            if LSE2[j,i] < botm[k]:
                k+=1
            deets.append([k,j,i,LSE2[j,i],Det_cond])
            deets_rivr.append([k,j,i,LSE2[j,i],Det_cond,LSE2[j,i]-0.5])
a = np.loadtxt('Wittenoom_recharge.dat')
for j in range(nrow):
    for i in range(ncol):
        if a[j,i] == 1.:
            k = np.argmin(abs(botm-LSE2[j,i]))
            if LSE2[j,i] < botm[k]:
                k+=1
            deets.append([k,j,i,LSE2[j,i],wit_cond])    
            deets_rivr.append([k,j,i,LSE2[j,i],wit_cond,LSE2[j,i]-0.5])
a = np.loadtxt('Bedrock_recharge.dat')
for j in range(nrow):
    for i in range(ncol):
        if a[j,i] == 1.:
            k = np.argmin(abs(botm-LSE2[j,i]))
            if LSE2[j,i] < botm[k]:
                k+=1
            k+=1
            deets.append([k,j,i,LSE2[j,i],bed_cond])  
            deets_rivr.append([k,j,i,LSE2[j,i],bed_cond,LSE2[j,i]-0.5])
spd = {}

spd[0] = deets
drn = flopy.modflow.ModflowDrn(mf, ipakcb=53, stress_period_data=spd)
#############ET parameters
spd = {}

spd[0] = deets_rivr
#riv = flopy.modflow.ModflowRiv(mf,stress_period_data = spd)


#######ET package
if firstrun == 'yes': # 'no'
    print('Dont use any ET for first run')
else:
    evt = flopy.modflow.ModflowEvt(mf,nevtop=3, ipakcb=53, surf=LSE2, evtr=0.0001, exdp=3.0) # evtr=0.0011

############

spd = {(0,0): ['save head', 'save budget']}
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd)

#spd = {(0,0): print['print head', 'print budget']}
#oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd,
#compact=True)


# - - - - - - - -
# Adding headobj file (used for PEST)

obs = []

roff = 0.0 # ml.dis.delc.array[0] / 2.0
coff = 0.0 # ml.dis.delr.array[0] / 2.0
# roff is the fractional offset from the centre of the cell in Y direction (between rows)
# coff is the fractional offset from the centre of the cell in the X direction (between cols)

obs_name_format = "j{1:02d}m{1:03d}n{2:03d}"
n = rc_df.well_name[0]
r = rc_df.row[0]
c = rc_df.column[0]
l = rc_df.layer[0]
wl = rc_df.rswl[0]

for n, r, c, l, wl in zip(rc_df.well_name, rc_df.row, rc_df.column, rc_df.layer, rc_df.rswl):
    
    name = obs_name_format.format(n, l, r, c) 
    print(name)
    obs_head_name = str(name)+"_h"
    
    # Create empty data shape to put the head data into.
    d = np.zeros((1, 2)) # 1 is referring to one time period
    print(d)
    
    # Now put the times into the first column of empty matrix.
    d[:,0] = 0 # Time is 0    
    d[:, 1] = wl

    # Now create a list of each obs - these are zero based and flopy converts them to 1-based (for modflow)
    obs.append(flopy.modflow.mfhob.HeadObservation(mf, obsname=obs_head_name, 
                                                   layer=l, row=r,
                                                   column=c, roff=roff, coff=coff,
                                                   time_series_data=d)) 
    
##? Do I need to do this...because the model takes so long to run I need to optimise    
flopy.modflow.mfhob.ModflowHob(mf, obs_data=obs, iuhobsv=600) 

# - - - - - - - - - -

# # # # # # # # # # 
# RUNNING THE MODEL

t1 = time.time()  
mf.write_input()
success, mfoutput = mf.run_model(silent=True, pause=False)
if not success:
    raise Exception("MODFLOW did not terminate normally.")
t2 = time.time()  
print("Time to run model: %2.2f" %((t2-t1)/60))
# # # # # # # # # #

plt.close('all')
    
################################################################################
# END MODEL RUN
###############################################################################

# -----------------------------------------------------------------------------
# Setting up data products 

#wells = np.loadtxt('Bore_Static_WL.txt',skiprows = 1,usecols=[2,3,5])
# kstpkper (tuple of ints) – A tuple containing the time step and stress period (kstp, kper). The kstp and kper values are zero based.

hdobj = bf.HeadFile('BO.hds', precision='single')
heads_data = hdobj.get_alldata()
heads_data.shape # (time period, layer, row, column)

strt = hdobj.get_data(kstpkper=(0, 0))       

# Setting up dataset of modelled heads at well locations

i=0

# Get the indexes for the wells to use for PEST in rc_df
nvals = []
mvals = []
jvals = []

mod = np.zeros(len(h_dat))
for i in range(len(h_dat)):
    print(h_dat["Sample Point"][i])
    n = np.argmin(abs(h_dat['Easting'][i]-x2))
    m = np.argmin(abs(h_dat['Northing'][i]-y2))
    nvals.append(n)
    mvals.append(m)
    jval = 0
    for j in range(nlay-1,-1,-1):
        if strt[j,m,n] > 0.:
            mod[i] = strt[j,m,n]
            jval = j
    jvals.append(jval) # Getting top value
    
    
cbb = bf.CellBudgetFile('BO.cbc')
cbb.list_records()
rec = cbb.get_data(kstpkper=(0,0), text='DRAINS')
dloss = np.zeros((nrow,ncol))
for i in range(len(rec[0])):
    if rec[0][i][1] != 0.:
        col = rec[0][i][0]%(ncol*nrow)%ncol
        row = int(rec[0][i][0]%(ncol*nrow)/ncol)
        dloss[row,col] = rec[0][i][1]


rec = cbb.get_data(kstpkper=(0,0), text='ET')
ET = rec[0][1]

        
RMSE = (mod -h_dat['(Calc) LEVEL WATER (COMB) (mAHD)'])**2.
RMSE = RMSE[mod>0.]
RMSE = np.average(RMSE)
RMSE = np.sqrt(RMSE)


# -----------------------------------------------------------------------------
# Plotting figures


# - - - - - - - - -
# Figure properties
extent = [A.xllc,A.xllc+A.cellsize*A.ncols,A.yllc,A.yllc+A.cellsize*A.nrows]



# - - - - - - - - -

#---------------------------------------------------
# Plot the difference between modelled and measured heads

############
plt.figure(figsize=(15,10))
############
plt.suptitle("Modelled versus measured heads")
# Plotting starting heads masked where cells are inactive
plt.imshow(np.ma.masked_where((IBOUND[12,:,:]==0), 
                              strt[12,:,:]), 
                              extent = extent)
cbar = plt.colorbar()
cbar.set_label('Head (m)')

for i in range(len(mod)):
    if mod[i] > 0.:
        if (mod[i] - h_dat['(Calc) LEVEL WATER (COMB) (mAHD)'][i]) < -5:
            plt.plot(h_dat['Easting'][i],h_dat['Northing'][i],'ro')
        elif (mod[i] - h_dat['(Calc) LEVEL WATER (COMB) (mAHD)'][i]) < 0.:
            plt.plot(h_dat['Easting'][i],h_dat['Northing'][i],'go')      
        else:
            plt.plot(h_dat['Easting'][i],h_dat['Northing'][i],'bo')
        plt.text(h_dat['Easting'][i],h_dat['Northing'][i],np.round(mod[i] 
        - h_dat['(Calc) LEVEL WATER (COMB) (mAHD)'][i],1))

# Make a little legend
plt.plot(h_dat['Easting'][i]+10000,h_dat['Northing'][i],'ro', 
         label="Modelled heads much less than measured")
plt.plot(h_dat['Easting'][i]+10000,h_dat['Northing'][i],'go',
         label="Modelled heads somewhat less than measured")      
plt.plot(h_dat['Easting'][i]+10000,h_dat['Northing'][i],'bo',
         label="Modelled heads greater than measured")


ax=plt.gca()
ax.set_xlim(extent[0], extent[1])
ax.set_ylim(extent[2], extent[3])

plt.legend()

name = os.path.join(figureDirectory, "mdl_vs_measured_heads")
plt.savefig(name, dpi = 500)




#---------------------------------------------------
# Plot transect 1

############
plt.figure()
############
locs = np.loadtxt('T1_locs.dat')  
T1_h = np.zeros(np.shape(locs)[0])

# Write transect data into a file and save
f = open('T1_modelled.dat','w')
f.write('Distance (m)    head(mAHD) \n')

for i in range(np.shape(locs)[0]):
    T1_h[i] = strt[11, int(locs[i,0]), int(locs[i,1])]
    f.write('%g %g \n' % (locs[i,2],T1_h[i]))
plt.plot(locs[:,2],T1_h)
f.close()

locs = np.loadtxt('T2_locs.dat')  
T1_h = np.zeros(np.shape(locs)[0])
f = open('T2_modelled.dat','w')
f.write('Distance (m)    head(mAHD) \n')
for i in range(np.shape(locs)[0]):
    T1_h[i] = strt[11,int(locs[i,0]),int(locs[i,1])]
    f.write('%g %g \n' % (locs[i,2],T1_h[i]))
plt.plot(locs[:,2],T1_h)
f.close()

locs = np.loadtxt('T3_locs.dat')  
T1_h = np.zeros(np.shape(locs)[0])
f = open('T3_modelled.dat','w')
f.write('Distance (m)    head(mAHD) \n')
for i in range(np.shape(locs)[0]):
    T1_h[i] = strt[11,int(locs[i,0]),int(locs[i,1])]
    f.write('%g %g \n' % (locs[i,2],T1_h[i]))
plt.plot(locs[:,2],T1_h)
f.close()

locs = np.loadtxt('T4_locs.dat')  
T1_h = np.zeros(np.shape(locs)[0])
f = open('T4_modelled.dat','w')
f.write('Distance (m)    head(mAHD) \n')
for i in range(np.shape(locs)[0]):
    T1_h[i] = strt[11,int(locs[i,0]),int(locs[i,1])]
    f.write('%g %g \n' % (locs[i,2],T1_h[i]))
plt.plot(locs[:,2],T1_h)
f.close()
