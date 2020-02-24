import arcpy
import numpy as np
from arcpy import env
from arcpy.sa import *
import math


env.workspace = 'Y:\MP_research\Geospatial\MP_Geospatial\Data'
env.overwriteOutput = True

precip_file = 'nc_map'
temp_file = 'nc_mat'
PET_file = 'nc_pet_daily'
BD_file = 'nc_bd'
Texture_file = 'nc_soiltext'

#Make NumPy arrays from datasets
precip = arcpy.RasterToNumPyArray(precip_file)
temp = arcpy.RasterToNumPyArray(temp_file)
PET = arcpy.RasterToNumPyArray(PET_file)
BD = arcpy.RasterToNumPyArray(BD_file)
Texture = arcpy.RasterToNumPyArray(Texture_file)

mask = (Texture != 0)

#Soil hydraulics values
#'clay (light)' was originally sandy clay
b1 = {1:4.05, 2:4.38, 3:4.90, 4:5.30, 5:5.39, 6:7.12,\
      7:7.75, 8:8.52, 9:10.4, 10:10.4, 11:11.4}
thetas1 = {1:0.395, 2:0.441, 3:0.435, 4:0.485, 5:0.451, 6:0.42,\
      7:0.477, 8:0.476, 9:0.426, 10:0.492, 11:0.482}
Ks1 = {1:176e-6, 2:156e-6, 3:34.1e-6, 4:7.2e-6, 5:7.0e-6, 6:6.3e-6,\
      7:1.70e-6, 8:2.5e-6, 9:2.2e-6, 10:1.0e-6, 11:1.3e-6}
Psis1 = {1:-0.121, 2:-0.09, 3:-0.218, 4:-0.786, 5:-0.478, 6:-0.299,\
      7:-0.356, 8:-0.63, 9:-0.153, 10:-0.490, 11:-0.405}
theta_w1 = {1:0.07, 2:0.075, 3:0.114, 4:0.1794, 5:0.1547, 6:0.1749,\
      7:0.2181, 8:0.2498, 9:0.2193, 10:0.2838, 11:0.2864}

#Get soil hydraulic value arrays from soil texture
[rows,cols] = Texture.shape
b = np.zeros(Texture.shape)
thetas = np.zeros(Texture.shape)
Ks = np.zeros(Texture.shape)
psis = np.zeros(Texture.shape)
theta_w = np.zeros(Texture.shape)
for row in range(rows):
    for col in range(cols):
        if mask[row,col]:
            b[row, col] = b1[Texture[row, col]]
            thetas[row, col] = thetas1[Texture[row, col]]
            Ks[row, col] = Ks1[Texture[row, col]]*1000*3600*24
            psis[row, col] = Psis1[Texture[row, col]]
            theta_w[row, col] = theta_w1[Texture[row, col]]
        else:
            b[row, col] = np.NaN
            thetas[row, col] = np.NaN
            Ks[row, col] = np.NaN
            psis[row, col] = np.NaN
            theta_w[row, col] = np.NaN
    
#Modeling setting
freq = 1/10   #return frequency between days (1/d)
N = 365  #modeling time span (d)
tday = range(N)
dt = 1
tt = np.arange(0,N,dt)  #time series
Nm = len(tt)

#Listed below are some constant parameters of my model. I didn't annotate the meanings of all of them as they do not matter much for this project.
#Plant properties
LAImax = 4
Zrmax = 600
k_intercept = 0.05
Time2maxLAI = 365
#Climatic condition
Ca = 400/1000000
RH = 0.6
#Contaminant properties
kd = 1e-5
TSCF = 0.5
xo = 1
tox = 0.01
SRL = 120
ma = 0.15

#The array to store phytoremediation efficiency result
PReff = np.zeros(precip.shape)
Rmv = np.zeros(precip.shape)
PReff[mask == False] = np.NaN
Rmv[mask == False] = np.NaN
#Modeling the phytoremediation process
for row in range(precip.shape[0]):
    for column in range(precip.shape[1]):
        if mask[row,column]:
            #Hydrology
            VPD = 0.611*math.exp(17.502*temp[row,column]/(249.91+temp[row,column]))*(1-RH)/101
            WUE = 0.625*Ca*0.25/VPD
            dep = (precip[row,column]*2.54*10/365)/freq
            R = precip_generator(freq,dep,N)
            R_interp = np.interp(tt,tday,R)
            ETo = PET[row,column]
            #Soil properties
            rhob = BD[row,column]*1000000
            Por = 1.01*thetas[row,column]
            sw = 0.1*(theta_w[row,column]/Por)
            s1 = thetas[row,column]/Por
            #Variables to be modeled
            s = [0.5*s1]
            LAI = [1]
            Zr = [50]
            Msh = [ma]
            Mrt = [50/SRL]
            x = [xo]
            UPxt = [0]
            An = []
            Re = []
            Tr = []
            Ev = []
            ET = []
            Is = []
            LQ = []
            xs = []
            UPx = []
            LEx = []
            for t in range(Nm):
                rs = max((s[t]-sw)/(s1-sw),0)
                xs.append(x[t]/(Por*s[t]+rhob*kd))
                #Water balance
                Tr.append(ETo*(0.33*LAI[t]+0.45)*(-2*s[t]**3+3*s[t]**2)-tox*UPxt[t])
                Ev.append(ETo*math.exp(-0.398*LAI[t])*rs)
                ET.append(min((Tr[t]+Ev[t]),ETo))
                An.append(12*Tr[t]*WUE/18)
                a_intercept = math.exp(-k_intercept*LAI[t])
                Is.append(min(a_intercept*R_interp[t], Ks[row,column]))
                LQ.append(Ks[row,column]*rs**(2*b[row,column]+3))
                s.append(max(s[t]+dt*(Is[t]-ET[t]-LQ[t])/(Por*Zr[t]),0.99*sw))
                s[t+1] = min(s[t+1],1)
                #Contaminant balance
                UPx.append(TSCF*ET[t]*xs[t])
                LEx.append(LQ[t]*xs[t])
                x.append(max(x[t]-dt*(UPx[t]+LEx[t])/Zr[t],0))
                UPxt.append(UPxt[t]+dt*UPx[t])
                #Plant growth
                Re.append((Msh[t]+Mrt[t])/Time2maxLAI)
                dMsh = max(0.75*(An[t]-Re[t])*2*dt,0)
                Msh.append(Msh[t]+dMsh)
                LAI.append(min(Msh[t+1]/ma,LAImax))
                dMrt = max(0.25*(An[t]-Re[t])*2*dt,0)
                Mrt.append(Mrt[t]+dMrt)
                Zr.append(min(Mrt[t+1]*SRL,Zrmax))
            arr_LEx = np.array(LEx)  #lists need to be converted to NumPy array for index-matching multiply
            Rmv[row,column] = 1-x[Nm]/xo
            PReff[row,column] = UPxt[Nm]/(sum(arr_LEx*dt)+UPxt[Nm])

#Get the coordinate of lowerleft corner of NC for georeference
descNC = arcpy.Describe('nc_map')
leftbnd = descNC.extent.XMin
bottombnd = descNC.extent.YMin
lowerleft = arcpy.Point(leftbnd,bottombnd)

LC = arcpy.RasterToNumPyArray('LCmask')
EffEval = PReff*LC
RmvEval = Rmv*LC
EffResult = arcpy.NumPyArrayToRaster(EffEval,lowerleft,1000,1000)
RmvResult = arcpy.NumPyArrayToRaster(RmvEval,lowerleft,1000,1000)
EffResult.save('EffEval.img')
RmvResult.save('RmvEval.img')

#%% Stochastic rainfall generation function
def precip_generator(freq,dep,N):
    Pr = np.zeros(N)
    r = np.random.rand(N)
    beta = freq
    tau = (-1/beta)*np.log(1-r+np.spacing(1))
    beta = 1/dep
    eta = (-1/beta)*np.log(1-r+np.spacing(1))
    tt = np.cumsum(tau)-tau[0]+1
    tti = np.floor(tt)
    tti = tti.astype(int)
    P = np.zeros(tti[-1])
    for i in range(N):
        P[tti[i]-1] = eta[i-1]
    Pr = P[:N]
    return Pr