import numpy as np
import math
from scipy import optimize as opt
from matplotlib import pyplot as plt
import os



def SSFun(eps, E_mod):
    return eps*E_mod

def SSxxFun(eps, V_mod):
    return -eps*V_mod

def FitCurve(eps, sig, fun, endval=0): #FC = fittedCurve, Cov = Covariance
    FC, cov=opt.curve_fit(fun, eps, sig)
    return FC


f = open('log.nanowire')                      #Opens the file
f1 = open('mean.dat')
fw = open('data.dat', 'w')
fw1 = open('Pois.dat','w')
print(os.path.realpath(f.name))
varBool = False                             #Just a bool variable
varBool1 = False                             #Just a bool variable
varBool2 = False                             #Just a bool variable
firstCase = True                            #Another bool variable
lx = 0                                      #Declaring a variable
pot = 0                                     #Declaring a variable
lxI = 0                                     #Declaring a variable
sigma = 3.345
a0scale = 1.24
scaleToeV = 2.61144742e22/(6.022e23)

ittEqui = []
tempEqui = []
EnergyEqui = []


EpsList = []
SigList = []
EpsXZList = []
elist = []
ediff = [0]
EMod = []
stepL = []

ic  = 0
iFrac = 0
fracture = False
temp1 = input('What is the temp: ')
rate1 = input('What is the rate: ')
BoxSize1 = input('What is the boxRate: ')

steps_frame = 1000

EpsXZListVMD = []
stepVMD = []

for line in f1:
    step = int(line.split()[0])*steps_frame
    XZ = float(line.split()[1])
    if step <10000: continue
    if firstCase:
        XZ0 = XZ
        firstCase = False
    stepVMD.append(step)
    EpsXZListVMD.append((XZ-XZ0)/XZ0)
firstCase = True
for line in f:                              #Loop through the file
    if (varBool):                           #True if prev line was 'Lx...'
        step = int(line.split()[0])        #Typecast step# to int
        temp = float(line.split()[1])        #Typecast temp to float
        toteng = float(line.split()[2])       #Typecast total Energy to float
        pressure = float(line.split()[3])       #Typecast Pressure to float
        volume = float(line.split()[4])       #Typecast volume to float
        f_pyy = float(line.split()[5])       #Typecast pressure in y to float
        Ly = float(line.split()[6])       #Typecast box size to float
        XZ = math.sqrt(volume/Ly)
#        if (potT < pot):                    #If lower than current lowest
#            lx = lxT
#            pot = potT
        if firstCase:                       #If first, finds the initial Lx
            firstCase = False
            L0 = Ly
            P0 = f_pyy
            XZ0 = XZ
        stepL.append(step)
        EpsList.append((Ly-L0)/L0)
        SigList.append(P0*1e-1-f_pyy*1e-1) #MPa
        elist.append(toteng)
        EpsXZList.append((XZ-XZ0)/XZ0)
        if len(elist)>20 and not fracture:
            if(SigList[-10]>SigList[-1]):
 #           if (np.absolute(elist[-2]-elist[-1])>5):
                fracture = True
                iFrac = ic-10
        if (ic+1)%5==0 and not fracture:
            EMod.append((sum(SigList[-6:])/sum(EpsList[-6:]))*1e-3)
        if (step == 110000):
            varBool = False
            break
        ic = ic+1
    if varBool2:
        ittEqui.append(int(line.split()[0]))
        tempEqui.append(float(line.split()[1]))
        EnergyEqui.append(float(line.split()[2]))
        if ittEqui[-1] == 10000:
            varBool2 = False
    if (varBool1 and line == 'Step Temp TotEng Press Volume f_pyy Ly \n'):            #Checks line
        varBool = True
    if (line == 'Step Temp TotEng Press Volume f_pyy Ly \n'):            #Checks line
        varBool1 = True
        varBool2 = True
iFrac2 = [i for i in range(len(stepVMD)) if stepVMD[i] < stepL[iFrac]][-1]
EpsList2 = [EpsList[i] for i in range(len(EpsList)) if stepL[i]%1000==0]
Emod2 = FitCurve(EpsList[1:iFrac],SigList[1:iFrac], SSFun) #Young's modulus
Vmod = FitCurve(EpsList2[1:iFrac2],EpsXZListVMD[1:iFrac2], SSxxFun) #Poisson's Ratio
for i,j in zip(EpsList, SigList):
    fw.write(str(i)+" "+str(j)+"\n")
fw1.write(str(Vmod[0]))
exit()
