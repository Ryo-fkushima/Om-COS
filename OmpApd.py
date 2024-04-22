# OmpApd.py
# version: 1.0.1 (Apr 22, 2024)
# author: Ryo Fukushima (Tohoku University) rpifukushima@gmail.com 
#
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate
import time
from numba import njit
import math
import matplotlib.animation as animation


rrr = 8.3143  # molar gas constant in J/K/mol

def gibbs(phi, temp): # excess Gibbs energy of omphacite
    return (11.4 * (temp - 1138) * phi**2) + (4317 * phi**6) # Carpenter et al. (1990, EJM)

@njit(cache = True)
def dgibbs(phi, temp): # 1st derivative of the excess Gibbs energy
    return (22.8 * (temp - 1138) * phi) + (25902 * phi**5)

def phi_m(temp): # phi value with the minimum Gibbs energy
    return (22.8 / 25902 * (1138 - temp))**(1/4)

def gibbs_m(temp): # minimum Gibbs energy
    return phi_m(temp)**2 * (7.6*temp - 8648.8) * (-1)

def gibbs_bup(phi, temp):
    return gibbs(phi, temp) + gibbs_m(temp) 

@njit(cache = True)
def PhaseFieldLoop2d(nsteps, ny, nx, dx, dy, p, dt, temp, aaa, periodi, periodj): # main loop (2d). Allen & Cahn (1979)

    if periodi == 1:
        ip_replace = 0
        im_replace = nx - 1
        
    if periodi == 0:
        ip_replace = nx - 1
        im_replace = 0
        
    if periodj == 1:
        jp_replace = 0
        jm_replace = ny - 1
        
    if periodj == 0:
        jp_replace = ny - 1
        jm_replace = 0    
        
    
    for t in range(nsteps):
        for j in range(ny):
            for i in range(nx):
                ip = i + 1
                im = i - 1
                jp = j + 1
                jm = j - 1
                if ip > nx - 1:
                    ip = ip_replace 
                if im < 0:
                    im = im_replace 
                if jp > ny - 1:
                    jp = jp_replace 
                if jm < 0:
                    jm = jm_replace 
                    
                nabx = (p[ip, j, t] - 2 * p[i, j, t] + p[im, j, t])/dx/dx
                naby = (p[i, jp, t] - 2 * p[i, j, t] + p[i, jm, t])/dy/dy
                
                p[i,j,(t+1)] = p[i,j, t] - (dgibbs(p[i,j,t], temp) - (aaa * aaa * (nabx + naby))) * dt / rrr / temp
        
        if (t + 1) % 50 == 0:
            print("nstep = ", t+1)
    
    return p
    

def CalcApd2d(**args): # main function for APD calculation (2d)

    time_start = time.time()
    
    omp_v, nx, ny, dx, dy = args["omp_v"], args["nx"], args["ny"], args["dx"], args["dy"]
    temp, aaa, intE, ini = args["temp"], args["aaa"], args["intE"], args["ini"]
    turb, dt, nsteps = args["turb"], args["dt"], args["nsteps"]
    periodi, periodj, inifield = args["periodi"], args["periodj"], args["inifield"]
    
    if aaa != 0:
        print("given gradient coef [m (J/mol)^0.5] = ", "{:.2g}".format(aaa))
        
    if aaa == 0:
        x = np.linspace(-phi_m(temp) + 0.00001, phi_m(temp) - 0.00001, num = 100, endpoint = True)
        y = (gibbs_bup(x, temp))**0.5
        S = integrate.simps(y,x)
        aaa = (intE * omp_v)/(np.sqrt(2) * S) # Cahn & Hilliard (1958)
        print("APB energy [J/m^2] = ", "{:.2g}".format(intE))
        print("calculated gradient coef [m (J/mol)^0.5] = ", "{:.2g}".format(aaa))
        
    print("temp [K] = ", "{:.3g}".format(temp))

    if ini == -1:
        print("user setting for the initial field")
        
    if ini == 0:
        print("sd of initial Gaussian = ", "{:.3g}".format(turb))

    if ini == 1:
        print("range of initial uniform distribution = ", "{:.3g}".format(turb))    
    
    p = np.zeros((nx, ny, (nsteps + 1)))
    
    if ini == -1:
        p[:,:,0] = inifield
    
    if ini == 0:
        p[:,:,0] = np.random.normal(loc = 0, scale = turb, size = (nx, ny))

    if ini == 1:
        p[:,:,0] = 0 + (np.random.rand(nx, ny) - 0.5) / 0.5 * turb
    
    PF = PhaseFieldLoop2d(nsteps, ny, nx, dx, dy, p, dt, temp, aaa, periodi, periodj)
    
    time_end = time.time()
    TotalTime = time_end - time_start
    print("PF calculation time = ", "{:.4f}".format(TotalTime), "s")
        
    return PF
    

@njit(cache = True)
def PhaseFieldLoop3d(nsteps, nz, ny, nx, dx, dy, dz, p, dt, temp, aaa, periodi, periodj, periodk): # main loop (3d)

    if periodi == 1:
        ip_replace = 0
        im_replace = nx - 1
        
    if periodi == 0:
        ip_replace = nx - 1
        im_replace = 0
        
    if periodj == 1:
        jp_replace = 0
        jm_replace = ny - 1
        
    if periodj == 0:
        jp_replace = ny - 1
        jm_replace = 0   
        
    if periodk == 1:
        kp_replace = 0
        km_replace = nz - 1
        
    if periodk == 0:
        kp_replace = nz - 1
        km_replace = 0       
        
    
    for t in range(nsteps):
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    ip = i + 1
                    im = i - 1
                    jp = j + 1
                    jm = j - 1
                    kp = k + 1
                    km = k - 1
                    if ip > nx - 1:
                        ip = ip_replace 
                    if im < 0:
                        im = im_replace 
                    if jp > ny - 1:
                        jp = jp_replace 
                    if jm < 0:
                        jm = jm_replace
                    if kp > nz - 1:
                        kp = kp_replace 
                    if km < 0:
                        km = km_replace     
                        
                    nabx = (p[ip, j, k, t] - 2 * p[i, j, k, t] + p[im, j, k, t])/dx/dx
                    naby = (p[i, jp, k, t] - 2 * p[i, j, k, t] + p[i, jm, k, t])/dy/dy
                    nabz = (p[i, j, kp, t] - 2 * p[i, j, k, t] + p[i, j, km, t])/dz/dz
                    
                    p[i,j,k,(t+1)] = p[i,j,k,t] - (dgibbs(p[i,j,k,t], temp) - (aaa * aaa * (nabx + naby + nabz))) * dt / rrr / temp
        
        if (t + 1) % 50 == 0:
            print("nstep = ", t+1)
    
    return p
    

def CalcApd3d(**args): # main function for APD calculation (3d)

    time_start = time.time()
    
    omp_v, nx, ny, nz, dx, dy, dz = args["omp_v"], args["nx"], args["ny"], args["nz"], args["dx"], args["dy"], args["dz"]
    temp, aaa, intE, ini = args["temp"], args["aaa"], args["intE"], args["ini"]
    turb, dt, nsteps = args["turb"], args["dt"], args["nsteps"]
    periodi, periodj, periodk, inifield = args["periodi"], args["periodj"], args["periodk"], args["inifield"]
    
    if aaa != 0:
        print("given gradient coef [m (J/mol)^0.5] = ", "{:.2g}".format(aaa))
        
    if aaa == 0:
        x = np.linspace(-phi_m(temp) + 0.00001, phi_m(temp) - 0.00001, num = 100, endpoint = True)
        y = (gibbs_bup(x, temp))**0.5
        S = integrate.simps(y,x)
        aaa = (intE * omp_v)/(np.sqrt(2) * S)
        print("APB energy [J/m^2] = ", "{:.2g}".format(intE))
        print("calculated gradient coef [m (J/mol)^0.5] = ", "{:.2g}".format(aaa))
        
    print("temp [K] = ", "{:.3g}".format(temp))

    if ini == -1:
        print("user setting for the initial field")
        
    if ini == 0:
        print("sd of initial Gaussian = ", "{:.3g}".format(turb))

    if ini == 1:
        print("range of initial uniform distribution = ", "{:.3g}".format(turb))    
    
    p = np.zeros((nx, ny, nz, (nsteps + 1)))
    
    if ini == -1:
        p[:,:,:,0] = inifield
    
    if ini == 0:
        p[:,:,:,0] = np.random.normal(loc = 0, scale = turb, size = (nx, ny, nz))

    if ini == 1:
        p[:,:,:,0] = 0 + (np.random.rand(nx, ny, nz) - 0.5) / 0.5 * turb
    
    PF = PhaseFieldLoop3d(nsteps, nz, ny, nx, dx, dy, dz, p, dt, temp, aaa, periodi, periodj, periodk)
    
    time_end = time.time()
    TotalTime = time_end - time_start
    print("PF calculation time = ", "{:.4f}".format(TotalTime), "s")
        
    return PF


@njit(cache = True)
def PhaseFieldLoop3d_CS(nsteps, nz, ny, nx, dx, dy, dz, p, p_save, zno, dt, temp, aaa, periodi, periodj, periodk): # main loop (3d but save only z = zno cross section)
	
    p_new = p
	
    if periodi == 1:
        ip_replace = 0
        im_replace = nx - 1
        
    if periodi == 0:
        ip_replace = nx - 1
        im_replace = 0
        
    if periodj == 1:
        jp_replace = 0
        jm_replace = ny - 1
        
    if periodj == 0:
        jp_replace = ny - 1
        jm_replace = 0   
        
    if periodk == 1:
        kp_replace = 0
        km_replace = nz - 1
        
    if periodk == 0:
        kp_replace = nz - 1
        km_replace = 0       
        
    
    for t in range(nsteps):
        for k in range(nz):
            for j in range(ny):
                for i in range(nx):
                    ip = i + 1
                    im = i - 1
                    jp = j + 1
                    jm = j - 1
                    kp = k + 1
                    km = k - 1
                    if ip > nx - 1:
                        ip = ip_replace 
                    if im < 0:
                        im = im_replace 
                    if jp > ny - 1:
                        jp = jp_replace 
                    if jm < 0:
                        jm = jm_replace
                    if kp > nz - 1:
                        kp = kp_replace 
                    if km < 0:
                        km = km_replace     
                        
                    nabx = (p[ip, j, k] - 2 * p[i, j, k] + p[im, j, k])/dx/dx
                    naby = (p[i, jp, k] - 2 * p[i, j, k] + p[i, jm, k])/dy/dy
                    nabz = (p[i, j, kp] - 2 * p[i, j, k] + p[i, j, km])/dz/dz
                    
                    p_new[i,j,k] = p[i,j,k] - (dgibbs(p[i,j,k], temp) - (aaa * aaa * (nabx + naby + nabz))) * dt / rrr / temp
                    
        p_save[:,:,(t + 1)] = p_new[:,:,zno]
        
        p = p_new
        
        if (t + 1) % 50 == 0:
            print("nstep = ", t+1)
    
    return p_save
 

def CalcApd3d_CS(**args): # main function for APD calculation (3d but save only a z cross section)

    time_start = time.time()
    
    omp_v, nx, ny, nz, dx, dy, dz = args["omp_v"], args["nx"], args["ny"], args["nz"], args["dx"], args["dy"], args["dz"]
    temp, aaa, intE, ini = args["temp"], args["aaa"], args["intE"], args["ini"]
    turb, dt, nsteps = args["turb"], args["dt"], args["nsteps"]
    periodi, periodj, periodk, inifield = args["periodi"], args["periodj"], args["periodk"], args["inifield"]
    zno = args["slicefor3d_cs"]
    
    if aaa != 0:
        print("given gradient coef [m (J/mol)^0.5] = ", "{:.2g}".format(aaa))
        
    if aaa == 0:
        x = np.linspace(-phi_m(temp) + 0.00001, phi_m(temp) - 0.00001, num = 100, endpoint = True)
        y = (gibbs_bup(x, temp))**0.5
        S = integrate.simps(y,x)
        aaa = (intE * omp_v)/(np.sqrt(2) * S)
        print("APB energy [J/m^2] = ", "{:.2g}".format(intE))
        print("calculated gradient coef [m (J/mol)^0.5] = ", "{:.2g}".format(aaa))
        
    print("temp [K] = ", "{:.3g}".format(temp))

    if ini == -1:
        print("user setting for the initial field")
        
    if ini == 0:
        print("sd of initial Gaussian = ", "{:.3g}".format(turb))

    if ini == 1:
        print("range of initial uniform distribution = ", "{:.3g}".format(turb))    
    
    p = np.zeros((nx, ny, nz))
    p_save = np.zeros((nx, ny, (nsteps + 1)))
    
    if ini == -1:
        p[:,:,:] = inifield
    
    if ini == 0:
        p[:,:,:] = np.random.normal(loc = 0, scale = turb, size = (nx, ny, nz))

    if ini == 1:
        p[:,:,:] = 0 + (np.random.rand(nx, ny, nz) - 0.5) / 0.5 * turb
        
    p_save[:,:,0] = p[:,:,zno]    
    
    PF = PhaseFieldLoop3d_CS(nsteps, nz, ny, nx, dx, dy, dz, p, p_save, zno, dt, temp, aaa, periodi, periodj, periodk)
    
    time_end = time.time()
    TotalTime = time_end - time_start
    print("PF calculation time = ", "{:.4f}".format(TotalTime), "s")
        
    return PF




def MobCalc2d(PF, **args): # calculate mobility and real time
    
    mob, premob, actEmob, temp, dt = args["mob"], args["premob"], args["actEmob"], args["temp"], args["dt"]
    
    if mob != 0:
        print("fixed mobility [mol/J/yr] = ", "{:.2g}".format(mob))

    if mob == 0:
        mob = premob * np.exp(-actEmob / rrr / temp)
        print("calculated mobility [mol/J/yr]", "{:.2g}".format(premob),"exp(-", "{:.2g}".format(actEmob),"/RT) = ", "{:.2g}".format(mob))
    
    TimeAst = np.linspace(0, (PF.shape[2] - 1), num = PF.shape[2]) * dt
    TimeReal = TimeAst / rrr / temp / mob
    
    return TimeReal
    

def MeanCalc2d(PF): # calculate mean order parameter
    
    MeanOrder = np.zeros(PF.shape[2])
    AbsPF = np.absolute(PF)
    
    for i in range(PF.shape[2]):
        
        MeanOrder[i] = np.mean(AbsPF[:,:,i])
        
    return MeanOrder    


def HistPlot2d(PF, select): # generate a histogram of order parameters
    
    AbsPF = np.absolute(PF)
    AbsPFc = AbsPF[:,:,select]
    AbsPFc1 = np.ravel(AbsPFc)
    
    plt.figure(figsize = (8,5))

    plt.xlabel("Order Parameter")
    
    plt.ylabel("Density")
    
    plt.title("Time step = %d" %select)
    
    plt.hist(AbsPFc1, bins=10, range=(0, 1), density = True)

    plt.show()
    
    return
    

@njit(cache = True)
def PFAbsRavel2d_jit(PF): # loop for HistAnime2d
	
	AbsPF_allTime = np.zeros(((PF.shape[0] * PF.shape[1]), PF.shape[2]))

	AbsPF = np.absolute(PF)
	
	for i in range(PF.shape[2]):
	
		AbsPF_allTime[:,i] = np.ravel(AbsPF[:,:,i])
	
	return AbsPF_allTime
	
    
def HistAnime2d(PF, BinsInput, speed, filename, **args): # make a gif image of histograms for calculated phase-fields

	plotNo = args["plotNo"]
	
	AbsPF_allTime = PFAbsRavel2d_jit(PF)
	
	fig = plt.figure()
	plt.ylim(0,AbsPF_allTime.shape[0])
	plt.tick_params(labelbottom=False, labelleft=False, labelright=False, labeltop=False, left = False, bottom = False)
	ax = fig.add_subplot(111)
	
	def PlotHistFrames(i):
	
		ax.cla()
		ax.title.set_text("Time step = %d" %i)
		ax.set_xlim(0,1)
		ax.set_ylim(0,AbsPF_allTime.shape[0])
		ax.set_xlabel("Order parameter")
		ax.set_ylabel("Fraction, %")
		ax.xaxis.set_major_locator(mpl.ticker.LinearLocator(11))
		ax.yaxis.set_major_formatter(mpl.ticker.PercentFormatter(AbsPF_allTime.shape[0],0,""))
		ax.yaxis.set_major_locator(mpl.ticker.LinearLocator(11))
		ax.hist(AbsPF_allTime[:,i], bins = BinsInput, range = (0,1), density = False, color = "dimgray")
		
	anim = animation.FuncAnimation(fig, PlotHistFrames, frames = plotNo, interval = speed)
	
	anim.save("%s.gif" %filename, writer = "pillow")
	plt.close()

	return	
	

@njit(cache = True)
def MakeDfMatrix2d_jit(PF2d, saturation, x1, yfactor, A): # loop for MakeDfMatrix2d
    
    #x1 = 0.1
    
    #yfactor = 0.0001
    
    #A = 6
    
    y1 = yfactor * A * (1 - x1) * (1 - x1) / (1 - yfactor)
    
    DfMatrix2d = np.zeros((PF2d.shape[0], PF2d.shape[1]))
    
    for i in range(PF2d.shape[0]):
        for j in range(PF2d.shape[1]):
                
            if np.absolute(PF2d[i,j]) > x1:
                DfMatrix2d[i,j] = A * (np.absolute(PF2d[i,j]) - x1)**2 + y1
                
            else:
                DfMatrix2d[i,j] = y1 - (y1 * ((np.absolute(PF2d[i,j]) - x1) / x1)**2)
                
    for i in range(PF2d.shape[0]):
        for j in range(PF2d.shape[1]):           
                
            if DfMatrix2d[i,j] > (saturation * (A * (1 - x1) * (1 - x1) + y1)):
                DfMatrix2d[i,j] = (saturation * (A * (1 - x1) * (1 - x1) + y1))
                 
            DfMatrix2d[i,j] = DfMatrix2d[i,j] / saturation / ((A * (1 - x1) * (1 - x1) + y1)) 
                
                    
    return DfMatrix2d     
    

    
def MakeDfMatrix2d(PF2d, **args): # convert a 2d phase-field into a gray-scale image
    
    saturation = args["saturation"]
    x1, yfactor, A = args["phiadjust"], args["Ifactor"], args["Icoef"]
    
    return MakeDfMatrix2d_jit(PF2d, saturation, x1, yfactor, A)    



@njit(cache = True)
def MakeBinMatrix2d_jit(PF2d, binthres): # loop for MakeBinMatrix2d
    
    BinMatrix2d = np.zeros((PF2d.shape[0], PF2d.shape[1]))
    
    for i in range(PF2d.shape[0]):
        for j in range(PF2d.shape[1]):
                
            if np.absolute(PF2d[i,j]) >  binthres:
                BinMatrix2d[i,j] = 1
                    
            else:
                BinMatrix2d[i,j] = 0
                    
    return BinMatrix2d    
    

def MakeBinMatrix2d(PF2d, **args): # convert a 2d phase-field into a binarized image
    
    binthres = args["binthres"]
    
    return MakeBinMatrix2d_jit(PF2d, binthres)         
    
    
def PlotApd2d(PF, plotmode, **args): # plot calculated phase-fields (basically not used)
    
    plotNo = args["plotNo"]
    
    
    if plotmode == -1:
        
        for i in plotNo:
            
            plt.figure()
            plt.title("Time step = %d" %i)
            plt.imshow(PF[:,:,i], cmap = "twilight")
            plt.colorbar(label = "Order parameter Q")
            plt.show()
    
    if plotmode == 0:
        
        for i in plotNo:
            
            plt.figure()
            plt.title("Time step = %d" %i)
            plt.imshow(PF[:,:,i], cmap = "twilight",vmin = -1, vmax = 1)
            plt.colorbar(label = "Order parameter Q")
            plt.show()
        
    
    if plotmode == 1:
        
        for i in plotNo:
           
           plt.figure()
           plt.title("Time step = %d" %i)
           plt.imshow(MakeDfMatrix2d(PF[:,:,i], **args), cmap = "gray", vmin = 0, vmax = 1) 
           plt.show()
           
        
    if plotmode == 2:
        
        for i in plotNo:
           
           plt.figure()
           plt.title("Time step = %d" %i)
           plt.imshow(MakeBinMatrix2d(PF[:,:,i], **args), cmap = "gray", vmin = 0, vmax = 1) 
           plt.show()
        
    return




def ApdAnime2d(PF, speed, plotmode, filename, **args): # make a gif image for calculated phase-fields
    
    plotNo = args["plotNo"]
    
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    
    if plotmode == -1:
        
        def PlotApdFrames(i):
        
            ax.cla()
            ax.title.set_text("Time step = %d" %i)
            ax.imshow(PF[:,:,i], cmap = "twilight")
    
        
        anim = animation.FuncAnimation(fig, PlotApdFrames, frames = plotNo, interval = speed)
            
    
    if plotmode == 0:
        
        fig.colorbar(ax.imshow(PF[:,:,0], cmap = "twilight", vmin = -1, vmax = 1))
        
        def PlotApdFrames(i):
        
            ax.cla()
            ax.title.set_text("Time step = %d" %i)
            ax.imshow(PF[:,:,i], cmap = "twilight", vmin = -1, vmax = 1)
    
        
        anim = animation.FuncAnimation(fig, PlotApdFrames, frames = plotNo, interval = speed)
            
        
    
    if plotmode == 1:
        
        def PlotApdFrames(i):
        
            ax.cla()
            ax.title.set_text("Time step = %d" %i)
            ax.imshow(MakeDfMatrix2d(PF[:,:,i], **args), cmap = "gray", vmin = 0, vmax = 1)
    
        
        anim = animation.FuncAnimation(fig, PlotApdFrames, frames = plotNo, interval = speed)
          
           
        
    if plotmode == 2:
        
       def PlotApdFrames(i):
       
           ax.cla()
           ax.title.set_text("Time step = %d" %i)
           ax.imshow(MakeBinMatrix2d(PF[:,:,i], **args), cmap = "gray", vmin = 0, vmax = 1)
   
       
       anim = animation.FuncAnimation(fig, PlotApdFrames, frames = plotNo, interval = speed)
           
        
    anim.save("%s.gif" %filename, writer = "pillow")    
    plt.close()    
        
    return







def FracCalc2d(PF, **args): # calculate fraction of P2/n phase in a 2d phase-field
    
    plotNo = args["plotNo"]
    
    Fraction = np.zeros(PF.shape[2])
    
    for i in plotNo:
        
        Fraction[i] = np.mean(MakeBinMatrix2d(PF[:,:,i], **args))
        
    return Fraction    



@njit(cache = True)
def LengthCount1d(PF1d): # calculate mean APD size along a line
    
    domainNo = 0
    
    PF1d_diff = np.zeros(len(PF1d))
    
    for i in range(len(PF1d)):
        
        if i == 0:
            
            im = len(PF1d) - 1 # periodic
            
        else:
            
            im = i - 1
        
        PF1d_diff[i] = PF1d[i] - PF1d[im] 
        
        if PF1d_diff[i]  == 1:
            
            domainNo = domainNo + 1
        
        
    LengthList = np.zeros(domainNo + 1)
    
    count = 0
    
    for i in range(len(PF1d)):
        
        if PF1d[i] == 1:
            
            LengthList[count] = LengthList[count] + 1
            
        if PF1d[i] == 0:
            
            if PF1d_diff[i] == -1:
                
                count = count + 1
                
    if LengthList[domainNo] == 0:
         
         LengthList = np.delete(LengthList, domainNo)
         
    else:
         
         LengthList[0] = LengthList[0] + LengthList[domainNo]
         
         LengthList = np.delete(LengthList, domainNo)


    return LengthList

    

def ApdMeasure2d(PF, **summary): # calculate mean APD size etc in a 2d phase-field
    
    plotNo, ApdMeasureWidth = summary["Visual"]["plotNo"], summary["Visual"]["ApdMeasureWidth"]
    delta = summary["Calc"]["dx"]
    
    MeasuredLineIndex = np.arange(0, PF.shape[0], ApdMeasureWidth)
    
    #ApdSizeList = np.empty(0)
    
    MeanAPD = np.zeros(PF.shape[2])
    HowManyAPD = np.zeros(PF.shape[2])
    sdAPDsize = np.zeros(PF.shape[2])
    
    for k in plotNo:
        
        ApdSizeList = np.empty(0)
        
        sliced_binary = MakeBinMatrix2d(PF[:,:,k], **summary["Visual"])
        
        for i in MeasuredLineIndex:
            
           ApdSizeList = np.append(ApdSizeList, LengthCount1d(sliced_binary[i,:])) 
        
        if len(ApdSizeList) == 0:
            
            MeanAPD[k] = 0
            sdAPDsize[k] = 0
        
        elif len(ApdSizeList) == 1:
            
            MeanAPD[k] = np.mean(ApdSizeList) * delta
            sdAPDsize[k] = 0
            
        else:
                
            MeanAPD[k] = np.mean(ApdSizeList) * delta
            sdAPDsize[k] = np.std(ApdSizeList, ddof = 1) * delta
        
        HowManyAPD[k] = len(ApdSizeList)
        
       
    return MeanAPD, HowManyAPD, sdAPDsize



      
def mean2d(x):
    
    y = np.sum(x) / np.size(x)
    
    return y



      
def corr2d(a,b):
    
    a = a - mean2d(a)
    b = b - mean2d(b)
    
    r = np.sum(a*b) / (np.sum(a*a) * np.sum(b*b))**(0.5)
    
    return r
 


def AcCalcX(PF, MaxOffset, plotNo, delta): # loop for AcCalcX_app. cf. Buscombe et al. (2010)
    
    PF_original = PF
    ACcurve = np.zeros(MaxOffset + 1)
    WaveLength = np.zeros(PF.shape[2])
    xxx = range(MaxOffset + 1)
    
    for k in plotNo:
        
        PF_layer = PF_original[:,0:(-(MaxOffset + 1)),k]
        ACcurve[0] = 1
        
        for i in range(MaxOffset):
            
            PF_compared = PF_original[:,i:(-(MaxOffset + 1) + i),k]
            ACcurve[i + 1] = corr2d(PF_layer, PF_compared)
        
        Fitted = interpolate.interp1d(ACcurve, xxx)   
        
        WaveLength[k] = Fitted(0.5) * delta * 2 * math.pi
        
            
    return WaveLength

    

    
def AcCalcXPlot(PF, MaxOffset, select, absswitch): # plot autocorrelation diagrams
    
    if absswitch == 0:
        PF_original = PF
        
    else:
        PF_original = np.absolute(PF)
        
    PF_layer = PF_original[:,0:(-(MaxOffset + 1)),select]
    ACcurve = np.zeros(MaxOffset + 1)
    ACcurve[0] = 1
    
    for i in range(MaxOffset):
        
        PF_compared = PF_original[:,i:(-(MaxOffset + 1) + i),select]
        ACcurve[i + 1] = corr2d(PF_layer, PF_compared)
        
    plt.figure()

    plt.xlabel("Lag, pix")
    
    plt.ylabel("Autocorrelation")
    
    plt.title("Time step = %d" %select)
    
    plt.plot(ACcurve, linestyle = None)

    plt.show()
    
    return
       


def AcCalcX_app(PF, **summary): # calculate autocorrelation and wavelength
    
    plotNo, MaxOffset = summary["Visual"]["plotNo"], summary["Visual"]["MaxOffset"]
    delta = summary["Calc"]["dx"] 
    
    WaveLength_original = AcCalcX(PF, MaxOffset, plotNo, delta)
    
    
    return WaveLength_original
    
    


def QuantEvalApd2d(PF, Name, WLswitch, **summary): # summary calculation

    time_start = time.time()
    
    ApdMeasureWidth, MaxOffset = summary["Visual"]["ApdMeasureWidth"], summary["Visual"]["MaxOffset"], 
    binthres = summary["Visual"]["binthres"]
    
    QuantData = pd.DataFrame(data = MobCalc2d(PF, **summary["Calc"]), columns = ["%s_time" %Name])
    
    print("time calculated")
    
    QuantData["%s_MeanOrder" % Name] = MeanCalc2d(PF)
    
    print("averaged degree of order calculated")
    
    QuantData["%s_pFraction" % Name] = FracCalc2d(PF, **summary["Visual"])
    
    print("Fraction of P-phase calculated with a threshold of phi =", "{:.3g}".format(binthres) )
    
    QuantData["%s_MeanApdSize" % Name], QuantData["%s_ApdNo" % Name], QuantData["%s_ApdSizeSd" % Name] = ApdMeasure2d(PF, **summary)
    
    print("Mean APD size calculated (interval: %d grids)" % ApdMeasureWidth)
    
    if WLswitch == 1:
    
    	QuantData["%s_WaveLength_original" % Name] = AcCalcX_app(PF, **summary)
    
    	print("Wavelength calculated with an offset of %d grids" % MaxOffset)
    
    time_end = time.time()
    TotalTime = time_end - time_start
    print("Summary calculation time = ", "{:.4f}".format(TotalTime), "s")
    
    return QuantData
    



