###################################
## fit of the daily radio signal ##
###################################
import numpy as np
import matplotlib.pyplot as plt
import sys
import os
cwd = os.getcwd()
classpath = cwd + '/../classes/'
utilspath = cwd + '/../utils/'
sys.path.append(utilspath)
sys.path.append(classpath)
from scipy.optimize import curve_fit
import datetime
import constant
import utils
import pickle
import dataset
import daydata
import argparse
import matplotlib as mpl
def fitwithexpo0(x, y, yerr, a, sigma,mu,c):
    try:
        popt, pcov = curve_fit(utils.expofunc0,x,y,sigma=yerr, p0=[a,sigma,mu,c])
    except RuntimeError:
        print("Error - curve_fit failed")
    #        popt, pcov = curve_fit(utils.expofunc,x,radio, p0=[a,sigma,mu,b,c])
    return [popt,pcov]




############################
##  argument parser       ##
############################
parser = argparse.ArgumentParser()
parser.add_argument("stname", type=str, nargs='?',default='popey', help="station name")

args = parser.parse_args()
stname = args.stname
stid = constant.GDstationidbyname[stname]
goodlistname = stname + '.txt'
[goodyear,goodday] = utils.getdaysfromlist(constant.listfolder + goodlistname)

temp = np.array([])
radio = np.array([])
nosuntime = np.array([])
nosunradio = np.array([])
nosuntemp = np.array([])
tsysall =  np.array([])
sigmatsysall =  np.array([])
sigmatsysstat =  np.array([])

datearray =[]
tempcorrmean = np.array([])
tempcorrsigma = np.array([])
ahmax = np.array([])
asimhmax = np.array([])
amax = np.array([])
adate = np.array([])
for y,d in zip(goodyear,goodday):
    outfolder = constant.resultfoldertest + '/t0p0/' + str(stid)
    datafilename = outfolder + '/data_' + str(int(y)) + '_' + str(int(d)) + '.pkl'
    if os.path.isfile(datafilename):
        datafile = open(datafilename,'rb')
    else:
        continue
    thedata = pickle.load(datafile)
    #load simulation
    timesim = thedata.timesim
    sim = thedata.sim
    maxsim = np.max(sim)
    #load data
    time = thedata.timedata
     ## period of the day when the sun is present
    timeofsunmax = timesim[np.argmax(sim)]/3600 + 3
    [poptsim,pcovsim] = fitwithexpo0(timesim/3600 +3, sim, 1, 10, 1.5,timeofsunmax,0)
#    print 'poptsim = ' , poptsim[2]
    t1 = utils.doytodate(int(y),int(d),int(timeofsunmax - 2), int(np.modf(timeofsunmax - 2)[0]))
    t2 = utils.doytodate(int(y),int(d),int(timeofsunmax + 2), int(np.modf(timeofsunmax + 2)[0]))
    nsradio = thedata.data[(thedata.timedata < t1)  | (thedata.timedata > t2) ]
    nstemp = thedata.tempLL[(thedata.timedata < t1)  | (thedata.timedata > t2) ]
    ## subtract the mean of the time when no sun is present in the data
    radio = thedata.data - np.mean(nsradio)
    temp = thedata.tempLL - np.mean(nstemp)
    ## get the file with temperature fit result
    fitfile = constant.resultfolder + str(stid) + '/tempfit_'+goodlistname[:-4] + '.npy'
#    fitfile = constant.resultfolder + str(stid) + '/tempfit.npy'
    tempfit = np.load(fitfile)
    pfit = np.poly1d(tempfit)
    radio = radio - pfit(temp)
    ## baseline corrected
    nsradioc = nsradio - np.mean(nsradio)
    nstemp = nstemp - np.mean(nstemp)
    nsradioc = nsradioc -  pfit(nstemp)
    print 'np.mean(nsradio) = ' , np.mean(nsradioc)
    ## uncertainty:
    uncert = np.std(nsradioc)
    tempcorrmean =np.append(tempcorrmean,np.mean(nsradioc))
    tempcorrsigma =np.append(tempcorrsigma,uncert)
    print '!!!!!!!!!!!!!!!uncert !!!!!!!!!!!!!! = ', uncert 


fig = plt.figure()
fig.suptitle(stname,fontsize=15, fontweight='bold')
plt.errorbar(np.arange(0,len(tempcorrmean),1), tempcorrmean, yerr=tempcorrsigma) 
plt.xlabel('meas nr')
plt.ylabel('std(radio corrected) [ADC]')
plt.show()
