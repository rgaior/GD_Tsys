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

def fitwithexpo2(x, y, yerr, a, sigma,mu,b,c,d):
    try:
        popt, pcov = curve_fit(utils.expofunc2,x,y,sigma=yerr, p0=[a,sigma,mu,b,c,d])
    except RuntimeError:
        print("Error - curve_fit failed")
        return [0.0,0.0]
    #        popt, pcov = curve_fit(utils.expofunc,x,radio, p0=[a,sigma,mu,b,c])
    return [popt,pcov]

def gethourarray(datetimearray):
    hourarray = np.array([])
    for d in datetimearray:
        hour = utils.timetohour(d.hour,d.minute)
        hourarray = np.append(hourarray,hour)
    return hourarray

def isgoodfit(stid,popt,tofmax):
    good = True
    # check time of max:
    deltamax = 1
    if popt[2] > tofmax + deltamax or  popt[2] < tofmax - deltamax:
        print 'popt[2] ' , popt[2]
        good = False
    if np.abs(popt[1]) > 1.5 or np.abs(popt[1]) < 0.3:
        print 'popt[1] ',  popt[1]
        good = False
    return good


############################
##  argument parser       ##
############################
parser = argparse.ArgumentParser()
parser.add_argument("stname", type=str, nargs='?',default='popey', help="station name")

args = parser.parse_args()
stname = args.stname

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

aeffuncertname = constant.aeffuncertfolder + '/aeffuncert_' + stname + '.npz'
aeffuncert = np.load(aeffuncertname)
#print aeffuncert['arr_0']
datearray =[]
ahmax = np.array([])
asimhmax = np.array([])
amax = np.array([])
adate = np.array([])
datafolder = '/t0p0/'
outfolder = constant.resultfoldertest + datafolder
#for y,d in zip(goodyear,goodday):
for y,d in zip(goodyear[:30],goodday[:30]):
    outf = outfolder + str(stid)
    datafilename = outf + '/data_' + str(int(y)) + '_' + str(int(d)) + '.pkl'
    print datafilename
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
    nsradioc = nsradioc -  pfit(nstemp)
    
    ## uncertainty:
    uncert = np.std(nsradioc)
    print '!!!!!!!!!!!!!!!uncert !!!!!!!!!!!!!! = ', uncert 
#    relaeffuncert = 2*aeffuncert['arr_1'][int(d)]
    relaeffuncert =0
    hourarray = gethourarray(time)
    datestring = t1.strftime('%d %b %Y')
    fig = plt.figure()
    fig.suptitle(stname + ' (' + t1.strftime('%d %b %Y')+')',fontweight='bold',fontsize=15)
    plt.plot(hourarray,radio,'.') 
    zero = np.zeros(len(hourarray))
    plt.plot(np.linspace(np.min(hourarray),np.max(hourarray),len(hourarray)),zero,lw=2)
    plt.xlabel('time [UTC]')
    plt.ylabel('baseline [adc]')
plt.show()
