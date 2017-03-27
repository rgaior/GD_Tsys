#######################################
## 
#######################################
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import sys
import os
cwd = os.getcwd()
classpath = cwd + '/../classes/'
utilspath = cwd + '/../utils/'
sys.path.append(utilspath)
sys.path.append(classpath)
import datetime
import constant
import utils
import pickle
import dataset
import daydata
import argparse

def func(x, a, b):
    return a*x  + b

############################
##  argument parser       ##
############################
parser = argparse.ArgumentParser()
parser.add_argument("stname", type=str, nargs='?',default='popey', help="station name")
parser.add_argument("-save", action='store_true', help="save the fit")
parser.add_argument("--deltatheta", type=int, nargs='?',default=0, help="deviation from central value")
parser.add_argument("--deltaphi", type=int, nargs='?',default=0, help="deviation from central value")
args = parser.parse_args()

stname = args.stname
savefit = args.save
stid = constant.GDstationidbyname[stname]
delt = args.deltatheta
delp = args.deltaphi
goodlistname = stname + '.txt'
#goodlistname = 'march2015popey.txt'
#goodlistname = 'janfevmarch2016popey.txt'
#goodlistname = 'novdec2015popey.txt'
[goodyear,goodday] = utils.getdaysfromlist(constant.listfolder + goodlistname)
#print goodday
temp = np.array([])
radio = np.array([])
nosuntime = np.array([])
fig = plt.figure(figsize=(8,8))
fig.suptitle(stname ,fontweight='bold',fontsize=15)
startday = 0
endday = 5
dayslope = np.array([])
dayerrslope = np.array([])
datearray = []
datafolder = '/t'+str(delt) + 'p'+str(delp)+'/'

for y,d in zip(goodyear,goodday):
#for y,d in zip(goodyear[startday:endday],goodday[startday:endday]):
    nosunradio = np.array([])
    nosuntemp = np.array([])
    print 'y = ' , y,  ' d = ' , d
    outfolder = constant.resultfoldertest + datafolder + str(stid)
    datafilename = outfolder + '/data_' + str(int(y)) + '_' + str(int(d)) + '.pkl'
    if os.path.isfile(datafilename):
        datafile = open(datafilename,'rb')
    else:
        continue
    thedata = pickle.load(datafile)
    time = thedata.timedata
    timeofsunmax = constant.exptime[stid]
    t1 = utils.doytodate(int(y),int(d),int(timeofsunmax - 2), int(np.modf(timeofsunmax - 2)[0]))
    t2 = utils.doytodate(int(y),int(d),int(timeofsunmax + 2), int(np.modf(timeofsunmax + 2)[0]))
    radio = thedata.data[ (thedata.timedata < t1)  | (thedata.timedata > t2) ]
    temp = thedata.tempLL[(thedata.timedata < t1)  | (thedata.timedata > t2) ]
    radio = radio-np.mean(radio)
    temp = temp-np.mean(temp)
    nosuntime = np.append(nosuntime,time[ (thedata.timedata < t1)  | (thedata.timedata > t2) ])
    nosunradio = np.append(nosunradio,radio)
    nosuntemp = np.append(nosuntemp,temp)
####################################################
    try:
        popt, pcov = curve_fit(func,nosuntemp,nosunradio,sigma=10)
    except RuntimeError:
        print("Error - curve_fit failed")
####################################################
#    popt, pcov = np.polyfit(nosuntemp,nosunradio,1,full=True,cov=True)
#    pfit = np.poly1d(fit)
    errslope = np.sqrt(np.diag(pcov))
    dayslope = np.append(dayslope,popt[0])
    dayerrslope = np.append(dayerrslope,errslope[0])
    datestring = t1.strftime('%d %b %Y')
    datearray.append(datestring)
    fig = plt.figure()
    fig.suptitle(datestring)
    plt.plot(nosuntemp,nosunradio,'.')
    plt.plot(nosuntemp,popt[0]*nosuntemp + popt[1])

#     radio = np.append(radio,time )
#     temp = np.append(temp, thedata.tempLL)

figres = plt.figure(figsize=(12,5))
plt.errorbar(np.arange(0, len(dayslope), 1), dayslope,yerr=dayerrslope)
plt.xlabel('meas. nr')
plt.ylabel('slope')
plt.show()
