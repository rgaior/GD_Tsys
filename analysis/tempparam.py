#######################################
## 
#######################################
import numpy as np
import matplotlib.pyplot as plt
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
nosunradio = np.array([])
nosuntemp = np.array([])
fig = plt.figure(figsize=(8,8))
fig.suptitle(stname ,fontweight='bold',fontsize=15)
startday = 0
endday = 5
datafolder = '/t'+str(delt) + 'p'+str(delp)+'/'
for y,d in zip(goodyear,goodday):
#for y,d in zip(goodyear[startday:endday],goodday[startday:endday]):
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


#     radio = np.append(radio,time )
#     temp = np.append(temp, thedata.tempLL)

[fit,cov] = np.polyfit(nosuntemp,nosunradio,1,cov=True)
pfit = np.poly1d(fit)
print fit, ' ' , np.sqrt(np.diag(cov))
plt.plot(nosuntemp,nosunradio,'b.')
plt.plot(nosuntemp,pfit(nosuntemp),'r',lw=2)
plt.xlabel('temperature at Los Leones [C]')
plt.ylabel('radio baseline [ADC]')
if savefit == True:
    outfilefit = constant.resultfolder + str(stid) + '/tempfit_'+goodlistname[:-4]
    print outfilefit
    np.save(outfilefit,fit)
plt.legend()
plt.show()
