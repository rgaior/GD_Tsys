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
args = parser.parse_args()

stname = args.stname
stid = constant.GDstationidbyname[stname]
goodlistname = stname + '.txt'
[goodyear,goodday] = utils.getdaysfromlist(constant.listfolder + goodlistname)
print goodday
temp = np.array([])
radio = np.array([])
nosuntime = np.array([])
nosunradio = np.array([])
nosuntemp = np.array([])
fig = plt.figure(figsize=(8,8))
fig.suptitle(stname ,fontweight='bold',fontsize=15)
for y,d in zip(goodyear,goodday):
    outfolder = constant.resultfolder + str(stid)
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

fit = np.polyfit(nosuntemp,nosunradio,1)
pfit = np.poly1d(fit)
print fit
plt.plot(nosuntime,nosunradio,'b.')
plt.ylabel('radio baseline [ADC]')
plt.show()
