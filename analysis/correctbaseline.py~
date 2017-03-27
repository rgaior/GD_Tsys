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
c_radio = np.array([])
nosuntime = np.array([])
nosunradio = np.array([])
nosuntemp = np.array([])
fig = plt.figure(figsize=(12,6))
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
    t1 = utils.doytodate(int(y),int(d),16,0)
    t2 = utils.doytodate(int(y),int(d),19,0)
    radio = thedata.data[ (thedata.timedata < t1)  | (thedata.timedata > t2) ]
    temp = thedata.tempLL[(thedata.timedata < t1)  | (thedata.timedata > t2) ]
    radio = radio-np.mean(radio)
    temp = temp-np.mean(temp)
    fitfile = constant.resultfolder + str(stid) + '/tempfit.npy'
    tempfit = np.load(fitfile)
    pfit = np.poly1d(tempfit)
    c_radio = np.append(c_radio,radio - pfit(temp))
    nosuntime = np.append(nosuntime,time[ (thedata.timedata < t1)  | (thedata.timedata > t2) ])

#plt.plot(nosuntime,c_radio)
bins= np.linspace(-500,500,200)
plt.hist(c_radio,bins=bins,log=True)
    
plt.show()
