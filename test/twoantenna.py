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
parser.add_argument("stnameref", type=str, nargs='?',default='luis', help="station name")
args = parser.parse_args()

stname = args.stname
stnameref = args.stnameref
stid = constant.GDstationidbyname[stname]
stidref = constant.GDstationidbyname[stnameref]
goodlistname = stname + '.txt'
#goodlistname = 'march2015popey.txt'
#goodlistname = 'janfevmarch2016popey.txt'
#goodlistname = 'novdec2015popey.txt'
[goodyear,goodday] = utils.getdaysfromlist(constant.listfolder + goodlistname)
#print goodday
temp = np.array([])
radio = np.array([])
radioref = np.array([])
nosuntime = np.array([])
nosunradio = np.array([])
nosunradioref = np.array([])
nosuntemp = np.array([])
fig = plt.figure(figsize=(8,8))
fig.suptitle(stname ,fontweight='bold',fontsize=15)
startday = 0
endday = 5
for y,d in zip(goodyear,goodday):
#for y,d in zip(goodyear[startday:endday],goodday[startday:endday]):
    print 'y = ' , y,  ' d = ' , d
    outfolder = constant.resultfolder + str(stid)
    outfolderref = constant.resultfolder + str(stidref)
    datafilename = outfolder + '/data_' + str(int(y)) + '_' + str(int(d)) + '.pkl'
    datafilenameref = outfolderref + '/data_' + str(int(y)) + '_' + str(int(d)) + '.pkl'
    if os.path.isfile(datafilename) == False or os.path.isfile(datafilenameref) == False:
        continue
    else:
        print 'aqui !'
        datafile = open(datafilename,'rb')
        datafileref = open(datafilenameref,'rb')
        thedata = pickle.load(datafile)
        thedataref = pickle.load(datafileref)
        time = thedata.timedata
        timeofsunmax = constant.exptime[stid]
        t1 = utils.doytodate(int(y),int(d),int(timeofsunmax - 2), int(np.modf(timeofsunmax - 2)[0]))
        t2 = utils.doytodate(int(y),int(d),int(timeofsunmax + 2), int(np.modf(timeofsunmax + 2)[0]))
        radio = thedata.data[ (thedata.timedata < t1)  | (thedata.timedata > t2) ]
        radior = thedataref.data[ (thedata.timedata < t1)  | (thedata.timedata > t2) ]
        temp = thedata.tempLL[(thedata.timedata < t1)  | (thedata.timedata > t2) ]
        radio = radio-np.mean(radio)
        radior = radior-np.mean(radior)
        temp = temp-np.mean(temp)
        nosuntime = np.append(nosuntime,time[ (thedata.timedata < t1)  | (thedata.timedata > t2) ])
        nosunradio = np.append(nosunradio,radio)
        nosunradioref = np.append(nosunradioref,radior)
        nosuntemp = np.append(nosuntemp,temp)



#     radio = np.append(radio,time )
#     temp = np.append(temp, thedata.tempLL)

plt.plot(nosunradioref,nosunradio)
# fit = np.polyfit(nosuntemp,nosunradio,1)
# pfit = np.poly1d(fit)
# print fit
# plt.plot(nosuntemp,nosunradio,'b.')
# plt.plot(nosuntemp,pfit(nosuntemp),'r',lw=2)
# plt.xlabel('temperature at Los Leones [C]')
# plt.ylabel('radio baseline [ADC]')
# if savefit == True:
#     outfilefit = constant.resultfolder + str(stid) + '/tempfit_'+goodlistname[:-4]
#     np.save(outfilefit,fit)
# plt.legend()
plt.show()
