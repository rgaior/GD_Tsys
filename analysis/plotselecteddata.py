#######################################
## make a selection of the day       ##
## usable for the analysis based on  ##
## the humidity data                 ##
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
parser.add_argument("daylist", type=str, nargs='?',default='march2015.txt', help="text file with list of year doy")
parser.add_argument("stname", type=str, nargs='?',default='popey', help="station name")
parser.add_argument("fluxorigin", type=str, nargs='?',default='nobeyama', help="nobeyama or canadian")
args = parser.parse_args()
listname = args.daylist
stname = args.stname
flux = args.fluxorigin
stid = constant.GDstationidbyname[stname]
[year,day] = utils.getdaysfromlist(constant.listfolder + listname)
goodlistname = stname + '.txt'
[goodyear,goodday] = utils.getdaysfromlist(constant.listfolder + goodlistname)

if flux == 'nobeyama':
    outfolder = constant.resultfolder2
elif flux == 'canadian':
    outfolder = constant.resultfolder
else:
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print 'choose a possible type of flux'
    print 'now exiting the script'
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    sys.exit()    

print goodday
humidity = np.array([])
time = np.array([])
radio = np.array([])
d_array = np.array([])
h_array = np.array([])
#for y,d in zip(year[:-2],day[:-2]):
fig = plt.figure(figsize=(12,6))
#fig.suptitle(stname + ' ' + constant.periods[listname],fontweight='bold',fontsize=15)
sel_time = np.array([])
nonsel_time = np.array([])
sel_data = np.array([])
nonsel_data = np.array([])
selcount = 0
nonselcount = 0
for y,d in zip(year,day):
    outf = outfolder + str(stid)
    datafilename = outf + '/data_' + str(int(y)) + '_' + str(int(d)) + '.pkl'
    if os.path.isfile(datafilename):
        datafile = open(datafilename,'rb')
    else:
        continue

#datafile = open(datafilename, 'rb')
    thedata = pickle.load(datafile)
#    for t in thedata.timedata:
#        time = np.append(time,10000*d + t.hour*100 + t.minute)
    time = np.append(time,thedata.timedata)
    humidity = np.append(humidity, thedata.humidity)
    radio = np.append(radio, thedata.data)
    d_array = np.append(d_array,d)
    h_array = np.append(h_array,np.mean(thedata.humidity))
    
    datafile.close()
    print y, ' ' ,  goodyear
    if (y,d) in zip(goodyear,goodday):
        sel_time = np.append(sel_time,thedata.timedata)
        sel_data = np.append(sel_data,thedata.data)
        if selcount == 0:
            plt.plot(thedata.timedata,thedata.data,'r',label='selected')
            selcount+=1
        else:
            plt.plot(thedata.timedata,thedata.data,'r')
    else:
#if y != goodyear[0] or int(d) not in goodday:
        nonsel_time = np.append(nonsel_time,thedata.timedata)
        nonsel_data = np.append(nonsel_data,thedata.data)
        if nonselcount == 0:
            plt.plot(thedata.timedata,thedata.data,'b',label='non selected')
            nonselcount+=1
        else:
            plt.plot(thedata.timedata,thedata.data,'b')
#plt.plot(d_array,h_array)

#plt.plot(sel_time, sel_data,'r.',label='selected')
#plt.plot(nonsel_time, nonsel_data,'b.',label='non selected')
plt.ylabel('radio baseline [ADC]')
plt.gcf().autofmt_xdate()
plt.legend(loc=2)
#plt.ylim(200,400)
plt.show()
