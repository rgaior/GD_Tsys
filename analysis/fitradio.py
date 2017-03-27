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


def fitwithexpo0(x, y, yerr, a, sigma,mu,c):
    try:
        popt, pcov = curve_fit(utils.expofunc0,x,y,sigma=yerr, p0=[a,sigma,mu,c])
    except RuntimeError:
        print("Error - curve_fit failed")
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
parser.add_argument("fluxorigin", type=str, nargs='?',default='nobeyama', help="canadian or nobeyama")
parser.add_argument("-save", action='store_true', help="save the fit")
parser.add_argument("-list", type=str, help="list of good days")
args = parser.parse_args()
stname = args.stname
savefit = args.save
flux = args.fluxorigin
if args.list:
    listfile = args.list
    goodlistname = listfile
else:
    goodlistname = stname + '.txt'

stid = constant.GDstationidbyname[stname]
[goodyear,goodday] = utils.getdaysfromlist(constant.listfolder + goodlistname)
print goodday

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

temp = np.array([])
radio = np.array([])
nosuntime = np.array([])
nosunradio = np.array([])
nosuntemp = np.array([])
tsysall =  np.array([])
sigmatsysall =  np.array([])
sigmatsysstat =  np.array([])
for y,d in zip(goodyear,goodday):
#for y,d in zip(goodyear[:10],goodday[:10]):
    outf = outfolder + str(stid)
    datafilename = outf + '/data_' + str(int(y)) + '_' + str(int(d)) + '.pkl'
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
    
    t1 = utils.doytodate(int(y),int(d),int(timeofsunmax - 2), int(np.modf(timeofsunmax - 2)[0]))
    t2 = utils.doytodate(int(y),int(d),int(timeofsunmax + 2), int(np.modf(timeofsunmax + 2)[0]))
    nsradio = thedata.data[(thedata.timedata < t1)  | (thedata.timedata > t2) ]
    nstemp = thedata.tempLL[(thedata.timedata < t1)  | (thedata.timedata > t2) ]
    ## subtract the mean of the time when no sun is present in the data
    radio = thedata.data - np.mean(nsradio)
    temp = thedata.tempLL - np.mean(nstemp)
    ## get the file with temperature fit result
    fitfile = constant.resultfolder + str(stid) + '/tempfit.npy'
    tempfit = np.load(fitfile)
    pfit = np.poly1d(tempfit)
    radio = radio - pfit(temp)
    ## baseline corrected
    nsradioc = nsradio - np.mean(nsradio)
    nsradioc = nsradioc -  pfit(nstemp)
    ## uncertainty:
    uncert = np.std(nsradioc)
    
    hourarray = gethourarray(time)
    
    [popt,pcov] = fitwithexpo0(hourarray, radio, uncert, 10, 1.5,timeofsunmax,0)
    perr = np.sqrt(np.diag(pcov))
    erronmax = perr[0]
    print 'erronmax = ', erronmax , 'uncert = ' , uncert
    rfit = utils.expofunc0(hourarray,popt[0],popt[1],popt[2],popt[3])
    
    goodfit = isgoodfit(stid,popt,timeofsunmax)   
    maxfit = popt[0]
    tsys = maxsim/(np.power(10,maxfit/500.) -1)
    hourarray = hourarray - 3
    hourarray = hourarray % 24
    #    rfit = np.roll(rfit,-sizetoroll)
    
    #     if goodfit == False:
#         fig = plt.figure(figsize=(12,6))
#         fig.suptitle(stname + ' (' + t1.strftime('%d %b %Y')+')',fontweight='bold',fontsize=15)
#         plt.plot(hourarray,radio,'b.',label='data')
#         plt.plot(hourarray,rfit,'r',lw=2,label='fit')
#         plt.xlabel('local time')
#         plt.ylabel('corrected baseline [ADC]')
    if goodfit == True:
#         fig = plt.figure(figsize=(12,6))
#         fig.suptitle(stname + ' (' + t1.strftime('%d %b %Y')+')',fontweight='bold',fontsize=15)
#         plt.plot(hourarray,radio,'b.',label='data')
#         plt.plot(hourarray,rfit,'r',lw=2,label='fit')
#         plt.xlabel('local time')
#         plt.ylabel('corrected baseline [ADC]')
        print popt
        tsysall = np.append(tsysall,tsys)
        
        #   fig2 = plt.figure(figsize=(12,6))
    #   plt.plot(time,radio,'b.',label='non selected')    
    #   plt.legend()
        relaeff = 0.1
        relfsun = 0.0
        sigdeltaP = erronmax/50
#        sigdeltaP = uncert/50
        deltaP = maxfit/50
        
        ln10on10 = np.log(10)/10
        reltsys = np.sqrt(relaeff**2 + relfsun**2 + ( (ln10on10* np.power(10,deltaP/10)) / (np.power(10,deltaP/10) - 1) )**2 *sigdeltaP**2)
        reltsysstat = np.sqrt( ( (ln10on10* np.power(10,deltaP/10)) / (np.power(10,deltaP/10) - 1) )**2 *sigdeltaP**2)

        sigmatsysall = np.append(sigmatsysall,reltsys*tsys)
        sigmatsysstat = np.append(sigmatsysstat,reltsysstat*tsys)
        
## compute weighted mean and error
oneoversigmasquare = 1/(sigmatsysall*sigmatsysall)
xoversigmasquare = tsysall*oneoversigmasquare
wmean = np.sum(xoversigmasquare)/np.sum(oneoversigmasquare)
errorwmean = np.sqrt(1/np.sum(oneoversigmasquare))
print 'final result is ... Tsys = ' , wmean , ' +- ',  errorwmean
x = np.arange(0,len(tsysall),1)
figres = plt.figure(figsize=(10,8)) 
plt.errorbar(x, tsysall, yerr=sigmatsysall,label='stat + sys')
plt.errorbar(x, tsysall, yerr=sigmatsysstat, lw=2, fmt='o',label='stat only')
#bins = np.linspace(0,300,150)
#plt.hist(tsysall,bins=bins)
plt.legend()
plt.show()
