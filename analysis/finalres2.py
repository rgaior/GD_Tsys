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
#mpl.rcParams['ytick.labelsize'] = 10
#mpl.rcParams['xtick.labelsize'] = 10

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
    if popt[0] < 0:
        good = False
    if popt[2] > tofmax + deltamax or  popt[2] < tofmax - deltamax:
#        print 'popt[2] ' , popt[2]
        good = False
    if np.abs(popt[1]) > 2 or np.abs(popt[1]) < 0.3:
#        print 'popt[1] ',  popt[1]
        good = False
    return good


############################
##  argument parser       ##
############################
parser = argparse.ArgumentParser()
parser.add_argument("stname", type=str, nargs='?',default='popey', help="station name")
parser.add_argument("fluxorigin", type=str, nargs='?',default='nobeyama', help="canadian or nobeyama")
parser.add_argument("--deltatheta", type=int, nargs='?',default='0', help="delta theta")
parser.add_argument("--deltaphi", type=int, nargs='?',default='0', help="delta phi")
parser.add_argument("-save", action='store_true', help="save the fit")
parser.add_argument("-list", type=str, help="list of good days")

args = parser.parse_args()
delt = args.deltatheta
delp = args.deltaphi
stname = args.stname
savefit = args.save
flux = args.fluxorigin
if args.list:
    listfile = args.list
    goodlistname = listfile
else:
    goodlistname = stname + '.txt'

print ' delt = ' ,delt ,' and delp = ', delp 
datafolder = '/t'+str(delt) + 'p'+str(delp)+'/'

stname = args.stname
stid = constant.GDstationidbyname[stname]
goodlistname = stname + '.txt'
[goodyear,goodday] = utils.getdaysfromlist(constant.listfolder + goodlistname)
if flux == 'nobeyama':
    outfolder = constant.resultfolder2
elif flux == 'canadian':
    outfolder = constant.resultfolder
elif flux == 'test':
    outfolder = constant.resultfoldertest + datafolder
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

aeffuncertname = constant.aeffuncertfolder + '/aeffuncert_' + stname + '.npz'
aeffuncert = np.load(aeffuncertname)
#print aeffuncert['arr_0']
datearray =[]
asimhmax = np.array([])
asimmax = np.array([])
ahmax = np.array([])
amax = np.array([])
aerrhmax = np.array([])
aerrmax = np.array([])
adate = np.array([])

outf = outfolder + str(stid)
for y,d in zip(goodyear[2::2],goodday[2::2]):
#for y,d in zip(goodyear[:30],goodday[:30]):
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
    fithourarray = hourarray[ np.where( (radio< 50) & (hourarray > timeofsunmax -4) & (hourarray < timeofsunmax + 4 ) )]
    fitradio = radio[np.where( (radio< 50) & (hourarray > timeofsunmax -4 ) & (hourarray < timeofsunmax + 4 ) ) ]

    [popt,pcov] = fitwithexpo2(fithourarray, fitradio, uncert, 15, 1.5,timeofsunmax,0,0,0)
#    [popt,pcov] = fitwithexpo2(hourarray, radio, uncert, 15, 1.5,timeofsunmax,0,0,0)
#    [popt,pcov] = fitwithexpo0(hourarray, radio, uncert, 10, 1.5,timeofsunmax,0)
    print '[popt,pcov] = ' , popt, ' ', type(pcov)
    if type(pcov)==type(1.1):
        continue
    perr = np.sqrt(np.diag(pcov))
    erronmax = perr[0]
    erronhmax = perr[2]
    
    erronmax += 5
    erronhmax += 0.2
    rfit = utils.expofunc2(fithourarray,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5])
#    rfit = utils.expofunc0(hourarray,popt[0],popt[1],popt[2],popt[3])
    goodfit = isgoodfit(stid,popt,timeofsunmax)
    maxfit = popt[0]
    adctop = 48.7
    tsys = maxsim/(np.power(10,maxfit/(10*adctop)) -1)
### here
#    hourarray = hourarray - 3
#    hourarray = hourarray % 24

    #    rfit = np.roll(rfit,-sizetoroll)

#     if goodfit == False:
#         fig = plt.figure(figsize=(12,6))
#         fig.suptitle(stname + ' (' + t1.strftime('%d %b %Y')+')',fontweight='bold',fontsize=15)
#         plt.plot(hourarray,radio,'b.',label='data')
#         plt.plot(fithourarray,rfit,'r',lw=2,label='fit')
#         plt.xlabel('local time')
#         plt.ylabel('corrected baseline [ADC]')
    if goodfit == True:
#         fig = plt.figure()
#         fig.suptitle(stname + ' (' + t1.strftime('%d %b %Y')+')',fontweight='bold',fontsize=15)
#         plt.plot(hourarray,radio,'b.',label='data')
#         plt.plot(fithourarray,rfit,'r',lw=2,label='fit')
#         plt.xlabel('time [UTC]')
        #         plt.ylabel('corrected baseline [ADC]')
        datestring = t1.strftime('%d %b %Y')
        print datestring , ' ', maxfit
        datearray.append(datestring)
        adate = np.append(adate,datetime.datetime(t1.year,t1.month,t1.day))
        asimhmax = np.append(asimhmax,poptsim[2])
        asimmax = np.append(asimmax, maxsim)
        aerrmax = np.append(aerrmax,erronmax)
        aerrhmax = np.append(aerrhmax,erronhmax)
        tsysall = np.append(tsysall,tsys)
        #   fig2 = plt.figure(figsize=(12,6))
        #   plt.plot(time,radio,'b.',label='non selected')
        #   plt.legend()
        relaeff = 0
#        relaeff = relaeffuncert
#        relfsun = 0.05
        relfsun = 0.0
        sigdeltaP = erronmax/adctop
        deltaP = maxfit/adctop
        ahmax = np.append(ahmax,popt[2])
        amax = np.append(amax,popt[0])
        ln10on10 = np.log(10)/10
        reltsys = np.sqrt(relaeff**2 + relfsun**2 + ( (ln10on10* np.power(10,deltaP/10)) / (np.power(10,deltaP/10) - 1) )**2 *sigdeltaP**2)
        reltsysstat = np.sqrt( ( (ln10on10* np.power(10,deltaP/10)) / (np.power(10,deltaP/10) - 1) )**2 *sigdeltaP**2)
        sigmatsysall = np.append(sigmatsysall,reltsys*tsys)
        sigmatsysstat = np.append(sigmatsysstat,reltsysstat*tsys)

## compute weighted mean and error
tant = 0
oneoversigmasquare = 1/(sigmatsysall*sigmatsysall)
xoversigmasquare = tsysall*oneoversigmasquare
wmean = np.sum(xoversigmasquare)/np.sum(oneoversigmasquare)
errorwmean = np.sqrt(1/np.sum(oneoversigmasquare))
tsys = wmean + tant
tsysall = tsysall + tant

print 'final result is ... Tsys = ' , tsys , ' +- ',  errorwmean
y = np.arange(0,len(tsysall),1)
x = np.arange(0,len(tsysall),1)
myticks = datearray
figres = plt.figure(figsize=(15,5))
figres.subplots_adjust(left=0.1, bottom=0.3, right=0.94, top=0.9,
                       wspace=None, hspace=None)
figres.suptitle(stname,fontsize=15,fontweight='bold')
#ax = plt.subplot(111)
theres = (tsys)*np.ones(len(y))
errtheres = errorwmean*np.ones(len(y))
ymin = np.min(tsysall-sigmatsysall)
ymax = np.max(tsysall+sigmatsysall)
#plt.errorbar(x,tsysall, yerr=sigmatsysall,lw=4, fmt='o',alpha=0.9,label='stat. and sys.')
plt.errorbar(x,tsysall, yerr=sigmatsysstat,color='b', fmt='o',lw=4, alpha=0.6)
#plt.errorbar(x,tsysall, yerr=sigmatsysstat,color='b', fmt='o',lw=4, alpha=0.6,label='stat. only')

plt.xticks(x,myticks,rotation='vertical',size=10)
plt.fill_between(x, theres-errtheres, theres +  errtheres,facecolor='red',alpha=0.5)
plt.ylabel('system temperature [K]')
plt.ylim(ymin,ymax)
plt.legend()
#plt.xlim(35,185)
#plt.ylim(y[0]-0.5,y[-1]+0.5)
#textstr = '\n $T_{sys}=%.1f \pm%.1f [K]$'%(tsys, errorwmean)
textstr = '$T_{sys}=%.1f \pm%.1f [K]$'%(tsys, errorwmean)
#textstr = '$testo$'
props = dict(boxstyle='round', facecolor='white')
#plt.gca().text(0.95, 0.9, textstr, transform=plt.gca().transAxes, fontsize=17, fontweight ='bold', verticalalignment='top',horizontalalignment='right', bbox=props)
plt.gca().text(0.05, 0.9, textstr, transform=plt.gca().transAxes, fontsize=17, fontweight ='bold', verticalalignment='top',horizontalalignment='left', bbox=props)

#bins = np.linspace(0,300,150)
#plt.hist(tsysall,bins=bins)

fighmax = plt.figure(figsize=(15,5))
#figres.suptitle(stname,fontsize=15,fontweight='bold')
fighmax.subplots_adjust(left=0.1, bottom=0.3, right=0.94, top=0.9,
                       wspace=None, hspace=None)
plt.errorbar(x,ahmax,yerr=aerrhmax, fmt='o',)
plt.plot(x,asimhmax,lw=2)
plt.ylabel('time of max [hour]')
plt.xticks(x,myticks,rotation='vertical',size=10)

figmax = plt.figure(figsize=(15,5))
#figres.suptitle(stname,fontsize=15,fontweight='bold')
figmax.subplots_adjust(left=0.1, bottom=0.3, right=0.94, top=0.9,
                       wspace=None, hspace=None)
plt.errorbar(x, amax, yerr=aerrmax,fmt='o',label='data')
t1 = 30
t2 = 50
t3 = 70
ex1 =adctop*10*np.log10(1+asimmax/t1)
ex2 =adctop*10*np.log10(1+asimmax/t2)
ex3 =adctop*10*np.log10(1+asimmax/t3)
plt.plot(x,ex1,lw=2,label='sim. T_sys = ' +str(t1) + ' K')
plt.plot(x,ex2,lw=2,label='sim. T_sys = ' +str(t2) + ' K')
plt.plot(x,ex3,lw=2,label='sim. T_sys = ' +str(t3) + ' K')
plt.xticks(x,myticks,rotation='vertical',size=10)
plt.ylabel('radio max [ADC]')
plt.legend()

resfile = outf + '/results.txt' 
fout = open(resfile,'w')
for i in range(len(tsysall)):    
    fout.write(str(adate[i].year) + ' ' + str(adate[i].month) + ' ' + str(adate[i].day) + ' ' + str(tsysall[i]) + ' ' + str(sigmatsysall[i])  + ' ' + str(ahmax[i]) + ' ' + str(aerrhmax[i]) + '\n' )

plt.show()
