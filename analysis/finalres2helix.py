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

def fitwithexpo1(x, y, yerr, a, sigma,mu,b,c):
    try:
        popt, pcov = curve_fit(utils.expofunc1,x,y,sigma=yerr, p0=[a,sigma,mu,b,c])
    except RuntimeError:
        print("Error - curve_fit failed")
        return [0.0,0.0]
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
    deltamax = 2
    if popt[0] < 0:
        print ' !!!!!!!!!!!!!!!!!! NEG PEAK !!!!!!!!!!!!!!!!!!!!!!!!', popt 
        good = False
    if popt[2] > tofmax + deltamax or  popt[2] < tofmax - deltamax:
        print ' !!!!!!!!!!!!!!!!!! Time of MAX !!!!!!!!!!!!!!!!!!!!!!!!', popt
#        print 'popt[2] ' , popt[2]
        good = False
    if np.abs(popt[1]) > 3 or np.abs(popt[1]) < 0.3:
        print ' !!!!!!!!!!!!!!!!!! SIGMA !!!!!!!!!!!!!!!!!!!!!!!!', popt
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
stid = constant.Helixstationidbyname[stname]
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

#aeffuncertname = constant.aeffuncertfolder + '/aeffuncert_' + stname + '.npz'
#aeffuncert = np.load(aeffuncertname)

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
for y,d in zip(goodyear,goodday):
#for y,d in zip(goodyear[:30],goodday[:30]):
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
    nsradioc = nsradio
    
    ## uncertainty:
    uncert = np.std(nsradioc)
    relaeffuncert =0
    hourarray = gethourarray(time)
    fithourarray = hourarray[ np.where((hourarray > timeofsunmax -6) & (hourarray < timeofsunmax + 6 ) )]
#    print nspopt[0]

    fitradio =  radio[np.where((hourarray > timeofsunmax - 6 ) & (hourarray < timeofsunmax + 6 ) ) ]

    [popt,pcov] = fitwithexpo2(fithourarray, fitradio, uncert, 15, 1.5,timeofsunmax,0,0,0)
#    [popt,pcov] = fitwithexpo2(hourarray, radio, uncert, 15, 1.5,timeofsunmax,0,0,0)
#    [popt,pcov] = fitwithexpo0(hourarray, radio, uncert, 10, 1.5,timeofsunmax,0)
    if type(pcov)==type(1.1):
        continue
    perr = np.sqrt(np.diag(pcov))
    erronmax = perr[0]
    erronhmax = perr[2]
    
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
#     fig = plt.figure(figsize=(12,6))
#     fig.suptitle(stname + ' (' + t1.strftime('%d %b %Y')+')',fontweight='bold',fontsize=15)
#     plt.plot(hourarray,radio,'b.',label='data')
#     plt.plot(fithourarray,rfit,'r',lw=2,label='fit')
#     plt.xlabel('local time')
#     plt.ylabel('corrected baseline [ADC]')
    if goodfit == True:
        nsfithourarray = hourarray[ np.where(hourarray < timeofsunmax -3)]
        nsfitradio = radio[np.where(hourarray < timeofsunmax - 3)]
#        [nspopt,nspcov] = fitwithexpo2(nsfithourarray, nsfitradio, uncert, 10, 1.5,5,0,0,0)
        [nspopt,nspcov] = fitwithexpo1(nsfithourarray, nsfitradio, uncert, 10, 1.5, 5, 0, 0)

#         fig = plt.figure()
#         fig.suptitle(stname + ' (' + t1.strftime('%d %b %Y')+')',fontweight='bold',fontsize=15)
#         plt.plot(hourarray,radio,'b.',label='data')
#         plt.plot(fithourarray,rfit,'r',lw=2,label='fit')
#         plt.xlabel('time [UTC]')
        #         plt.ylabel('corrected baseline [ADC]')

        if type(nspopt) == type(1.0):
            nsmod = 0
        else:
#            nsrfit = utils.expofunc2(nsfithourarray,nspopt[0],nspopt[1],nspopt[2],nspopt[3],nspopt[4],nspopt[5])
            nsrfit = utils.expofunc1(nsfithourarray,nspopt[0],nspopt[1],nspopt[2],nspopt[3],nspopt[4])
#            plt.plot(nsfithourarray,nsrfit,'g',lw=2,label='fit')
            nsmod = nspopt[0]
        print ' !!!!!!!!!!!!!!!!! nsmod = ', nsmod


        

        datestring = t1.strftime('%d %b %Y')
#        print datestring , ' ', maxfit
        erronmaxbyhand = 15
        erronhmaxbyhand = 0.5
        datearray.append(datestring)
        adate = np.append(adate,datetime.datetime(t1.year,t1.month,t1.day))
        asimhmax = np.append(asimhmax,poptsim[2])
        asimmax = np.append(asimmax, maxsim)
        aerrmax = np.append(aerrmax,erronmax +erronmaxbyhand)
        aerrhmax = np.append(aerrhmax,erronhmax +erronhmaxbyhand)
        tsysall = np.append(tsysall,tsys)
        #   fig2 = plt.figure(figsize=(12,6))
        #   plt.plot(time,radio,'b.',label='non selected')
        #   plt.legend()
        relaeff = 0
#        relaeff = relaeffuncert
#        relfsun = 0.05
        relfsun = 0.0
#        sigdeltaP = (erronmax + np.absolute(nsmod))/adctop
        sigdeltaP = (erronmax + erronmaxbyhand)/adctop
#        sigdeltaP = erronmax/adctop
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
#plt.errorbar(x,tsysall, yerr=sigmatsysstat,color='y', fmt='o',lw=4, alpha=0.9,label='stat. only')
plt.errorbar(x,tsysall, yerr=sigmatsysstat,color='b', fmt='o',lw=4, alpha=0.9)
plt.fill_between(x, theres-errtheres, theres +  errtheres,facecolor='red',alpha=0.5)
plt.xlim(np.min(x) - 0.5, np.max(x) + 0.5 )
plt.xticks(x,myticks,rotation='vertical',size=10)
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
plt.xlim(np.min(x) - 0.5, np.max(x) + 0.5 )
plt.xticks(x,myticks,rotation='vertical',size=10)

figmax = plt.figure(figsize=(15,5))
#figres.suptitle(stname,fontsize=15,fontweight='bold')
figmax.subplots_adjust(left=0.1, bottom=0.3, right=0.94, top=0.9,
                       wspace=None, hspace=None)
plt.errorbar(x, amax, yerr=aerrmax,fmt='o',label='data')
t1 = 100
t2 = 150
t3 = 200
ex1 =adctop*10*np.log10(1+asimmax/t1)
ex2 =adctop*10*np.log10(1+asimmax/t2)
ex3 =adctop*10*np.log10(1+asimmax/t3)
plt.plot(x,ex1,lw=2,label='sim. T_sys = ' +str(t1) + ' K')
plt.plot(x,ex2,lw=2,label='sim. T_sys = ' +str(t2) + ' K')
plt.plot(x,ex3,lw=2,label='sim. T_sys = ' +str(t3) + ' K')
plt.xlim(np.min(x) - 0.5, np.max(x) + 0.5 )
plt.xticks(x,myticks,rotation='vertical',size=10)
plt.ylabel('radio max [ADC]')
plt.legend()

resfile = outf + '/results.txt' 
fout = open(resfile,'w')
fighmax.savefig('/Users/romain/work/Auger/EASIER/LPNHE/notes/tex/tsyshelix/plots/' + 'hmax_'+ stname+'.png')
figmax.savefig('/Users/romain/work/Auger/EASIER/LPNHE/notes/tex/tsyshelix/plots/' + 'max_'+ stname+'.png')
figres.savefig('/Users/romain/work/Auger/EASIER/LPNHE/notes/tex/tsyshelix/plots/' + 'tsys_'+ stname+'.png')

####################################################
######### put the three figures together ###########
####################################################
figall, (ax1,ax2,ax3) = plt.subplots(3,1,figsize=(12,10),sharex=True)
ax1.errorbar(x,tsysall, yerr=sigmatsysstat,color='b', fmt='o',lw=4, alpha=0.9)
ax1.fill_between(x, theres-errtheres, theres +  errtheres,facecolor='red',alpha=0.5)
ax1.set_ylabel('system temperature [K]')
ax1.set_ylim(ymin,ymax)
#ax1.text(2, 100, r'an equation: $E=mc^2$', fontsize=15)
ax1.text(2, ymax-20, textstr, fontsize=17, fontweight ='bold', verticalalignment='top',horizontalalignment='left', bbox=props)

ax2.errorbar(x, amax, yerr=aerrmax,fmt='o',label='data')
ax2.plot(x,ex1,lw=2,label='sim. T_sys = ' +str(t1) + ' K')
ax2.plot(x,ex2,lw=2,label='sim. T_sys = ' +str(t2) + ' K')
ax2.plot(x,ex3,lw=2,label='sim. T_sys = ' +str(t3) + ' K')
ax2.set_ylim(np.min(amax - aerrmax), np.max(amax + aerrmax) + 20)
ax2.set_ylabel('radio max [ADC]')
ax2.legend(ncol=2)

ax3.errorbar(x,ahmax,yerr=aerrhmax, fmt='o',)
ax3.plot(x,asimhmax,lw=2)
ax3.set_ylabel('time of max [hour]')


ax3.set_xlim(np.min(x) - 0.5, np.max(x) + 0.5 )
plt.xticks(x,myticks,rotation='vertical',size=10)

figall.savefig('/Users/romain/work/Auger/EASIER/LPNHE/notes/tex/tsyshelix/plots/' + 'all_'+ stname+'.png')


for i in range(len(tsysall)):    
    fout.write(str(adate[i].year) + ' ' + str(adate[i].month) + ' ' + str(adate[i].day) + ' ' + str(tsysall[i]) + ' ' + str(sigmatsysall[i])  + ' ' + str(ahmax[i]) + ' ' + str(aerrhmax[i]) + '\n' )

plt.show()

