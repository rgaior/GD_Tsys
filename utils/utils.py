import numpy as np
import pickle
import scipy.signal as signal
###################################
### data loading/reading ##########
###################################
import ROOT as r
from root_numpy import root2rec
def loaddata(filename):
    a = root2rec(filename)
    return a

def loadparams(filename):
    fin = open(filename,'r')
    params = {}
    for l in fin:
        par =  np.asarray(l.split()[4:],dtype=float)
        stnr = int(l.split()[1])
        params[stnr] = par
    return params

def readtempfile(filename):
    f = open(filename,'r')
    hour = np.array([])
    min = np.array([])
    temp = np.array([])
    for l in f:
        lsplit = l.split()
        hour = np.append(hour,float(lsplit[0]))
        min = np.append(min,float(lsplit[1]))
        temp = np.append(temp,float(lsplit[2]))
    return [hour,min,temp]

def getdaysfromlist(filename):
    f = open(filename,'r')
    year = np.array([])
    day = np.array([])
    for l in f:
        lsplit = l.split()
        year = np.append(year,int(lsplit[0]))
        day = np.append(day,int(lsplit[1]))
    return [year,day]

def makeprofile(x,y,nrofbins,xmin,xmax):
#    from root_numpy import fill_profile
    p = r.TProfile("test","test",int(nrofbins),float(xmin),float(xmax))
#    p = r.TProfile(xmin,xmax)
    for xel,yel in zip(x,y):
        p.Fill(xel,yel)
        
    xout = np.array([])
    yout = np.array([])
    erryout = np.array([])
    errxout = np.array([])
    binsize = float((xmax-xmin))/nrofbins
    for ibin in range(1,nrofbins+1,1):
        ix = xmin + (ibin-1 + 0.5)*binsize
        xout = np.append(xout,ix)
        iy = p.GetBinContent(ibin)
        yout = np.append(yout,iy)
        ierr = p.GetBinError(ibin)
        erryout = np.append(erryout,ierr)
        errxout = np.append(errxout,binsize/2)
    return [xout,yout,errxout,erryout]


def readtwocolfile(file):
    f = open(file,'r')
    t = np.array([])
    controlperiod = np.array([])
    for l in f:
        t = np.append(t,float(l.split()[0]))
        controlperiod = np.append(controlperiod,int(l.split()[1]))
    f.close()
    return [t,controlperiod]

##################################
### conversion function ##########
##################################
from datetime import date
import datetime
def gpstodate(gpssecond):
    return  date.fromtimestamp(gpssecond+315964800)
import datetime
def gpstodatetime(gpssecond):
    return  datetime.datetime.fromtimestamp(gpssecond+315964800)

def tstamptodatetime(tstamp):
    return  datetime.datetime.utcfromtimestamp(tstamp)
def nptstamptodatetime(tstamp):
    date = np.array([])
    for t in tstamp:
        d = tstamptodatetime(t)
        date = np.append(date,d)
#        print t, ' ' , d
    return  date
 
def datettotimestamp(dt):
    if not isinstance(dt, datetime.date):
        print 'you should give a datetime.date or datetime.datetime instance in argument'
        return 
    elif not isinstance(dt, datetime.datetime):
        dt = datetime.datetime(dt.year,dt.month,dt.day) 
    timestamp = (dt - datetime.datetime(1970, 1, 1)).total_seconds()
    return timestamp
    
def datestringtodate(date):
#    print date[:4]
    y = int(date[:4])
    m = int(date[4:6])
    d = int(date[6:8])
    print 'y = ' ,y, ' m = ', m , ' d = ',d
    thedate = datetime.date(y,m,d)
    return thedate



def doytodate(year,doy,hour=None,minute=None):
    if hour == None and minute == None:
        date = datetime.datetime.strptime(str(year)+ ' '+str(doy), '%Y %j')
    elif minute == None:
        date = datetime.datetime.strptime(str(year)+ ' '+str(doy) + ' '+str(hour), '%Y %j %H')
    else:
        date = datetime.datetime.strptime(str(year)+ ' '+str(doy) + ' '+str(hour)+ ' ' +str(minute) , '%Y %j %H %M')
    return date
    
def datetodoy(date):
    year = date.year
    day = int(date.strftime('%j'))
    return (year,day)

def doytoUTC(year,doy,hour=None,minute=None):
    date = doytodate(year,doy,hour,minute)
    tstamp = datettotimestamp(date)
    return tstamp

def UTCtodoy(utc):
    date = tstamptodatetime(utc)
    return datetodoy(date)

def hhmmtosecond(hhmm):
    hh = hhmm/100
    mm = hhmm % 100
    sec = hh*3600 + mm*60
    return sec

def sectohhmm(sec):
    h = int(sec/3600)
    m = int(sec%3600)
    hhmm = h*100+m
    return hhmm

def hhmmtohour(hhmm):
    hh = hhmm/100
    mm = hhmm % 100
    sec = hh*3600 + mm*60
    h = hh + mm.astype(float)/60
    return h

def hourtohhmm(hours):
    hh = hours/100
    hh = hh.astype(int)
    mm = (hours-hh*100)*0.6
    hhmm = hh*100 + mm
    return hhmm

def timetohour(hour,min):
    return hour + float(min)/60.


###########################################
####   data selection/correction      #####
###########################################

def correctwithfct(data, x, function):
#    print function
    newdata = data - function(x)
    return newdata

def correctwithpoly(data,x,poly):
#    print function
    newdata = data - poly(x)
    return newdata

def kinkfcn(params, x, data):                                                                            
    a1 = params['a1'].value 
    b1 = params['b1'].value                                                                    
    a2 = params['a2'].value
    b2 = params['b2'].value
    t = params['t'].value
    x1 = x[x<t]
    x2 = x[x>=t]
    y1 = a1*x1 + b1
    y2 = a2*x2 + b2
    model = np.array([])
    model = np.append(model,y1)
    model = np.append(model,y2)
    return model - data
###############################################
####              filtering               #####
###############################################

def lowpass(amp, sampling, order, fcut):
    Nyfreq = sampling/2
    ratiofcut = float(fcut)/Nyfreq
    b, a = signal.butter(order, ratiofcut, 'low')
    filtered = signal.filtfilt(b, a, amp)
    return filtered

def lowpasshard(amp, sampling, fcut):
    fft = np.fft.rfft(amp)
    freq = np.fft.rfftfreq(len(fft),float(1./sampling))
    Nyfreq = sampling/2
    #    print 'Nyfreq = ' , Nyfreq, 'fcut = ', fcut
    min = np.min(np.absolute(fft))
    ratiofcut = float(fcut)/float(Nyfreq)
    size = len(fft)
    newpass = fft[:int(ratiofcut*size)]
    sizeofzeros = size - len(newpass)
    newcut = np.zeros(sizeofzeros)
    newfft = np.append(newpass,newcut)
    out = np.fft.irfft(newfft)
    return out.real

def highpass(amp, sampling, order, fcut):
    Nyfreq = sampling/2
    ratiofcut = float(fcut)/Nyfreq
    b, a = signal.butter(order, ratiofcut, 'high')
    filtered = signal.filtfilt(b, a, amp)
    return filtered

def highpasshard(amp, sampling, fcut):
    fft = np.fft.rfft(amp)
    freq = np.fft.rfftfreq(len(fft),float(1./sampling))
    Nyfreq = sampling/2
    min = np.min(np.absolute(fft))
    ratiofcut = float(fcut)/Nyfreq
    size = len(fft)
    newpass = fft[int(ratiofcut*size):]
    sizeofzeros = size - len(newpass)
    newcut = np.zeros(sizeofzeros)
    newfft = np.append(newpass,newcut)
    out = np.fft.irfft(newfft)
    return out.real


def slidingwindow(y,bins,option=None):
    window = np.ones(bins)/bins
    if option is not None:
        if option.lower() not in ['full','same','valid']:
            print 'invalid option, check your sliding window'
    if option == None:
        return np.convolve(y,window,'same')
    else:
        return np.convolve(y,window,option)

def matchedfilter(y,sig,option=None):
    filtered = signal.correlate(y,sig, mode='full')
#    filtered = signal.correlate(y,sig, mode='valid')
    return filtered


def getxforconv(x,size):
    diff = len(x) - size
    if (diff==0):
        return x
    else:
        return x[diff/2:-diff/2]


def gauss(x,a,b,c):
    g = a*np.exp( -((x-b)**2)/(2*c**2) )
    return g

def issimilarfit(fit1,fit2,tol):
    sim = True
    for f1,f2 in zip(fit1,fit2):
        if f2 > (1+tol)*f1 or f2 < (1-tol)*f1:
            sim = False
    return sim

##############################
## sun transit function   ####
##############################
import pickle
def getsunmax(date):
    sundatafile = '/Users/romain/work/Auger/EASIER/LPSC/monitoring/data/sundata/sundata.pkl'
    pkl_file = open(sundatafile, 'rb')
    sunflux = pickle.load(pkl_file)
    pkl_file.close()
    flux = sunflux[date]
    print 'flux = ', flux
    return flux

def suntemptoadc(tsys,temp):
    adc =10*np.log10( (tsys+temp)/tsys )*50
    return adc 

def adctotsys(adc,temp):
    adctodb = 50
    print adc
    print temp
    tsys = temp/( np.power(10, adc/(adctodb*10)) - 1 )
    print tsys
    return tsys

def expofunc0(x, a, sigma,mu,c):
    return a*np.exp(-((x - mu)/sigma)**2) + c

def expofunc1(x, a, sigma,mu,b,c):
    return a*np.exp(-((x - mu)/sigma)**2) + b*x + c

def expofunc(x, a, sigma,mu,b,c):
    return a*np.exp(-((x - mu)/sigma)**2) + b*x + c

def expofunc2(x, a, sigma,mu,b,c, d):
    return a*np.exp(-((x - mu)/sigma)**2) + b*x**2 + c*x + d

def expofunctwolinear(x, a, sigma, mu, p00, p01, p10, p11):
    hmax = 1700
    res = np.array([])
    res = a*np.exp(-((x[x<hmax] - mu)/sigma)**2) + p00 + p01*x[x<hmax]
    res = np.append(res,a*np.exp(-((x[x>=hmax] - mu)/sigma)**2) + p10 + p11*x[x>=hmax])
    return res
        

def getexpectedtemp(file):
    #    file format: y m d maxadc
    f = open(file,'r')
    day = {}
    for l in f:
        ls = l.split()
        d = (int(ls[0]),int(ls[1]),int(ls[2]))
        temp = float(ls[3])
        day[d] = temp
    return day


def checkfit(file, tol, date,timeofmax,width):
    f = open(file,'rb')
    dict = pickle.load(f)
    timeofmaxth = dict[date][0] + 300
    widthth = dict[date][1]
    ok = True
#    print 'timeofmaxth = ', timeofmaxth , ' widthth = ', widthth
#    print 'timeofmax = ', timeofmax , ' width = ', width
    if width <widthth -tol*widthth or  width > widthth + tol*widthth:
        ok = False
    if timeofmax <timeofmaxth -tol*timeofmaxth or  timeofmax > timeofmaxth + tol*timeofmaxth:
        ok = False
    return ok

def getexpected(file,date):
    f = open(file,'rb')
    dict = pickle.load(f)
    timeofmaxth = dict[date][0] + 300
    widthth = dict[date][1]
#    print 'timeofmaxth = ', timeofmaxth , ' widthth = ', widthth
    return [timeofmaxth,widthth]


def getthetaphifromfilename(fname):
    filename = fname.split('/')[-1]
    theta = getangle(filename,'theta')
    phi = getangle(filename,'phi')
    return [theta,phi]
#     print filename
#     centertheta = isolatestr(filename,'theta_','_')
#     centerphi = isolatestr(filename,'phi_','_')
#     print centertheta
#     print centerphi
#    theta_-20_-2.5_phi_0_-2.5_expparam_2016_orteguina
    
def getangle(fname,thetaorphi):
    if thetaorphi == 'theta':
        ada = isolatestr(fname,'theta','_phi')
    elif thetaorphi == 'phi':
        ada = isolatestr(fname,'_phi','_exp')
    print ada
    a = isolatestr(ada,'_','_')
    da = ada[ada.rfind('_')+1:]
    angle = float(a) + float(da)
    return angle

def isolatestr(value, a, b):
    # Find and validate before-part.
    pos_a = value.find(a)
    if pos_a == -1: return "vvv"
    # Find and validate after part.
#    pos_b = value.find(b)
    pos_b = value.find(b,pos_a+len(a))
    if pos_b == -1: return "bbb"
    # Return middle part.
    adjusted_pos_a = pos_a + len(a)
    if adjusted_pos_a >= pos_b: return "nnn"
    return value[adjusted_pos_a:pos_b]

def between(value, a, b):
    # Find and validate before-part.
    pos_a = value.find(a)
    if pos_a == -1: return ""
    # Find and validate after part.
    pos_b = value.rfind(b)
    if pos_b == -1: return ""
    # Return middle part.
    adjusted_pos_a = pos_a + len(a)
    if adjusted_pos_a >= pos_b: return ""
    return value[adjusted_pos_a:pos_b]


#def minutedif(min1,nin2):
    
