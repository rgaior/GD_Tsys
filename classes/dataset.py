import utils
import numpy as np
import datetime 
from scipy.optimize import curve_fit
import radiotemp
datafolder = '/Users/romain/work/Auger/EASIER/LPSC/monitoring/data/'
correctionfile = datafolder + 'newcorrections.txt' #missing a few stations: chuletas, chape and gato is strange
class Dataset:
    def __init__(self, name = '',date1 =None,date2 =None):
        self.name = name
        self.date1 = None
        self.date2 = None
        if date1 is not None:
            self.date1 = date1
        if date2 is not None:
            self.date2 = date2
        self.data = []
        self.radio = np.array([])
        self.varradio = np.array([])
        #corrected radio after a fit
        self.radioc= np.array([])
        self.tempLL = np.array([])
        self.pressureLL = np.array([])
        self.tempelec = np.array([])
        self.humLL = np.array([])
        self.hhmm = np.array([])
        self.date = np.array([])
        self.time = np.array([])
        self.solarV = np.array([])
        self.UBV = np.array([])
        # voltage corrected value of radio power
        self.powerinput = np.array([])
        self.solarI = np.array([])
        self.period = np.array([])
        self.fit = np.array([])
        self.params = [] # array of tuple for correction parameters
#        self.params = np.array([])

    def loaddata(self):
        #first get all the tree data
        self.data = utils.loaddata(self.name)
        

    def loadperiod(self,pname):
        timeperiod = utils.readtwocolfile(pname)
        #        self.period = np.interp(self.getdata('Time'),timeperiod[0],timeperiod[1])
        self.period = np.interp(self.getdata('Time'),timeperiod[0],timeperiod[1])
        print timeperiod[0], ' ', timeperiod[1]
        print self.getdata('Time')
#        print self.time[(self.period==1)]

    def loadparams(self):
        self.params = utils.loadparams(correctionfile)
#        print self.params

    def selectleafs(self,station=None):
        if self.date1 is None:
            basiccond = np.where(
#                (self.data['LSId']==station)
                (np.invert(np.isnan(self.data['FDLLOutsideTemp']) ))
                &(np.invert(np.isinf(self.data['FDLLOutsideTemp']) ))
                &(np.invert(np.isnan(self.data['LLHumidity']) ))
                &(np.invert(np.isinf(self.data['LLHumidity']) ))
                )
        elif self.date2 is None:
            basiccond = np.where(
#                (self.data['LSId']==station)
                (self.data['yymmdd'] > date1)
                &(np.invert(np.isnan(self.data['FDLLOutsideTemp']) ))
                &(np.invert(np.isinf(self.data['FDLLOutsideTemp']) ))
                &(np.invert(np.isnan(self.data['LLHumidity']) ))
                &(np.invert(np.isinf(self.data['LLHumidity']) ))
                )
        else:
            basiccond = np.where(
#                (self.data['LSId']==station)
                (self.data['yymmdd'] > self.date1)
                & (self.data['yymmdd'] < self.date2)
#                & (self.data['FDLLOutsideTemp']>-10)
                &(np.invert(np.isnan(self.data['FDLLOutsideTemp']) ))
                &(np.invert(np.isinf(self.data['FDLLOutsideTemp']) ))
                &(np.invert(np.isnan(self.data['LLHumidity']) ))
                &(np.invert(np.isinf(self.data['LLHumidity']) ))
                )

        #then, get the data we are interested in:
        self.tempLL = self.getdata('FDLLOutsideTemp',basiccond)
        self.pressureLL = self.getdata('LLPressure',basiccond)
        self.tempelec = self.getdata('TempElec',basiccond)
        self.humLL = self.getdata('LLHumidity',basiccond)
        self.radio = self.getdata('Anode',basiccond,0)
        self.varradio = self.getdata('VarianceAnode',basiccond,0)
        self.time = self.getdata('Time',basiccond)
        self.hhmm = self.getdata('hhmm',basiccond)
        self.solarV = self.getdata('SolarPanelV',basiccond)
        self.solarI = self.getdata('SolarPanelI',basiccond)
        self.UBV = self.getdata('UB24V',basiccond)
        self.powerinput = self.getdata('PowerInput',basiccond)
        self.date = utils.nptstamptodatetime(self.time)
        if len(self.period) != 0:
            print 'get the period'
            self.period = self.period[basiccond]

    def freedata(self):
        self.data = []

    def correctsaturation(self):
        radio = self.radio
        radioshift = np.roll(radio,-1)
        newradio = np.array([])
        size = len(radio)
        jumplimit = 400
        for i in range(size-1):
            r = radio[i]
            rs = radio[i+1]
            if r-rs > jumplimit:
                newrs = rs + 655
                radio[i+1] = newrs
            elif r-rs < -jumplimit:
                newrs = rs - 655
                radio[i+1] = newrs
            else:
                radio[i+1] = radio[i+1]
        self.radio = radio

    def correctbaseline(self):
        radio = self.radio
        newradio = np.array([])
        limit = 200
        for amp in radio:
            if amp < limit:
                newamp = amp+655
            else:
                newamp = amp
            newradio = np.append(newradio,newamp)
        self.radio = newradio

    def correctUBV(self,params):
        #params: [Vreg (or Vcst), Vshift, Vth, Slope]
        vreg = params[0]
        vshift = params[1]
        vth = params[2]
        slope = params[3]
        radio = np.array(self.radio)
        ubv = np.array(self.UBV)
        vcorrplus = slope*(vreg-24)
        vcorrminus = vcorrplus + slope*(ubv-vth)
        radio[ubv>vth] = radio[ubv>vth] + 511.5*(vcorrplus + vshift)
        radio[ubv<vth] = radio[ubv<vth] + 511.5*(vcorrminus[ubv<vth] + vshift)
        self.radio = radio
        return radio


    def filterHP(self, timeinday):
        if (len(self.time) %2 == 1):
            self.tempLL = self.tempLL[:-1]
            self.pressureLL = self.pressureLL[:-1]
            self.humLL = self.humLL[:-1]
            self.radio = self.radio[:-1]
            self.varradio = self.varradio[:-1]
            self.radioc = self.radioc[:-1]
            self.time = self.time[:-1]
            self.date = self.date[:-1]
            self.hhmm = self.hhmm[:-1]
            self.period = self.period[:-1]
            self.solarV = self.solarV[:-1]
            self.solarI = self.solarI[:-1]    
            self.UBV = self.UBV[:-1]
            self.powerinput = self.powerinput[:-1]
#        print 'self.time[1] = ' , self.time[1] , 'self.time[0] = ', self.time[0]
        deltat = self.time[1] - self.time[0]
        tcut = timeinday*24
        fcut = 1./(tcut*3600)
        self.tempLL = utils.highpass(self.tempLL,1./(deltat),2,fcut)
        self.pressureLL = utils.highpass(self.pressureLL,1./(deltat),2,fcut)
        self.radio = utils.highpass(self.radio,1./(deltat),2,fcut)
#        self.UBV = utils.highpass(self.UBV,1./(deltat),2,fcut)
        self.powerinput = utils.highpass(self.powerinput,1./(deltat),2,fcut)
        
    def slidingwindow(self, nbins):
        self.radio = utils.slidingwindow(self.radio,nbins,'valid')
        self.powerinput = utils.slidingwindow(self.powerinput,nbins,'valid')
        size = len(self.radio)
        self.tempLL = utils.getxforconv(self.tempLL,size)
        self.date = utils.getxforconv(self.date,size)
        self.time = utils.getxforconv(self.time,size)
        self.humLL = utils.getxforconv(self.humLL,size)
        self.hhmm = utils.getxforconv(self.hhmm,size)
        self.period =utils.getxforconv(self.period,size)

#     def tempcorrectioncond(self,fitorder, cond):
#         temp = self.tempLL[cond]
#         radio = self.radio[cond]
#         fit = np.polyfit(temp,radio,fitorder)
#         self.fit = fit
#         pfit = np.poly1d(fit)
#         self.radioc = self.radio-pfit(self.tempLL)

    def tempcorrection(self,fitorder,corr=None):
        if (corr == None or corr==False):
            fit = np.polyfit(self.tempLL,self.radio,fitorder)
            self.fit = fit
            pfit = np.poly1d(fit)
            self.radioc = self.radio-pfit(self.tempLL)

        elif corr==True:
            fit = np.polyfit(self.tempLL,self.radioc,fitorder)
            self.fit = fit
            pfit = np.poly1d(fit)
            self.radioc = self.radioc-pfit(self.tempLL)


    def tempcorrectionwithfit(self,fit):
        pfit = np.poly1d(fit)
        self.fit = fit
        self.radioc = self.radio-pfit(self.tempLL)

        
    def getdata(self, x, cond=None, inarray=None):
        if inarray is not None:
            xdata = self.data[x][:,inarray]
        else:
            xdata = self.data[x]
        if cond is not None:
            xdata = xdata[cond]
        return xdata


    def getXYdata(self,x,y,date1=None,date2=None):
        if date1 is None and date2 is None:
            xdata = self.data[x]
            ydata = self.data[y]
        elif date2 is None or date1 is None:
            date1 = utils.datestringtodate(date1)
            xdata = self.data[x][(self.data.Time > date1)]
            ydata = self.data[y][(self.data.Time > date1)]
        else:
            date1 = utils.datestringtodate(date1)
            date2 = utils.datestringtodate(date2)
            print type(date1)
            date1ts = utils.datettotimestamp(date1)
            date2ts = utils.datettotimestamp(date2)
            xdata = self.data[x][(self.data.Time > date1ts) & (self.data.Time < date2ts)]
            ydata = self.data[y][(self.data.Time > date1ts) & (self.data.Time < date2ts)]

        return [xdata,ydata]

    def getnewdatawithid(self, station):
        newd = Dataset()
        newd.data = self.data[self.data.LSId == station]
        newd.date1 = self.date1
        newd.date2 = self.date2
        if len(newd.data) == 0:
            print ' apparently this Id is not in the tree '
            print ' will return an empty data set'
        return newd

    def getnewdatasetcond(self, cond):
        #newd = self.getnewdataset()
        newd = Dataset()
        newd.tempLL = self.tempLL[cond]
        newd.pressureLL = self.pressureLL[cond]
        newd.radio = self.radio[cond]
        newd.varradio = self.varradio[cond]
        if (len(self.radioc) !=0):
            newd.radioc = self.radioc[cond]
        newd.tempelec = self.tempelec[cond]
        newd.humLL = self.humLL[cond]
        newd.presesureLL = self.pressureLL[cond]
        newd.time = self.time[cond]
        newd.hhmm = self.hhmm[cond]
        newd.solarV = self.solarV[cond]
        newd.solarI = self.solarI[cond]
        newd.UBV = self.UBV[cond]
        newd.powerinput = self.powerinput[cond]
        newd.date = self.date[cond]
        if (len(self.period) !=0):
            newd.period = self.period[cond]
        return newd
        

    def getnewdataset(self, t1=None, t2=None):
        newd = Dataset()
        timecond = np.where(self.time !=0)
        if t1 !=None and t2!=None:
            timecond = np.where(
                (self.time > t1) 
                & (self.time  < t2)
                )
        elif t1 != None:
            print 'ici'
            timecond = np.where(self.time > t1)
        elif t2 != None:
            print 'icic'
            timecond = np.where(self.time  < t2)
            #        else
        #then, get the data we are interested in:
        newd.tempLL = self.tempLL[timecond]
        newd.pressureLL = self.pressureLL[timecond]  
        newd.radio = self.radio[timecond]
        newd.varradio = self.varradio[timecond]
        if (len(self.radioc) !=0):
            newd.radioc = self.radioc[timecond]
        newd.tempelec = self.tempelec[timecond]
        newd.humLL = self.humLL[timecond]
        newd.time = self.time[timecond]
        newd.hhmm = self.hhmm[timecond]
        newd.solarV = self.solarV[timecond]
        newd.solarI = self.solarI[timecond]
        newd.date = self.date[timecond]
        newd.UBV = self.UBV[timecond]
        newd.powerinput = self.powerinput[timecond]
        if (len(self.period) !=0):
            newd.period = self.period[timecond]
        return newd


    def getfakedailybaseline(self,typeofbl=None,maxadc=None,corr=None):
        #        self.tempcorrection(1)
        dataarray = []
        nrofday = int((self.time[-1] - self.time[0])/(24*60*60))
        firstday = utils.tstamptodatetime(self.time[0])
        lastday = utils.tstamptodatetime(self.time[-1])
        firstday = datetime.datetime(firstday.year,firstday.month,firstday.day,0,0,0)
        hour = 0
        minute = 0      
        for d in range(0,nrofday,1):
            t0 = firstday + datetime.timedelta(days=d)
            t1 = firstday + datetime.timedelta(days=d+1)
            t0 = utils.datettotimestamp(t0)
            t1 = utils.datettotimestamp(t1)
            daydata = self.getnewdataset(t0,t1)
            if len(daydata.radio) > 215: 
                dataarray.append(daydata)

        sizeofdat = len(dataarray)
        maxlen = 216 # maximum nr of point in one day (1 point every 400s)
        # 1) create a average spectrum
        maxfreq = np.fft.rfftfreq(maxlen,400)
        specarray = np.ndarray(shape= (len(dataarray), len(maxfreq) ) )
        phasearray =np.ndarray(shape= (len(dataarray), len(maxfreq) ) )
        count = 0        
        for dat in dataarray:
            if corr == None or corr==True:
                radiodata = dat.radioc
            else:
                radiodata = dat.radio
            fft = np.fft.rfft(radiodata)
            spec = np.absolute(fft)
            phase = np.angle(fft)
            freq = np.fft.rfftfreq(len(radiodata),400)
            spec = np.interp(maxfreq,freq,spec)
            phase = np.interp(maxfreq,freq,phase)
            specarray[count] = spec
            phasearray[count] = phase
            count +=1
        meanspec = np.mean(specarray,axis=0)
        stdspec = np.std(specarray,axis=0)
        # we set the number of possible baseline as the possible combination 
        # of spectrum/phase from real data
        nrfake = sizeofdat*sizeofdat    
        # we have implemented two possibilities to produce fake baseline:
        # - either we draw a random spectrum around the mean spectrum and then 
        # assume one of the phase of the data
        # -  either we combine one spectrum with one phase.
        # 2) draw random spectrum:
        radtemp = []
        fakebls = []
        random = False
        combine = False
        if typeofbl.lower()=='random':
            for i in range(nrfake):
                fakespec = np.array([])
                for m,s in zip(meanspec,stdspec):
                    specpoint = np.random.normal(m,s)
                    fakespec = np.append(fakespec,specpoint)
                phaseindex = int(np.random.uniform(0,sizeofdat))
                fakephase = phasearray[phaseindex]
                fakefft = fakespec*np.exp(1J*fakephase)
                fakebl = np.fft.irfft(fakefft)
                fakebls.append(fakebl)
        # 3) draw a particular phase
        elif typeofbl.lower() == 'combine' or typeofbl==None:
            for i in range(sizeofdat):
                specindex = i
                for j in range(i, i+sizeofdat):
                    phaseindex = j%sizeofdat
                    fakespec = specarray[specindex]
                    fakephase = phasearray[phaseindex]
                    fakefft = fakespec*np.exp(1J*fakephase)
                    fakebl = np.fft.irfft(fakefft)
                    fakebls.append(fakebl)
        if maxadc == None:
            maxadc = 0
        #produce fake signal
        time = np.linspace(0,24,len(fakebls[0]))
        a = maxadc
        b = 17
        c = 2
        sig = utils.gauss(time,a,b,c)
        count = 0 
        for bl in fakebls:
            rt = radiotemp.Radiotemp()
            bl = bl - sig
            fakebls[count] = bl
            rt.radio  = bl
            rt.corr  = corr
            radtemp.append(rt)
            count +=1
        return radtemp
        
    def gettruedailybaseline(self,maxadc=None,corr=None):
        dataarray = []
        nrofday = int((self.time[-1] - self.time[0])/(24*60*60))
        firstday = utils.tstamptodatetime(self.time[0])
        lastday = utils.tstamptodatetime(self.time[-1])
        firstday = datetime.datetime(firstday.year,firstday.month,firstday.day,0,0,0)
        hour = 0
        minute = 0      
        length = 216
        for d in range(0,nrofday,1):
            t0 = firstday + datetime.timedelta(days=d)
            t1 = firstday + datetime.timedelta(days=d+1)
            t0 = utils.datettotimestamp(t0)
            t1 = utils.datettotimestamp(t1)
            daydata = self.getnewdataset(t0,t1)
            if len(daydata.radio)== length and len(daydata.tempLL)== length: 
                dataarray.append(daydata)
        count = 0
        truebls = []
        truetemps = []
        radtemp = []
        #produce fake signal
        if maxadc!=None:
            time = np.linspace(0,24,length)
            a = maxadc
            b = 17
            c = 2
            sig = utils.gauss(time,a,b,c)
        for dat in dataarray:
            rt = radiotemp.Radiotemp()
            if corr == None or corr==True:
                radiodata = dat.radioc
#                print 'dat.radio = ', dat.radio
                rt.corr = True
            else:
                radiodata = dat.radio
                rt.corr = False
            truebl = radiodata
            if maxadc!=None:
                truebl = truebl - sig
            truetemps.append(dat.tempLL)
            rt.radio = truebl
            rt.temp = dat.tempLL
            rt.date = dat.date
#            print 'rt.date = ' , rt.date
            radtemp.append(rt)
        return radtemp


    def getfulldays(self):
        nrofday = int((self.time[-1] - self.time[0])/(24*60*60))
        firstday = utils.tstamptodatetime(self.time[0])
        lastday = utils.tstamptodatetime(self.time[-1])
        firstday = datetime.datetime(firstday.year,firstday.month,firstday.day,0,0,0)
        hour = 0
        minute = 0      
        length = 216
        days = []
        for d in range(0,nrofday,1):
            t0date = firstday + datetime.timedelta(days=d)
            t1date = firstday + datetime.timedelta(days=d+1)
            
            t0 = utils.datettotimestamp(t0date)
            t1 = utils.datettotimestamp(t1date)
            daydata = self.getnewdataset(t0,t1)
            if len(daydata.radio)== length and len(daydata.tempLL)== length: 
                days.append(t0date)
        return days


    def getradiowithtimearray(self, t):
        radio  = np.interp(t,self.time,self.radio)
        temp  = np.interp(t,self.time,self.tempLL)
        return [temp,radio]


    def fitwithexpo0(self, a, sigma,mu,c):
        x = utils.hhmmtohour(self.hhmm[(self.radioc<60) & (self.radioc > -60)])
        hhmm = self.hhmm[(self.radioc<60) & (self.radioc > -60)]
        radio = self.radioc[ (self.radioc<60) & (self.radioc > -60)]
        try:
            popt, pcov = curve_fit(utils.expofunc0,hhmm,radio, p0=[a,sigma,mu,c])
        except RuntimeError:
            print("Error - curve_fit failed")
        #        popt, pcov = curve_fit(utils.expofunc,x,radio, p0=[a,sigma,mu,b,c])
        return [hhmm,popt]

    def fitwithexpotwolinear(self, a, sigma,mu,p00,p01,p10,p11):
        x = utils.hhmmtohour(self.hhmm[(self.radioc<60) & (self.radioc > -60)])
        hhmm = self.hhmm[(self.radioc<60) & (self.radioc > -60)]
        radio = self.radioc[ (self.radioc<60) & (self.radioc > -60)]
        try:
            popt, pcov = curve_fit(utils.expofunctwolinear,hhmm,radio, p0=[a,sigma,mu,p00,p01,p10,p11])
        except RuntimeError:
            print("Error - curve_fit failed")
        #        popt, pcov = curve_fit(utils.expofunc,x,radio, p0=[a,sigma,mu,b,c])
        return [hhmm,popt]

    def fitwithexpo(self, a, sigma,mu,b,c):
        x = utils.hhmmtohour(self.hhmm[(self.radioc<60) & (self.radioc > -60)])
        hhmm = self.hhmm[(self.radioc<60) & (self.radioc > -60)]
        radio = self.radioc[ (self.radioc<60) & (self.radioc > -60)]
        try:
            popt, pcov = curve_fit(utils.expofunc,hhmm,radio, p0=[a,sigma,mu,b,c])
        except RuntimeError:
            print("Error - curve_fit failed")
        #        popt, pcov = curve_fit(utils.expofunc,x,radio, p0=[a,sigma,mu,b,c])
        return [hhmm,popt]

    def fitwithexpo2(self, a, sigma,mu,b,c,d):
        x = utils.hhmmtohour(self.hhmm[(self.radioc<60) & (self.radioc > -60)])
        x= x*100
        hhmm = self.hhmm[(self.radioc<60) & (self.radioc > -60)]
        radio = self.radioc[ (self.radioc<60) & (self.radioc > -60)]
        try:
#            popt, pcov = curve_fit(utils.expofunc2,hhmm,radio, p0=[a,sigma,mu,b,c,d])
            popt, pcov = curve_fit(utils.expofunc2,x,radio, p0=[a,sigma,mu,b,c,d])
            return [hhmm,popt,True]
        except (RuntimeError, TypeError, NameError):
#        except RunTimeError:
#        except Improper input:
            print("Error - curve_fit failed")
            return [hhmm,[],False]
#        except TypeError:
#        except Improper input:
#            print("Error - curve_fit failed")
#            return [hhmm,[],False]

#        popt, pcov = curve_fit(utils.expofunc2,x,radio, p0=[a,sigma,mu,b,c,d])
#        return [hhmm,popt]
