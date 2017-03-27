import utils
import numpy as np
import datetime 
from scipy.optimize import curve_fit
import radiotemp
class Daydata:
    def __init__(self, name = '', year=0, doy=0, timedata=[], data=[], temp=[], humidity=[], timesim=[], sim=[]):
        self.name = name
        self.timedata = timedata
        self.data = data
        self.tempLL = temp
        self.humidity =  humidity
        self.timesim = timesim
        self.sim = sim
        
