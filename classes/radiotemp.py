import utils
import numpy as np
import datetime 

class Radiotemp:
    def __init__(self):
        self.radio = np.array([])
        self.temp = np.array([])
        self.fit = []
        self.corr = False
        self.date = None
        self.signalmax = None
        
