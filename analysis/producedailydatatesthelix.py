###########################################
## produce python type files             ##
## with radio, temp, humidity simulation ##
## out of the root file produced with    ##
## Corinne's code.                       ##
###########################################
import numpy as np
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
parser.add_argument("deltatheta", type=int, nargs='?',default=0, help="deviation from central value")
parser.add_argument("deltaphi", type=int, nargs='?',default=0, help="deviation from central value")
parser.add_argument("--filename", type=str, nargs='?',default='alldata.root', help="root file")
parser.add_argument("--stname", type=str, nargs='?',default='popey', help="station name")
parser.add_argument("--fluxorigin", type=str, nargs='?',default='nobeyama', help="canadian or nobeyama")
args = parser.parse_args()
delt = args.deltatheta
delp = args.deltaphi
fname = args.filename
stname = args.stname
flux = args.fluxorigin
stid = constant.Helixstationidbyname[stname]

print ' delt = ' ,delt ,' and delp = ', delp 
datafolder = '/t'+str(delt) + 'p'+str(delp)+'/'
############################
## load the root file     ##
############################
file = constant.datafolder + 'GIGADuck/' + fname
data = dataset.Dataset(file)
data.loaddata()
data = data.getnewdatawithid(stid)
data.selectleafs()
data.correctbaseline()

if flux == 'nobeyama':
    simulationfolder = constant.simulationfolder2
    outfolder = constant.resultfolder2 
elif flux == 'canadian':
    simulationfolder = constant.simulationfolder
    outfolder = constant.resultfolder 
elif flux == 'test':
    simulationfolder = constant.simulationfoldertest 
    outfolder = constant.resultfoldertest + datafolder
else:
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    print 'choose a possible type of flux'
    print 'now exiting the script'
    print '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
    sys.exit()    

years = [2016,2017] 
for y in years:
    for d in range(1,constant.doy[y],1):
        sdatelim1 = utils.doytodate(y,d,hour=3)
        sdatelim2 = utils.doytodate(y,d+1,hour=3)
        cond = np.where((data.date > sdatelim1) & (data.date < sdatelim2) )
        daydataset = data.getnewdatasetcond(cond)
        if (len(daydataset.date) == 0):
            continue
        oritheta = constant.Helixantennageo[stname][0]
        oriphi = constant.Helixantennageo[stname][1]


        if flux=='test':
            simfile = simulationfolder + str(stid) + '/exptemp_' + stname + '_' +str(delt) + '_' + str(delp) + '_' + str(y) + '_' + str(d) + '.txt'
            [h,m,temp] = utils.readtempfile(simfile)
            sec = 3600*h + 60*m
            ddata = daydata.Daydata('simpledata', y, d, daydataset.date, daydataset.radio, daydataset.tempLL, daydataset.humLL, sec, temp)

        ######################################
        ## load the sun simulation file     ##
        ######################################
        ## name example: exptemp_popey_2015_354.txt 
        else:
            simfile = simulationfolder + str(stid) + '/exptemp_' + stname + '_' + str(y) + '_' + str(d) + '.txt'
            [h,m,temp] = utils.readtempfile(simfile)
            sec = 3600*h + 60*m
            ddata = daydata.Daydata('simpledata', y, d, daydataset.date, daydataset.radio, daydataset.tempLL, daydataset.humLL, sec, temp)


        #######################################
        ## save in pkl format the daily data ##
        #######################################
        outf = outfolder + str(stid)
        outfilename = outf + '/data_' + str(y) + '_' + str(d) + '.pkl'
        out = open(outfilename,'wb')
        pickle.dump(ddata,out)
        out.close()
