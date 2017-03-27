listname = 'fevmarch2016.txt'
f = open(listname,'w')
year = 2016
dstart = 32
dstop = 90
for d in range(dstart,dstop+1,1):
    f.write(str(year) + ' ' + str(d) + '\n')

