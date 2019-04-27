import numpy as np 
import sys
import os
import fromexo
import glob

directory=sys.argv[1]
filenames=sorted(glob.glob(directory+'/*.csv'))
startline=1000
cols=[]
f=open('averages.csv','a+')

for file in filenames:
	a=str(file)
	print(a)
	radii,wt_local,wt_total=fromexo.rel_shells(file,startline)
	status=fromexo.plot_rel_contributions(radii,wt_local)
	
	array=fromexo.weightedavg1(file,startline,wt_local,wt_total)
	cols=fromexo.read_cols(file)
	for i in range(len(list(array[1,:]))):
		f.write(str(file)[len(directory):]+','+str(cols[i])+','+str(sum(array[:,i]))+'\n')
	print(status)

f.close()
print('\n\nPlease see averages.csv for results.\n')