import numpy as np 
import sys
import os
import fromexo

file='0.05_0.91_0.06_0.47.dat'
startline=1000
cols=[]

radii,wt_local,wt_total=fromexo.rel_shells(file,startline)
status=fromexo.plot_rel_contributions(radii,wt_local)
print(status)
array=fromexo.weightedavg1(file,startline,wt_local,wt_total)
cols=fromexo.read_cols(file)
for i in range(len(list(array[1,:]))):
	print(cols[i], '\t', sum(array[:,i]))
#print(cols)
