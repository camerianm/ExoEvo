em = ('An', 'CaPrv', 'cEn', 'Coe', 'En', 'Fa', 'FeAki', 'FeCf', 'FePpv', 'FePrv', 'FeRwd', 'FeWds', 'Fo', 'Fs',
	'Maj', 'hpcEn', 'hpcFs', 'Ky', 'MgAki', 'MgCf', 'MgPpv', 'MgPrv', 'MgRwd',  'MgWds', 'Per',
	'Qz', 'Seif', 'Spl', 'Sti', 'Wus')
m = {}
kdefault=5.0

def generate_endmembers():
	import os
	n=0
	for e in em:
		m[e] = {}
		m[e]['k'] = kdefault
		m[e]['cpgrid'] = os.path.join('CPgrid/', (e + '_cpgrid.csv'))
		m[e]['alphagrid'] = os.path.join('alphagrid/', (e + '_alphagrid.csv'))
		try:
			a=open(m[e]['cpgrid'],'r')
			a=open(m[e]['alphagrid'],'r')
		except:
			print('Please make grid files for %s' % e)
	return m
m=generate_endmembers()




# Thermal conductivities that are not default values
var='k'
m['An'][var]=1.71544
m['Fo'][var]=5.10448
m['En'][var]=4.3932 #NC sample from Horai 1971 p.1988
m['Maj'][var] = 3.175656 #Pyrope value, not Maj, from Horai et al
m['Qz'][var] = 7.686008
m['Ky'][var] = 14.154472
m['Spl'][var] = 9.47676

# ExoPlex Names 
exonames={'C2/c':['hpcEn'], 'Wus':['Per','Wus'], 'Pv':['MgPrv','FePrv'],
	'O':['Fo','Fa'], 'Wad':['MgWds','FeWds'], 'Ring':['MgRwd','FeRwd'],
	'Opx':['En','Fs'], 'Cpx':['cEn'], 'Aki':['MgAki','FeAki'],
	'Gt_maj':['Maj'], 'Ppv':['MgPpv','FePpv'], 'CF':['MgCf','FeCf'],
	'ca-pv':['CaPrv'], 'cfs':['hpcFs'], 'Sp':['Spl'], 'st':['Sti'],
	'q':['Qz'], 'an':['An'], 'coe':['Coe'], 'ky':['Ky'], 'seif':['Seif']}



#FINAL OUTPUT.
minerals=m

