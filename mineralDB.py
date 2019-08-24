# The sample ExoPlex run provided here was produced with the following requirements:
  # *all* Fe must be distributed to core - i.e., NO Fe in mantle minerals
  # Na is also excluded from mantle phases
  # Ca, Al, Si, Mg, and O are included in mantle phases
  # (NOTE: ExoPlex file name shows, in order: Ca/Mg, Si/Mg, Al/Mg, Fe/Mg)

# Most previous studies have assumed strictly Mg silicates, so we'll start there re: material parameters.
# However, ExoPlex outputs may also include (relatively minor contributions from) Ca and/or Al phases of
# akimotoite, orthopyroxene, clinopyroxene, perovskite, and post-perovskite. Stixrude database abbreviations
# for these additional phases are in comments next to their Mg analogs.

#Unless otherwise indicated, thermal conductivities from: DOI:10.1029/JB076i005p01278

#Prospective additions
#SI from http://www.pnas.org/cgi/doi/10.1073/pnas.1110594108 
#For wad/ring/o: P/T expression for thermal conductivity on page 10 of DOI:10.1016/j.pepi.2004.03.005
#For periclase: thermal conductivity as a function of P/T - DOI: 10.1073/pnas.0907194107 
#For O (Fo93): 10.1016/j.pepi.2003.10.010
#For Mg perovskite: P/T/XFe dependencies in follow-up studies from 2017AGUFMMR31A0426T
#Also for perovskite: XFe dependencies in: 10.1073/pnas.1718557115 
minerals = {}
minerals['C2/c'] = {}
minerals['Wus'] = {}
minerals['Pv'] = {}
minerals['an'] = {}
minerals['O'] = {}
minerals['Wad'] = {}
minerals['Ring'] = {}
minerals['Opx'] = {}
minerals['Cpx'] = {}
minerals['Aki'] = {}
minerals['Gt_maj'] = {}
minerals['Ppv'] = {}
minerals['CF'] = {}
minerals['st'] = {}
minerals['q'] = {}
minerals['ca-pv'] = {}
minerals['cfs'] = {}
minerals['coe'] = {}
minerals['ky'] = {}
minerals['seif'] = {}
minerals['Sp'] = {}

kdefault = 5.00
alphadefault=1.0e-6
Ev_default=300.0e3

minerals['C2/c']['name'] ='C2/c'
minerals['C2/c']['stix'] ='hpcEn' #More info at DOI: 10.2138/am-2019-6740
minerals['C2/c']['k'] = kdefault #NEED
minerals['C2/c']['water'] = 714.0/(1e6) #DOI: 10.1007/s00410-002-0365-6

minerals['Wus']['name'] ='Periclase'
minerals['Wus']['stix'] ='Per'
minerals['Wus']['k'] = kdefault #NEED
minerals['Wus']['water'] = 0.0075/100.

minerals['Pv']['name'] ='Mg-Perovskite'
minerals['Pv']['stix'] ='MgPrv' #others: AlPrv, ca-pv included
minerals['Pv']['k'] = kdefault #NEED
minerals['Pv']['water'] = 0.001/100.

minerals['an']['name'] ='Anorthite'
minerals['an']['stix'] ='An'
minerals['an']['k'] = 1.71544 
minerals['an']['Ev'] = 648.0e3 #DOI: 10.1029/2000JB900223
minerals['an']['water'] = 0.051/100.

minerals['O']['name'] ='Forsterite'
minerals['O']['stix'] ='Fo'
minerals['O']['k'] = 5.10448 #5.858
minerals['O']['Ev'] = 261.0e3 #DOI: 10.1029/2007JB005100
minerals['O']['water'] = 0.12/100.

minerals['Wad']['name'] ='Mg-Wadsleyite'
minerals['Wad']['stix'] ='MgWds'
minerals['Wad']['k'] = kdefault #NEED
minerals['Wad']['Ev'] = 261.0e3 #DOI: 10.1029/2007JB005100
minerals['Wad']['water'] = 2.4/100.

minerals['Ring']['name'] ='Ringwoodite'
minerals['Ring']['stix'] ='MgRwd'
minerals['Ring']['k'] = kdefault #NEED
minerals['Ring']['Ev'] = 261.0e3 #DOI: 10.1029/2007JB005100
minerals['Ring']['water'] = 2.5/100.

minerals['Opx']['name'] ='Orthopyroxene/En'
minerals['Opx']['stix'] ='En' #others: MgTs, oDi. k0,k1,k2,k3=1332.63600,-9604.704,-18164480.000,2233202400.
minerals['Opx']['k'] = 4.3932 #NC sample from Horai 1971 p.1988
minerals['Opx']['Ev'] = 420.0e3 #DOI: 10.1002/jgrb.50284  Alt DOI: 10.1002/2017JB014400
minerals['Opx']['water'] = 0.15/100.  #pure enstatite is lower (199e-6) in DOI: 10.1016/j.epsl.2005.04.022 

minerals['Cpx']['name'] ='Clinopyroxene/cEn'
minerals['Cpx']['stix'] ='cEn' #others: CaTs, Di
minerals['Cpx']['k'] = kdefault #NEED
minerals['Cpx']['Ev'] = 560.0e3 #DOI: 10.1029/2001JB000333
minerals['Cpx']['water'] = 0.08/100.

minerals['Aki']['name'] ='Akimotoite'
minerals['Aki']['stix'] ='MgAki' #others: AlAki
minerals['Aki']['k'] = kdefault #NEED
minerals['Aki']['water'] = 443.0/(1e6) #DOI: 10.1016/S0012-821X(00)00244-2

minerals['Gt_maj']['name'] ='Majoritic Garnet'
minerals['Gt_maj']['stix'] ='Maj'
minerals['Gt_maj']['k'] = 9.77 #This is room temp value though... DOI: 10.1016/S0012-821X(03)00630-7
minerals['Gt_maj']['water'] = 0.0675/100.

minerals['Ppv']['name'] ='Post-perovskite'
minerals['Ppv']['stix'] ='MgPpv' #others: AlPpv
minerals['Ppv']['k'] = kdefault #NEED
minerals['Ppv']['water'] = 0.001/100.

minerals['CF']['name'] ='Mg-Clinoferrosilite'
minerals['CF']['stix'] ='MgCf'
minerals['CF']['k'] = kdefault #NEED

minerals['st']['name'] ='Stishovite'
minerals['st']['stix'] ='Sti'
minerals['st']['k'] = 7.82 #Median of mantle geotherm values reported in DOI: 10.1029/2011JB009119
minerals['st']['water'] = 100./(1.0e6) #DOI: 10.1029/2002JB002053

minerals['q']['name'] ='Quartz'
minerals['q']['stix'] ='Qz'
minerals['q']['k'] = 7.686008

minerals['ca-pv']['name'] ='Ca-perovskite'
minerals['ca-pv']['stix'] ='CaPrv'
minerals['ca-pv']['k'] = kdefault #NEED
minerals['ca-pv']['water'] = 0.004/100.

minerals['cfs']['name'] ='hpcFs'
minerals['cfs']['stix'] ='hpcFs'
minerals['cfs']['k'] = kdefault #NEED

minerals['coe']['name'] ='Coesite'
minerals['coe']['stix'] ='Coe'
minerals['coe']['k'] = kdefault #NEED

minerals['ky']['name'] ='Kyanite'
minerals['ky']['stix'] ='Ky'
minerals['ky']['k'] = 14.154472
minerals['ky']['water'] = 100./(1.0e6) #DOI: 10.2138/rmg.2006.62.8

minerals['seif']['name'] ='Seifertite'
minerals['seif']['stix'] ='Seif'
minerals['seif']['k'] = kdefault #NEED

minerals['Sp']['name'] ='Spinel'
minerals['Sp']['stix'] ='Spl'
minerals['Sp']['k'] = 9.47676
minerals['Sp']['water'] = 0.2/100.


minerals['sample_bulk'] = {'C2/c':5.605938492, 'Wus':0.196424301, 'Pv':58.03824705, 'an':0.00, \
               'O':0.249338793, 'Wad':0.072264906, 'Ring':0.028707673, 'Opx':14.88882685, \
               'Cpx':1.099284717, 'Aki':0.0000703828958617, 'Gt_maj':9.763623743, 'Ppv':6.440039009, \
               'CF':0.00, 'st':0.00, 'q':0.00, 'ca-pv':3.617239801, \
               'cfs':0.00, 'coe':0.00, 'ky':0.00, 'seif':0.00}
