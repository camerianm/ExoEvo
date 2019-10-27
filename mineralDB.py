# The sample ExoPlex run provided here was produced with the following requirements:
# *all* Fe must be distributed to core - i.e., NO Fe in mantle minerals
# Na is also excluded from mantle phases
# Ca, Al, Si, Mg, and O are included in mantle phases
# (NOTE: ExoPlex file name shows, in order: Ca/Mg, Si/Mg, Al/Mg, Fe/Mg)


# Most previous studies have assumed strictly Mg silicates, so we'll start there re: material parameters.
# However, ExoPlex outputs may also include (relatively minor contributions from) Ca and/or Al phases of
# akimotoite, orthopyroxene, clinopyroxene, perovskite, and post-perovskite. Stixrude database abbreviations
# for these additional phases are in comments next to their Mg analogs.

# updated thermal conductivities: table 1 of https://science.sciencemag.org/content/283/5408/1699/tab-figures-data
# Unless otherwise indicated, thermal conductivities from: DOI:10.1029/JB076i005p01278

# Prospective additions
# For wad/ring/o: P/T expression for thermal conductivity on page 10 of DOI:10.1016/j.pepi.2004.03.005
# For O (Fo93): DOI: 10.1016/j.pepi.2003.10.010
# For Mg perovskite: P/T/XFe dependencies in follow-up studies from 2017AGUFMMR31A0426T
# Also for perovskite: XFe dependencies in: DOI: 10.1073/pnas.1718557115
# Anorthite rheology:  DOI: 10.1029/2000JB900223
# for thermal conductivity: doi.org/10.1016/j.epsl.2015.06.050
# For Pv: 10.7  # DOI:10.1038/s41598-017-05523-6 - apparently increases by 5x from 0-100GPa
# For Pv and Wus, k(P,T): 10.1073/pnas.1110594108
# For Wus, k(P,T): doi: 10.1073/pnas.0907194107 
# for Sti: Median of geotherm values reported in DOI: 10.1103/PhysRevB.96.195201 is 16.
# For Pv, radiative portion of k(T): DOI: 10.1126/science.1164609
# For Wus: 33.5 at ambient, per DOI: 10.1029/RF003p0105, p.118-124
# For Wus: 55.2 at ambiet from DOI: 10.1126/science.283.5408.1699
# For Pv, k(T=ambient)~5 DOI: 10.1016/j.epsl.2012.09.002
# For Pv, Ppv, and Wus: k(pm, T) in DOI: 10.1016/j.epsl.2012.09.002 (table 2 and equation 6)

# Necessary changes!!!
# minerals['CF']['name'] = 'MgAl2O4, Ca-ferrite (mfer)' #NEED TO CHANGE TO THIS and get data files!!! Ca-ferrite phase of MgAl2SiO4 spinel!! Nonzero in Exoplex

kdefault = 5.00
alphadefault = 1.0E-6
Ev_default = 300.0e3
grainsize=15. #microns

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

minerals['C2/c']['name'] = 'C2/c'
minerals['C2/c']['stix'] = 'hpcEn'  # More info at DOI: 10.2138/am-2019-6740
minerals['C2/c']['k'] =  4.9  # DOI: 10.1126/science.283.5408.1699
#minerals['C2/c']['Ev'] = 560.0e3  # treats same as cpx, as often done for O/Ring/Wad. DOI: 10.1029/2001JB000333
minerals['C2/c']['water'] = 714.0E-6  # DOI: 10.1007/s00410-002-0365-6

minerals['Wus']['name'] = 'Periclase'
minerals['Wus']['stix'] = 'Per'
minerals['Wus']['k'] = 	20.0  # lower-mantle range from 10.1103/PhysRevLett.104.208501
minerals['Wus']['water'] = 0.0075E-2

minerals['Pv']['name'] = 'Mg-Perovskite'
minerals['Pv']['stix'] = 'MgPrv'  # others: AlPrv, ca-pv included
minerals['Pv']['k'] = 4.7   # DOI: 10.1126/science.283.5408.1699
#minerals['Pv']['Ev'] = 370.0e3  # DOI: 10.1016/j.icarus.2013.03.013
#minerals['Pv']['visc0'] = 2.5e11  # DOI: 10.1016/j.icarus.2013.03.013, normalized to form in Foley & Smye 2017, assuming Tref ~1600
minerals['Pv']['water'] = 0.001E-2

minerals['an']['name'] = 'Anorthite'
minerals['an']['stix'] = 'An'
minerals['an']['k'] = 2.4  # intermediate value from DOI: 10.1029/RF003p0105, p.118-124
#minerals['an']['Ev'] = 467.0e3  # 10.1029/2000JB900223
minerals['an']['water'] = 0.051E-2

minerals['O']['name'] = 'Forsterite'
minerals['O']['stix'] = 'Fo'
minerals['O']['k'] = 5.2  # DOI: 10.1126/science.283.5408.1699
minerals['O']['Ev'] = 300.0e3  # Possibly replace, as this is mantle olivine; DOI: 10.1029/2007JB005100
minerals['O']['visc0'] = 4.0e10 # DOI: 10.1002/jgrb.50284, Mantle olivine: 10**5.25 DOI: 10.1029/2007JB005100
minerals['O']['water'] = 0.12E-2

minerals['Wad']['name'] = 'Mg-Wadsleyite'
minerals['Wad']['stix'] = 'MgWds'
minerals['Wad']['k'] = 8.10  # DOI:10.1016/j.pepi.2004.03.005
# minerals['Wad']['Ev'] = minerals['O']['Ev'] #261.0e3  # O/Wad/Ring treated as same; DOI: 10.1029/2007JB005100
# minerals['Wad']['visc0'] = minerals['O']['visc0']  # DOI: 10.1029/2007JB005100
minerals['Wad']['water'] = 2.4E-2

minerals['Ring']['name'] = 'Ringwoodite'
minerals['Ring']['stix'] = 'MgRwd'
minerals['Ring']['k'] = 10.0   # DOI: 10.1016/S0012-821X(03)00630-7. in 10.1126/science.283.5408.1699 is est. 7.7
# minerals['Ring']['Ev'] = minerals['O']['Ev'] #261.0e3  # O/Wad/Ring treated as same; DOI: 10.1029/2007JB005100
# minerals['Ring']['visc0'] = minerals['O']['visc0']  # DOI: 10.1029/2007JB005100
minerals['Ring']['water'] = 2.5E-2

minerals['Opx']['name'] = 'Orthopyroxene/En'
minerals['Opx']['stix'] = 'En'  # others: MgTs, oDi.
minerals['Opx']['k'] = 4.5  # DOI: 10.1126/science.283.5408.1699
minerals['Opx']['Ev'] = 360.0e3  # DOI: 10.1002/jgrb.50284  Alt DOI, using water content: 10.1002/2017JB014400
minerals['Opx']['water'] = 0.15E-2  # pure enstatite is lower (199E-6) in DOI: 10.1016/j.epsl.2005.04.022

minerals['Cpx']['name'] = 'Clinopyroxene/cEn'
minerals['Cpx']['stix'] = 'cEn'  # others: CaTs, Di
minerals['Cpx']['k'] = 4.9  # DOI: 10.1126/science.283.5408.1699
# minerals['Cpx']['Ev'] = 560.0e3  # DOI: 10.1029/2001JB000333
minerals['Cpx']['water'] = 0.08E-2

minerals['Aki']['name'] = 'Akimotoite'
minerals['Aki']['stix'] = 'MgAki'  # others: AlAki
minerals['Aki']['k'] = minerals['Cpx']['k']  # NEED
minerals['Aki']['water'] = 443.0E-6  # DOI: 10.1016/S0012-821X(00)00244-2

minerals['Gt_maj']['name'] = 'Majoritic Garnet'
minerals['Gt_maj']['stix'] = 'Maj'
minerals['Gt_maj']['k'] = 4.48  # Assuming 3% pyrope. Room temp value. DOI: 10.1016/S0012-821X(03)00630-7
minerals['Gt_maj']['water'] = 0.0675E-2

minerals['Ppv']['name'] = 'Post-perovskite'
minerals['Ppv']['stix'] = 'MgPpv'  # others: AlPpv
minerals['Ppv']['k'] = kdefault  # NEED
minerals['Ppv']['water'] = 0.001E-2

minerals['st']['name'] = 'Stishovite'
minerals['st']['stix'] = 'Sti'
minerals['st']['k'] =   8.6  # DOI: 10.1126/science.283.5408.1699 
minerals['st']['water'] = 100.E-6  # DOI: 10.1029/2002JB002053

minerals['q']['name'] = 'Quartz'
minerals['q']['stix'] = 'Qz'
minerals['q']['k'] = 7.7  # DOI: 10.1029/2005GC001053, DOI: 10.1029/JB076i005p01278

minerals['ca-pv']['name'] = 'Ca-perovskite'
minerals['ca-pv']['stix'] = 'CaPrv'
minerals['ca-pv']['k'] = minerals['Pv']['k']  # NEED - going to assume same as Mg-pv for now.
minerals['ca-pv']['water'] = 0.004E-2

minerals['cfs']['name'] = 'hpcFs'
minerals['cfs']['stix'] = 'hpcFs'
minerals['cfs']['k'] = minerals['C2/c']['k'] - 0.5  # NEED - assuming reduction of 0.5, from 10.1126/science.283.5408.1699

minerals['coe']['name'] = 'Coesite'
minerals['coe']['stix'] = 'Coe'
minerals['coe']['k'] = 8.0  # 4GPa value from DOI: 10.1016/0031-9201(78)90036-5

minerals['ky']['name'] = 'Kyanite'
minerals['ky']['stix'] = 'Ky'
minerals['ky']['k'] = 14.16  #DOI: 10.1016/0012-821X(69)90186-1, DOI: 10.1029/RF003p0105, p.118-124
minerals['ky']['water'] = 100.E-6  # DOI: 10.2138/rmg.2006.62.8

minerals['seif']['name'] = 'Seifertite'
minerals['seif']['stix'] = 'Seif'
minerals['seif']['k'] = 16.  # Geotherm value reported in DOI: 10.1103/PhysRevB.96.195201, w/CaCl2 structure

minerals['Sp']['name'] = 'Spinel'
minerals['Sp']['stix'] = 'Spl'
minerals['Sp']['k'] = 9.48  # Horai 1971 reported in DOI: 10.1029/RF003p0105, p.118-124
minerals['Sp']['water'] = 0.2E-2

minerals['CF']['name'] = 'Mg-Clinoferrosilite' #NEED TO CHANGE THIS!!! Ca-ferrite phase of MgAl2SiO4 !! Nonzero in Exoplex
minerals['CF']['stix'] = 'MgCf'
minerals['CF']['k'] = minerals['Sp']['k']  # NEED TO CHANGE - unsure how phase transition changes this

# A composition that can be imported to test modules.
minerals['sample_bulk'] = {'C2/c': 5.605938492, 'Wus': 0.196424301, 'Pv': 58.03824705, 'an': 0.00,
                           'O': 0.249338793, 'Wad': 0.072264906, 'Ring': 0.028707673, 'Opx': 14.88882685,
                           'Cpx': 1.099284717, 'Aki': 0.0000703828958617, 'Gt_maj': 9.763623743, 'Ppv': 6.440039009,
                           'CF': 0.00, 'st': 0.00, 'q': 0.00, 'ca-pv': 3.617239801,
                           'cfs': 0.00, 'coe': 0.00, 'ky': 0.00, 'seif': 0.00}
