# The sample ExoPlex run provided here was produced with the following requirements:
  # *all* Fe must be distributed to core - i.e., NO Fe in mantle minerals
  # Na is also excluded from mantle phases
  # Ca, Al, Si, Mg, and O are included in mantle phases
  # (NOTE: ExoPlex file name shows, in order: Ca/Mg, Si/Mg, Al/Mg, Fe/Mg)

# Most previous studies have assumed strictly Mg silicates, so we'll start there re: material parameters.
# However, the data provided may also include (relatively minor contributions from) Ca and/or Al phases of
# akimotoite, orthopyroxene, clinopyroxene, perovskite, and post-perovskite. Stixrude database abbreviations
# for these additional phases are in comments next to their Mg analogs.

# 'Cp' coefficients below were approximated by querying ENKI Portal's (not yet public) ThermoEngine interface for
# Stixrude database, to find Cp across T range(1000,3000,5) at P=750000, then applying scipy's curve_fit function to the
# Berman database's equation structure.

# Code snippet is below. Cps is a numpy array containing Cp for each mineral at each T queried.
# Column 0 is T; all others are mineral Cps.

# from scipy.optimize import curve_fit

# def berman_coeffs(Tp,k0,k1,k2,k3,k4,k5,k6):
#     return k0 + k1*Tp**(-0.5) + k2*Tp**(-2) + k3*Tp**(-3) + k4*Tp**(-1) + k5*Tp + k6*Tp**2
#     return molcp

# xdata = Cps[:,0]
# ncols = len(Cps[0])
# for i in range(1,ncols):
#     ydata = Cps[:,i]
#     popt,pcov = curve_fit(berman_coeffs, xdata, ydata)
#     plt.plot(xdata, berman_coeffs(xdata, *popt))
#     print(list(mins)[i-1],popt)




#Note that for the majority of these, P =  1bar, T = range(295, 2300, 5). Quartz has P=7.0GPa to maintain stability.

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


minerals['C2/c']['name'] ='C2/c'
minerals['C2/c']['stix'] ='hpcEn'
minerals['C2/c']['MW'] =  200.7774
minerals['C2/c']['Cp'] = [2.27829771e+02, 1.52787017e+03, -7.61645263e+06, 1.00000000e+00, -3.58789979e+04, 4.85420496e-03, 2.70149911e-08]
minerals['C2/c']['alpha'] = [-2.61577690e-03, 3.53969405e-04, 1.04893756e-02, -8.65247622e-02, 1.18526795e-11, 2.81205315e+02, 5.86425467e-05]

minerals['Wus']['name'] ='Periclase'
minerals['Wus']['stix'] ='Per'
minerals['Wus']['MW'] =  40.3044
minerals['Wus']['Cp'] = [4.77942048e+01, 1.72585091e+02, -1.63216436e+06, 1.00000000e+00, -4.47609499e+03, 7.10543684e-04, 5.08410582e-08]
minerals['Wus']['alpha'] = [-1.44801427e-02, 1.95818108e-03, 5.64031450e-02, -4.28913402e-01, 5.19454687e-11, 1.44140702e+03, 3.26476299e-04]

minerals['Pv']['name'] ='Mg-Perovskite'
minerals['Pv']['stix'] ='MgPrv' #others: AlPrv, ca-pv included
minerals['Pv']['MW'] =  100.3887
minerals['Pv']['Cp'] = [1.16712336e+02, 6.86749361e+02, -4.41730797e+06, 1.00000000e+00, -1.82392373e+04, 2.58132717e-03, 2.31526744e-07]
minerals['Pv']['alpha'] = [-1.38469328e-03, 1.84920943e-04, 5.74338354e-03, -5.15642676e-02, 6.01570064e-12, 1.82434593e+02, 3.03435551e-05]

minerals['an']['name'] ='Anorthite'
minerals['an']['stix'] ='An'
minerals['an']['MW'] =  278.20928
minerals['an']['Cp'] = [3.12970082e+02, 7.72690260e+02, -8.02001680e+06, 1.00000000e+00, -1.76746147e+04, 1.36761594e-03, -6.89745353e-08]
minerals['an']['alpha'] = [-9.34009553e-05, 1.14531855e-05, 5.32822973e-04, -7.11315798e-03, 2.40775934e-13, 1.85829201e+01, 1.65498628e-06]

minerals['O']['name'] ='Forsterite'
minerals['O']['stix'] ='Fo'
minerals['O']['MW'] =  140.6931
minerals['O']['Cp'] = [1.63794657e+02, 7.60006106e+02, -5.14171950e+06, 1.00000000e+00, -1.77597894e+04, 1.69537901e-03, -1.18970042e-08]
minerals['O']['alpha'] = [-6.04458372e-03, 8.19818391e-04, 2.36261823e-02, -1.82109450e-01, 2.35199643e-11, 6.02596070e+02, 1.36718160e-04]

minerals['Wad']['name'] ='Mg-Wadsleyite'
minerals['Wad']['stix'] ='MgWds'
minerals['Wad']['MW'] =  140.6931
minerals['Wad']['Cp'] = [1.62608641e+02, 8.74374078e+02, -5.45585633e+06, 1.00000000e+00, -2.09787785e+04, 2.26948100e-03, 4.90548272e-08]
minerals['Wad']['alpha'] = [-6.94031620e-03, 9.40252015e-04, 2.71192718e-02, -2.08918944e-01, 2.64001620e-11, 6.99514807e+02, 1.56783099e-04]


minerals['Ring']['name'] ='Ringwoodite'
minerals['Ring']['stix'] ='MgRwd'
minerals['Ring']['MW'] =  140.6931
minerals['Ring']['Cp'] = [1.61832368e+02, 9.13111471e+02, -5.37106987e+06, 1.00000000e+00, -2.16086958e+04, 1.98722871e-03, 1.09870168e-08]
minerals['Ring']['alpha'] = [-3.50493058e-03, 4.74685331e-04, 1.38079188e-02, -1.09078175e-01, 1.40961876e-11, 3.66484899e+02, 7.90239255e-05]


minerals['Opx']['name'] ='Orthopyroxene/En'
minerals['Opx']['stix'] ='En' #others: MgTs, oDi. k0,k1,k2,k3=1332.63600,-9604.704,-18164480.000,2233202400.
minerals['Opx']['MW'] =  200.7774
minerals['Opx']['Cp'] = [2.37841407e+02, 7.95686393e+02, -6.69645970e+06, 1.00000000e+00, -1.82513768e+04, 1.19949468e-03, -6.89957310e-08]
minerals['Opx']['alpha'] = [-2.67039901e-01, 3.55042125e-02, 1.05013416e+00, -8.07058690e+00, 7.81239473e-10, 2.89208074e+04, 5.87291138e-03]

minerals['Cpx']['name'] ='Clinopyroxene/cEn'
minerals['Cpx']['stix'] ='cEn' #others: CaTs, Di
minerals['Cpx']['MW'] =  200.7774
minerals['Cpx']['Cp'] = [2.33547673e+02, 1.11335531e+03, -7.33924996e+06, 1.00000000e+00, -2.59553590e+04, 2.73309035e-03, -2.26478533e-08]
minerals['Cpx']['alpha'] = [-9.41030115e-03, 1.27559276e-03, 3.66676806e-02, -2.80013398e-01, 3.55277002e-11, 9.32810911e+02, 2.12826234e-04]

minerals['Aki']['name'] ='Akimotoite'
minerals['Aki']['stix'] ='MgAki' #others: AlAki
minerals['Aki']['MW'] =  100.3887
minerals['Aki']['Cp'] = [1.14175783e+02, 7.63880116e+02, -3.93415025e+06, 1.00000000e+00, -1.82609813e+04, 1.60804810e-03, 2.16096009e-08]
minerals['Aki']['alpha'] = [-9.26069286e-03, 1.25273272e-03, 3.61124291e-02, -2.76755648e-01, 3.40386513e-11, 9.45128496e+02, 2.08901657e-04]

minerals['Gt_maj']['name'] ='Majoritic Garnet'
minerals['Gt_maj']['stix'] ='Maj'
minerals['Gt_maj']['MW'] =  401.5548
minerals['Gt_maj']['Cp'] = [4.68059234e+02, 2.17334178e+03, -1.47345964e+07, 1.00000000e+00, -5.09106023e+04, 5.79519759e-03, -3.87668241e-09]
minerals['Gt_maj']['alpha'] = [-8.50252056e-04, 1.13512741e-04, 3.57311272e-03, -3.25100293e-02, 3.84877785e-12, 1.04944226e+02, 1.85348164e-05]

minerals['Ppv']['name'] ='Post-perovskite'
minerals['Ppv']['stix'] ='MgPpv' #others: AlPpv
minerals['Ppv']['MW'] =  100.3887
minerals['Ppv']['Cp'] = [1.30959324e+02, -3.22278673e+01, -6.00013296e+06, 1.00000000e+00, -6.41307679e+03, 1.23564196e-03, 7.12128547e-07]
minerals['Ppv']['alpha'] = [-3.55817160e-03, 4.81760018e-04, 1.41660048e-02, -1.15190024e-01, 1.52841297e-11, 3.85063056e+02, 8.00169513e-05]

minerals['CF']['name'] ='Mg-Clinoferrosilite'
minerals['CF']['stix'] ='hpcEn'
minerals['CF']['MW'] =  200.7774
minerals['CF']['Cp'] = [2.27829771e+02, 1.52787017e+03, -7.61645263e+06, 1.00000000e+00, -3.58789979e+04, 4.85420496e-03, 2.70149911e-08]
minerals['CF']['alpha'] = [-2.61577690e-03, 3.53969405e-04, 1.04893756e-02, -8.65247622e-02, 1.18526795e-11, 2.81205315e+02, 5.86425467e-05]

minerals['st']['name'] ='Stishovite'
minerals['st']['stix'] ='Sti'
minerals['st']['MW'] =  60.0843
minerals['st']['Cp'] = [1.77337094e+02, -5.40968943e+03, -1.50571616e+07, 1.00000000e+00, 9.12219000e+04, -1.47817928e-02, 1.51752147e-06]
minerals['st']['alpha'] = [-2.57120795e-03, 3.46825997e-04, 1.02395506e-02, -8.45796039e-02, 1.04177113e-11, 3.20394308e+02, 5.76283687e-05]

minerals['q']['name'] ='Quartz'
minerals['q']['stix'] ='Qz'
minerals['q']['MW'] =  60.0843
minerals['q']['Cp'] = [7.27929443e+01, 1.41834097e+02, -1.77003676e+06, 1.00000000e+00, -3.28665849e+03, 3.93611628e-04, -4.97032464e-09]
minerals['q']['alpha'] =  [-1.10333968e+00, 1.30654569e-01, 4.84809451e+00, -4.51987051e+01, 1.58594476e-09, -9.00568725e+00,  2.03775761e-02]

minerals['ca-pv']['name'] ='Ca-perovskite'
minerals['ca-pv']['stix'] ='CaPrv'
minerals['ca-pv']['MW'] =  116.16369999999999
minerals['ca-pv']['Cp'] = [1.24349825e+02, 2.48712015e+02, -5.01188090e+06, 1.00000000e+00, -9.70025110e+03, 2.39294589e-03, 4.55523051e-07]
minerals['ca-pv']['alpha'] = [-4.46111814e-04, 5.71553141e-05, 2.14734460e-03, -2.43125032e-02, 2.49897612e-12, 7.52788329e+01, 8.83620546e-06]

minerals['cfs']['name'] ='Mg-CaFerrite'
minerals['cfs']['stix'] ='MgCf'
minerals['cfs']['MW'] =  142.26568
minerals['cfs']['Cp'] = [1.63133719e+02, 8.64879302e+02, -5.60792628e+06, 1.00000000e+00, -2.12353915e+04, 3.33196186e-03, 1.27048791e-07]
minerals['cfs']['alpha'] = [-7.94480612e-04, 1.05126732e-04, 3.42630290e-03, -3.28753535e-02, 3.53780452e-12, 1.09080706e+02, 1.70274898e-05]

minerals['coe']['name'] ='Coesite'
minerals['coe']['stix'] ='Coe'
minerals['coe']['MW'] =  60.0843
minerals['coe']['Cp'] = [7.08868885e+01, 2.69991265e+02, -2.05729299e+06, 1.00000000e+00, -6.18231418e+03, 4.56727462e-04, -2.45946238e-08]
minerals['coe']['alpha'] = [-9.08521507e-05, 1.07959229e-05, 5.17026250e-04, -7.18119832e-03, 1.74087375e-13, 2.49111482e+01, 1.56157181e-06]

minerals['ky']['name'] ='Kyanite'
minerals['ky']['stix'] ='Ky'
minerals['ky']['MW'] =  162.04558
minerals['ky']['Cp'] = [1.76920777e+02, 1.57335327e+03, -5.94551970e+06, 1.00000000e+00, -3.65052249e+04, 3.64528932e-03, -7.58649057e-08]
minerals['ky']['alpha'] = [-4.09598255e-04, 5.17459467e-05, 1.95265974e-03, -2.26307637e-02, 1.67843065e-12, 9.06108873e+01, 8.07660465e-06]

minerals['seif']['name'] ='Seifertite'
minerals['seif']['stix'] ='Seif'
minerals['seif']['MW'] =  60.0843
minerals['seif']['Cp'] = [6.34587322e+01, 8.49086478e+02, -2.15902907e+06, 1.00000000e+00, -2.07714244e+04, 1.46163518e-03, 4.26150099e-08]
minerals['seif']['alpha'] = [-2.58952716e-03, 3.49439141e-04, 1.03061983e-02, -8.52119364e-02, 1.05181772e-11, 3.27518651e+02, 5.80840856e-05]

