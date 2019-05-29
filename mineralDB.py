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




# Note that P =  75 Gpa! T = range(295, 2300, 5). Quartz has P=7.0GPa to maintain stability.

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
minerals['C2/c']['Cp'] = [1025.6591486257257, 11818.859998847893, -63858153.6438232, 9537641908.433262, -205319.0534887796, 0.045437236874565544, -2.111618540270154e-06]
minerals['C2/c']['alpha'] = [-7.524002910306186e-05, 8.10449193704788e-06, 0.0004836297309670938, -0.00830681672925755, 7.748481275327078e-14, 49.36108832996416, 1.1053498606618068e-06]

minerals['Wus']['name'] ='Periclase'
minerals['Wus']['stix'] ='Per'
minerals['Wus']['MW'] =  40.3044
minerals['Wus']['Cp'] = [404.1398695719276, 44869.36043558041, 29990967.637979966, 1.0, -794433.0149694802, 0.14134504004314302, -1.1068623955655294e-05]
minerals['Wus']['alpha'] = [-9.748451697348211e-05, 1.0857594825196061e-05, 0.0005954520777660942, -0.009532372315951535, 1.3784959755845453e-13, 49.97322932651001, 1.5210190930104927e-06]

minerals['Pv']['name'] ='Mg-Perovskite'
minerals['Pv']['stix'] ='MgPrv' #others: AlPrv, ca-pv included
minerals['Pv']['MW'] =  100.3887
minerals['Pv']['Cp'] = [899.3598454956191, 19517.500144376067, -52912592.33872181, 9855795887.22111, -359219.63672784844, 0.06898234289952984, -1.9503652440605634e-06]
minerals['Pv']['alpha'] = [-9.133936829520618e-05, 9.541290348745544e-06, 0.0006169712891312184, -0.011436710571413709, 1.3358871308035114e-13, 78.71995579467153, 1.2619521673784443e-06]

minerals['an']['name'] ='Anorthite'
minerals['an']['stix'] ='An'
minerals['an']['MW'] =  278.20928
minerals['an']['Cp'] = [1168.3508978168977, -645.542922365265, -56683776.32172916, 5759983464.917573, 26282.707235418482, 0.001934799996160609, -1.404496123879176e-07]
minerals['an']['alpha'] = [-1.9292084877503827e-05, 2.2989302995853087e-06, 0.00011787399209816958, -0.0017631959684502553, 1.93345958807218e-14, 5.967395164005401, 3.3099674732601117e-07]

minerals['O']['name'] ='Forsterite'
minerals['O']['stix'] ='Fo'
minerals['O']['MW'] =  140.6931
minerals['O']['Cp'] = [1142.613289120178, 4997.614795172775, -68659557.8520133, 8636718934.346085, -75450.88927443516, 0.020139386466534997, -1.0865270895455417e-06]
minerals['O']['alpha'] = [-5.197625793496891e-05, 5.79353681272302e-06, 0.00031951659319754403, -0.005088606242304254, 6.330813381822515e-14, 25.232587126150463, 8.08430089196178e-07]

minerals['Wad']['name'] ='Mg-Wadsleyite'
minerals['Wad']['stix'] ='MgWds'
minerals['Wad']['MW'] =  140.6931
minerals['Wad']['Cp'] = [1075.3552534361452, 8911.543813314505, -66927874.55737444, 9328683125.041737, -150220.69439370383, 0.03305958376620849, -1.4862358369130071e-06]
minerals['Wad']['alpha'] = [-6.826734114616092e-05, 7.485956146880651e-06, 0.0004249473699989218, -0.00700046985048004, 8.91725544824556e-14, 39.01269479748892, 1.0363679319427074e-06]

minerals['Ring']['name'] ='Ringwoodite'
minerals['Ring']['stix'] ='MgRwd'
minerals['Ring']['MW'] =  140.6931
minerals['Ring']['Cp'] = [1075.01438635049, 8897.416097804395, -66939873.79552747, 9314924361.648382, -149440.86884689835, 0.03033770315724921, -1.6931541829694122e-06]
minerals['Ring']['alpha'] = [-5.588778864316124e-05, 6.110643071810583e-06, 0.00034843470568069807, -0.005741989493801437, 7.376528578430402e-14, 31.930218401203568, 8.428594238629225e-07]

minerals['Opx']['name'] ='Orthopyroxene/En'
minerals['Opx']['stix'] ='En' #others: MgTs, oDi. k0,k1,k2,k3=1332.63600,-9604.704,-18164480.000,2233202400.
minerals['Opx']['MW'] =  200.7774
minerals['Opx']['Cp'] = [1214.158594385337, 1013.2694281324275, -65779608.77063802, 7279953182.908563, -1731.6555519755727, 0.005891467670361558, -5.503872077918571e-07]
minerals['Opx']['alpha'] = [-1.0702829405423318e-05, 1.2121990178240237e-06, 6.575624029969593e-05, -0.001009420386786074, 1.56774617595697e-14, 4.093781110879472, 1.6786385610939955e-07]

minerals['Cpx']['name'] ='Clinopyroxene/cEn'
minerals['Cpx']['stix'] ='cEn' #others: CaTs, Di
minerals['Cpx']['MW'] =  200.7774
minerals['Cpx']['Cp'] = [449.8353103364512, 42429.07452343824, 26423252.8371617, 1.0, -747245.7014922032, 0.12727792268538396, -1.1496652726011655e-05]
minerals['Cpx']['alpha'] = [-5.501726047795653e-05, 6.1431377572672424e-06, 0.00033846221081888593, -0.005405861373938266, 6.267813891751682e-14, 27.022867147457866, 8.600911751833464e-07]

minerals['Aki']['name'] ='Akimotoite'
minerals['Aki']['stix'] ='MgAki' #others: AlAki
minerals['Aki']['MW'] =  100.3887
minerals['Aki']['Cp'] = [998.6645275719708, 13418.052560732784, -62349202.070575505, 9713208005.277187, -236364.4297871703, 0.04159006234315155, -2.4297753797098807e-06]
minerals['Aki']['alpha'] = [-5.404710345393571e-05, 5.80776415746769e-06, 0.0003440809050149698, -0.005894435254963551, 7.83021792523853e-14, 36.08084231453707, 7.906137663323928e-07]

minerals['Gt_maj']['name'] ='Majoritic Garnet'
minerals['Gt_maj']['stix'] ='Maj'
minerals['Gt_maj']['MW'] =  401.5548
minerals['Gt_maj']['Cp'] = [1141.7685759486499, 5127.450312161002, -68661630.24266943, 8666443057.619387, -77998.42038041295, 0.02286024906167644, -1.0434690058795255e-06]
minerals['Gt_maj']['alpha'] = [-6.146137149815503e-05, 6.874900451835035e-06, 0.0003776047464721723, -0.006021140547371027, 6.999691815797051e-14, 29.99014130702474, 9.641809923501479e-07]

minerals['Ppv']['name'] ='Post-perovskite'
minerals['Ppv']['stix'] ='MgPpv' #others: AlPpv
minerals['Ppv']['MW'] =  100.3887
minerals['Ppv']['Cp'] = [895.3302719765878, 20543.894102909813, -50681148.52582984, 10062237091.437113, -390918.8032447353, 0.07567844089396092, 1.0530071790140821e-06]
minerals['Ppv']['alpha'] = [-0.0001191823449498798, 1.2530023017625056e-05, 0.0007974829257734519, -0.014852473945544686, 2.5240932078223114e-13, 106.42025966958569, 1.6634831990833895e-06]

minerals['CF']['name'] ='Mg-Clinoferrosilite'
minerals['CF']['stix'] ='hpcEn'
minerals['CF']['MW'] =  200.7774
minerals['CF']['Cp'] = [1025.6591486257257, 11818.859998847893, -63858153.6438232, 9537641908.433262, -205319.0534887796, 0.045437236874565544, -2.111618540270154e-06]
minerals['CF']['alpha'] = [-7.524002910306186e-05, 8.10449193704788e-06, 0.0004836297309670938, -0.00830681672925755, 7.748481275327078e-14, 49.36108832996416, 1.1053498606618068e-06]

minerals['st']['name'] ='Stishovite'
minerals['st']['stix'] ='Sti'
minerals['st']['MW'] =  60.0843
minerals['st']['Cp'] = [914.8902078654878, 20731.212506242013, -57224613.25249578, 12227880713.553852, -419641.82764810056, 0.04674418538561956, -1.465283744873362e-06]
minerals['st']['alpha'] = [-0.0013218113958139801, 0.00016012889159311175, 0.005778318370049411, -0.055403067800412625, 2.040019587669454e-12, 294.4001222037483, 2.513291048030152e-05]

minerals['q']['name'] ='Quartz'
minerals['q']['stix'] ='Qz'
minerals['q']['MW'] =  60.0843
minerals['q']['Cp'] = [1260.6552743762795, -1325.9479149251931, -56949165.60169142, 5522224005.254361, 38474.73291542475, 0.002576804634352612, 1.255335082664611e-07]
minerals['q']['alpha'] = [1.573380089083486e-06, -2.9885055482167204e-07, -4.225874046274854e-06, 2.12830045446341e-05, 1.1804525328325757e-14, -0.039344823103704074, -6.705094118445598e-08]

minerals['ca-pv']['name'] ='Ca-perovskite'
minerals['ca-pv']['stix'] ='CaPrv'
minerals['ca-pv']['MW'] =  116.16369999999999
minerals['ca-pv']['Cp'] = [877.4557029144429, 11269.028594503166, -54251091.539402746, 8560729336.840203, -205536.52450404645, 0.05275059800009518, 7.461726916994674e-07]
minerals['ca-pv']['alpha'] = [-0.00011772684720937997, 1.2613628615810695e-05, 0.0007564427651031725, -0.01315127462742632, 1.900349897194678e-13, 82.55681799415672, 1.7083550550759398e-06]

minerals['cfs']['name'] ='Mg-CaFerrite'
minerals['cfs']['stix'] ='MgCf'
minerals['cfs']['MW'] =  142.26568
minerals['cfs']['Cp'] = [1043.1745839243615, 10048.350186517271, -64928080.39937413, 9351435113.578207, -173159.7025387903, 0.04363623581802041, -1.2416303128882374e-06]
minerals['cfs']['alpha'] = [-9.356112139477609e-05, 1.0209920759169271e-05, 0.000589595304718903, -0.00989479691341634, 1.1221236636980947e-13, 57.0042428886838, 1.4093480328399917e-06]

minerals['coe']['name'] ='Coesite'
minerals['coe']['stix'] ='Coe'
minerals['coe']['MW'] =  60.0843
minerals['coe']['Cp'] = [585.9100643338426, 35077.094784737914, 17192141.733436067, 1.0, -612995.4713988506, 0.10283896194268785, -9.993717932096531e-06]
minerals['coe']['alpha'] = [-2.362409802201591e-05, 2.727625364070623e-06, 0.00014328256839896603, -0.002209175625594023, 2.149929233862014e-14, 9.543746048703294, 3.9117097161230223e-07]

minerals['ky']['name'] ='Kyanite'
minerals['ky']['stix'] ='Ky'
minerals['ky']['MW'] =  162.04558
minerals['ky']['Cp'] = [912.0212952876877, 17743.58109049697, -55376522.70437774, 9683789332.242794, -318920.1741378065, 0.05410270509796997, -3.6780172354479827e-06]
minerals['ky']['alpha'] = [-6.256496142973576e-05, 6.60485387107993e-06, 0.0004207719307798186, -0.0077075357253996325, 4.978715367723395e-14, 51.13766676277885, 8.887607397133211e-07]

minerals['seif']['name'] ='Seifertite'
minerals['seif']['stix'] ='Seif'
minerals['seif']['MW'] =  60.0843
minerals['seif']['Cp'] = [503.78603941078865, 42921.05103118808, -5762887.685564069, 8222242786.756776, -818110.691226836, 0.10835852033113826, -7.259333760842503e-06]
minerals['seif']['alpha'] = [-2.250941618188277e-05, 1.6699725699524992e-06, 0.0002497153648994858, -0.006692479848909264, 6.023492058393578e-14, 62.23424054591963, 1.0580399857488087e-07]
