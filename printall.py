# A module to ease printing
# Recommended import style:
#   import printall as prnt
#   from printall import Pe as Pe
#   from printall import Pf as Pf

def Pe(n):
	a=(format(n, '.4e'))
	return a #format(n, '.4e')

	
def Pf(n):
	a=(format(n, '.4f'))
	return a #float(format(n, '.4f'))


def unchanging(params, composition):
    print()
    print('Planet mass:\t', Pf(params['Mpl']), 'Me')
    print('Planet mass:\t', Pe(params['Mp']), 'kg')
    print('Core mass:\t', Pe(params['Mc']), 'kg')
    print('Core mass frac:\t\t', Pf(params['CMF']))
    print()
    print('Planet radius:\t', Pf(params['Rpl']), 'Re')
    print('Planet radius:\t', Pf(params['Rp']/1000), 'km')
    print('Core radius:\t', Pf(params['Rc']/1000), 'km')
    print('Core radius frac:\t', Pf(params['CRF']))
    print()
    print('Mantle depth:\t', Pf(params['d']/1000), 'km')
    print('Mantle volume:\t', Pe(params['Vm']/(1000*1000*1000)), 'km3')
    print('Mantle density:\t', Pf(params['pm']), 'kg m-3')
    print('Surface grav:\t', Pf(params['g']), 'm s-2')
    print('Surface area:\t', Pe(params['Sa']/(1000*1000)), 'km2')
    print()
    print('Core-mantle boundary pressure:\t', Pf(params['Pcmb']), 'Gpa')
    print('Reference pressure:\t', Pf(params['Pref']), 'Gpa')
    print('Radiogenic abundance:\t', Pf(params['Qpl']), 'x Earth')
    print('Starting CMB temp:\t', Pf(params['Tcmb']), 'K')
    print('Viscosity prefactor:\t', Pe(params['visc0']), "Pa s")
    print('Activation energy:\t', Pe(params['Ev']), 'J/mol')
    print()
    
    print('Mantle composition:')
    for i in composition:
        if composition[i]>0:
            print(i, '\t', Pf(100*composition[i]), '%')
    print()
    print()
    return

def options(component):
    print('\n')
    
    if component == "options":
        print('Options for option() function:')
        print('* \'minerals\'')
        print('* \'mode\'')
        
    if component == "mineral":
        print('Options for mineral variable:')
        print('* \'forsterite\'')
        print('* \'fayalite\'')
        print('* \'orthoenstatite\'')
        print('* \'clinoenstatite\'')
        print('* \'periclase\'')
        print('* \'corundum\'')
        print('* \'spinel\'')
        print('* \'diopside\'')
        print('* \'diamond\'')
        print('* \'ca-al pyroxene\'')

    if component == "mode":
        print('Options for mode variable:')
        print('* \'dynamic\' - changes alpha, Cp, and k at each temperature')
        print('* \'static\' - alpha, Cp, and k are values from Foley & Smye 2018, DOI:10.1089/ast.2017.1695')
        #print('* \'dorn\' - Benchmark case: Dorn, et al. 2018, DOI:10.1051/0004-6361/201731513')
        #print('* \'foley\' - Benchmark case: Foley and Smye 2018, DOI:10.1089/ast.2017.1695')
        #print('* \'korenaga\' - Benchmark case: Korenaga 2006, DOI:10.1029/164GM03')

    if component == "outputs":
        print('Column numbers in Evolution output array:')
        print('0  * \'time in Ga\'')
        print('1  * \'temperature in C\'')
        print('2  * \'Rayleigh number\'')
        print('3  * \'radiogenic heat produced\'')
        print('4  * \'heat loss through surface\'')
        print('5  * \'Urey ratio = production/loss\'')
        
    print('\n')
    return
