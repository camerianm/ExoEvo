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


def unchanging(params):
	print('Planet mass:\t', Pe(params['Mp']), 'kg')
	print('Core mass:\t', Pe(params['Mc']), 'kg')
	print('Planet radius:\t', Pe(params['Rp']), 'm')
	print('Core radius:\t', Pe(params['Rc']), 'm')
	print('Mantle depth:\t', Pe(params['d']), 'm')
	print('Surface area:\t', Pe(params['Sa']), 'm2')
	print('Density:\t', Pf(params['pm']), 'kg m-3')
	print('Surf gravity:\t', Pf(params['g']), 'm s-2')
	print('Visc constant:\t', Pf(params['c1']))
	print('Activ. energy:\t', Pe(params['Ev']), 'J/mol')
	print('Visc prefactor:\t', Pe(params['visc0']), "Pa s")
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
