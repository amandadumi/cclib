import copy

import numpy

from cclib.parser.utils import convertor
from cclib.parser.utils import find_package
from cclib.method import Volume



_found_pyscf = find_package("pyscf")
if _found_pyscf:
    from cclib.bridge import cclib2pyquante

    pymol = cclib2pyquante.makepyquante(ccdata)

    # ao_on_grid = pymol.eval_ao(eval_name=GOval_cart,vol)
    # if pymol.mo_coeff == 1: # then ROHF or RHF calc (alpha = beta)
    #     mo_a_on_grid = mo_b_on_grid = ao_on_grid.dot(pymol.mo_coeff[0])
    # else pymol.mo_coeff ==2:
    #     mo_a_on_grid = ao_on_grid.dot(pymol.mo_coeff[0])
    #     mo_b_on_grid = ao_on_grid.dot(pymol.mo_coeff[1])
    # mo_on_grid = mo_a_on_grid + mo_b_on_grid
    # return pygrid

def electrondensity(ccdata, vol, mocoeffslist):
    if len(mocoeffslist) == 2:
        return electrondensity_spin(ccdata, volume, [mocoeffslist[0]]) + electrondensity_spin(
            ccdata, volume, [mocoeffslist[1]]
        )
    else:
        edens = electrondensity_spin(ccdata, volume, [mocoeffslist[0]])
        edens.data *= 2
        return edens


def electrondensity_spin(ccdata, volume, mocoeffslist):
    ao_on_grid = pymol.eval_ao(eval_name=GOval_cart,vol)
    density = copy.copy(volume)
    density.data = numpy.zeros(density.data.shape, "d")
    for mocoeffs in mocoeffslist: # are we looking at alpha or beta orbitals?
        for mocoeff in mocoeffs: #for each mo in the alpha set & beta set maybe?
            mo_on_grid = ao_on_grid.dot(ccdata.mocoeff[0,mocoeffslist])
            wavefn = mo_on_grid.T*mo_on_grid # probably not this easy....
    density.data = wavefn**2
### pyscf: evaluate aos on a grid.
### pull out mo coeffs and dot to get ao on grid. 

