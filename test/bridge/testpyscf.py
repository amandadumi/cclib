# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest

import numpy as np
import sys
from cclib.bridge import cclib2pyscf
from cclib.parser.utils import find_package


class PyscfTest(unittest.TestCase):
    self.data, self.logfile = getdatafile(
            "Gaussian", "basicGaussian16", ["dvb_un_sp.fchk"]
        )
    datadir = os.path.abspath(
            os.path.join(os.path.dirname(__file__), "..", "..", "data")
        )

    """Tests for the cclib2pyscf bridge in cclib."""

    def setUp(self):
        super(PyscfTest, self).setUp()        
        if not find_package('pyscf'):
            raise ImportError('Must install pyscf to run this test')

    def test_makepyscf(self):
        import pyscf
        from pyscf import scf

        atomnos = np.array([1, 8, 1], "i")
        atomcoords = np.array([[-1, 1, 0], [0, 0, 0], [1, 1, 0]], "f")
        pyscfmol = cclib2pyscf.makepyscf(atomcoords, atomnos)
        pyscfmol.basis = "6-31G**"
        pyscfmol.cart = True
        pyscfmol.verbose = 0
        pyscfmol.build()

        mhf = pyscfmol.HF(conv_tol=1e-6)
        en = mhf.kernel()
        ref = -75.824754602
        assert abs(en - ref) < 1.0e-6

    def test_makepyscf_from_molden(self):
        pyscfmol, pyscf_data = cclib2pyscf.makepyscf_from_molden(self.data)
        assert np.allclose(pyscf_data['mo_energy'], self.data['mo_energy'])
        # check first MO coefficient 
        assert pyscfmol['mocoeffs'][0] == self.data['mo_energy'][0][0]
        # check a random middle MO coefficient 
        assert pyscfmol['mocoeffs'][10] == self.data['mo_energy'][0][10]  
        # test behavior of parser if datatype isn't present.


if __name__ == "__main__":
    unittest.main()
    PyscfTest.test_makepyscf_from_molden()
