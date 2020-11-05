# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

import unittest
import numpy as np
from cclib.bridge import cclib2pyscf
from cclib.parser.utils import find_package, convertor
from ..test_data import getdatafile


class PyscfTest(unittest.TestCase):
    """Tests for the cclib2pyscf bridge in cclib."""

    def setUp(self):
        super(PyscfTest, self).setUp()        
        if not find_package('pyscf'):
            raise ImportError('Must install pyscf to run this test')
        self.data, self.logfile = getdatafile("GAMESS", "basicGAMESS-US2018", ["water_mp2.out"])
    
    def test_makepyscf(self):
        import pyscf
        from pyscf import scf
        refen = self.data.scfenergies[-1] # value in eVs
        pyscfmol = cclib2pyscf.makepyscf(self.data)
        pyscfmol.basis = "STO-3G"
        pyscfmol.cart = True
        pyscfmol.verbose = 5
        pyscfmol.build()
        mhf = pyscfmol.HF(conv_tol=1e-6)
        en = mhf.kernel()
        assert abs(convertor(en,'hartree','eV')- refen) < 1.0e-6


if __name__ == "__main__":
    unittest.main()
