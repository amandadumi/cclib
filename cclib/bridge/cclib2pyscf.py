# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PySCF (https://github.com/pyscf/pyscf)."""

from cclib.parser.utils import find_package
from cclib.io import ccwrite
import tempfile

_found_pyscf = find_package("pyscf")
if _found_pyscf:
    from pyscf import gto, tools


def _check_pyscf(found_pyscf):
    if not found_pyscf:
        raise ImportError("You must install `pyscf` to use this function")


def makepyscf(atomcoords, atomnos, charge=0, mult=1):
    """Create a Pyscf Molecule."""
    _check_pyscf(_found_pyscf)
    mol = gto.Mole(
        atom=[["{}".format(atomnos[i]), atomcoords[i]] for i in range(len(atomcoords))],
        unit="Angstrom",
        charge=charge,
        multiplicity=mult
    )

    return mol


def makepyscf_from_molden(data):
    """
    Creates an molecule object and additional pyscf data formatted for PySCF.
    Parameters:
    ----
    data: cclib data object from parsed output file

    Returns:
    ----
    mol: a pyscf molecule object
        a molecule object that consists of available data in Molden file.
    pyscf_data: dictionary
        a dictionary of additional values returned in an appropriate format for use in PySCF. 
        If available in ccdata, the following objects are available in the dictionary, otherwise, values are None
    """
    _check_pyscf(_found_pyscf)
    temp = tempfile.NamedTemporaryFile(mode="w+t")
    pyscf_data = {}
    # taking advantage of pyscf molden parser as it accounts for basis normalization specifications.
    mldn_obj = ccwrite(data, outputtype="molden", returnstr=True)
    temp.write(mldn_obj)
    temp.seek(0)
    (
        mol,
        pyscf_data["mo_energy"],
        pyscf_data["mo_coeff"],
        pyscf_data["mo_occ"],
        pyscf_data["irrep_labels"],
        pyscf_data["spins"],
    ) = tools.molden.load(temp.name)
    temp.close()
    return mol, pyscf_data


del find_package
