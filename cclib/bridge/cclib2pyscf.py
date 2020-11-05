# -*- coding: utf-8 -*-
#
# Copyright (c) 2019, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Bridge for using cclib data in PySCF (https://github.com/pyscf/pyscf)."""

from cclib.parser.utils import find_package, PeriodicTable
import numpy as np

l_sym2num = {"S": 0, "P": 1, "D": 2, "F": 3, "G": 4}


class MissingAttributeError(Exception):
    pass

_found_pyscf = find_package("pyscf")
if _found_pyscf:
    from pyscf import gto, tools


def _check_pyscf(found_pyscf):
    if not found_pyscf:
        raise ImportError("You must install `pyscf` to use this function")


def makepyscf(data, charge=0, mult=1):
    """Create a Pyscf Molecule."""
    _check_pyscf(_found_pyscf)
    # if hasattr(data, "gbasis"):
    #     basis = 
    inputattrs = data.__dict__
    print(inputattrs['atomcoords'])
    print([["{}".format(inputattrs['atomnos'][i]),inputattrs['atomcoords'][i]] for i in range(len(inputattrs['atomcoords']))])
    mol = gto.Mole(
        atom=[
            ["{}".format(data.atomnos[i]), data.atomcoords[-1][i]]
            for i in range(data.natom)
        ],
        atom=[["{}".format(inputattrs['atomnos'][i]),inputattrs['atomcoords'][i]] for i in range(len(inputattrs['atomcoords']))],
        unit="Angstrom",
        charge=charge,
        # basis = basis,
        multiplicity=mult  
    )
    inputattr = data.__dict__
    pt = PeriodicTable()
    if "gbasis" in inputattr:
        basis = {}  # object for internal PySCF format
        uatoms, uatoms_idx = np.unique(
            data.atomnos, return_index=True
        )  # find unique atoms
        for idx, i in enumerate(uatoms_idx):
            curr_atom_basis = data.gbasis[i]
            for jdx, j in enumerate(curr_atom_basis):
                curr_l = j[0]
                curr_e_prim = j[1]
                new_list = [l_sym2num["{}".format(curr_l)]]
                new_list += curr_e_prim
                if not "{}".format(pt.element[uatoms[idx]]) in basis:
                    basis["{}".format(pt.element[uatoms[idx]])] = [new_list]
                else:
                    basis["{}".format(pt.element[uatoms[idx]])].append(new_list)
        mol.basis = basis
    return mol



def makepyscf_with_mo(data):
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
    # temp = tempfile.NamedTemporaryFile(mode="w+t")
    # Check required attributes.
    required_attrs = {"atomcoords", "atomnos"}
    missing = [x for x in required_attrs if not hasattr(data, x)]
    if missing:
        missing = " ".join(missing)
        raise MissingAttributeError(
            "Could not create pyscf molecule due to missing attribute: {}".format(missing)
        )

    pyscf_data = {}
    # # taking advantage of pyscf molden parser as it accounts for basis normalization specifications.
    # mldn_obj = ccwrite(data, outputtype="molden", returnstr=True)
    # temp.write(mldn_obj)
    # temp.seek(0)
    # (
    #     mol,
    #     pyscf_data["mo_energy"],
    #     pyscf_data["mo_coeff"],
    #     pyscf_data["mo_occ"],
    #     pyscf_data["irrep_labels"],
    #     pyscf_data["spins"],
    # ) = tools.molden.load(temp.name)

    #mo energies
    required_attrs = {"moenergies","mo_syms","mo_syms",""}
    # pyscf needs 1-d array, thus flatten if data is of rank 2 array in case 
    # of unrestricted.
    pyscf_data['mo_energy'] = data['moenergies'].flatten()
    # TODO: check on the form of symmetry labels in PySCF
    pyscf_data['irrep_labels'] = data['mo_syms'].flatten()
    #mo_occ
    # cclib doesn't directorly store occupation, but can create
    # ground state occupation from knowing whether it is restricted or not. 
    pyscf_data['mo_occ'] = np.zeros_like(pyscf_data['mo_energy'])
    if np.shape(data['homo'])[0] == 1:
            pyscf_data['mo_occ'][:data['homos'][0]] = 1
            pyscf_data['spins'][:data['homos'][0]] = 'Alpha'
    elif np.shape(data['homo'])[0] == 2:
            pyscf_data['mo_occ'][data['nmo']:data['nmo']+data['homos'][1]] = 1 
            pyscf_data['spins'][data['nmo']:data['nmo']+data['homos'][1]] = 'Beta'

    #mo_coeffs

    temp.close()
    return mol, pyscf_data


>>>>>>> 4aa2d381 (intermediate progress)
del find_package
