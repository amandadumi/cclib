# Copyright (c) 2024, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test unrestrictied single point logfiles in cclib"""

import numpy
from skip import skipForParser


class GenericSPunTest:
    """Generic unrestricted single point unittest"""

    def testnatom(self, data) -> None:
        """Is the number of atoms equal to 20?"""
        assert data._ccCollection._parsed_data[0].natom == 20

    def testatomnos(self, data) -> None:
        """Are the atomnos correct?"""
        assert numpy.all(
            [
                numpy.issubdtype(atomno, numpy.signedinteger)
                for atomno in data._ccCollection._parsed_data[0].atomnos
            ]
        )
        assert data._ccCollection._parsed_data[0].atomnos.shape == (20,)
        assert (
            sum(data._ccCollection._parsed_data[0].atomnos == 6)
            + sum(data._ccCollection._parsed_data[0].atomnos == 1)
            == 20
        )

    @skipForParser("ADF", "???")
    @skipForParser(
        "DALTON",
        "DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now",
    )
    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser("Gaussian", "The parser is still being developed for version 2")
    @skipForParser("Jaguar", "???")
    @skipForParser("Molcas", "Length is zero for some reason")
    @skipForParser("Molpro", "???")
    @skipForParser("Psi4", "The parser is still being developed for version 2")
    @skipForParser("Turbomole", "???")
    def testatomcharges(self, data) -> None:
        """Are atomic charges consistent with natom?"""
        for atomcharge_type in data._ccCollection._parsed_data[0].atomcharges:
            charges = data._ccCollection._parsed_data[0].atomcharges[atomcharge_type]
            natom = data._ccCollection._parsed_data[0].natom
            assert (
                len(charges) == natom
            ), f"len(atomcharges['{atomcharge_type}']) = {len(charges)}, natom = {natom}"

    @skipForParser("ADF", "???")
    @skipForParser(
        "DALTON",
        "DALTON has a very low accuracy for the printed values of all populations (2 decimals rounded in a weird way), so let it slide for now",
    )
    @skipForParser("FChk", "The parser is still being developed so we skip this test")
    @skipForParser("Gaussian", "The parser is still being developed for version 2")
    @skipForParser("Jaguar", "???")
    @skipForParser("Molcas", "Length is zero for some reason")
    @skipForParser("Molpro", "???")
    @skipForParser("Psi4", "The parser is still being developed for version 2")
    @skipForParser("Turbomole", "???")
    def testatomcharges_mulliken(self, data) -> None:
        """Do Mulliken atomic charges sum to positive one?"""
        charges = data.atomcharges["mulliken"]
        assert abs(sum(charges) - 1.0) < 1.0e-2

    def testatomcoords(self, data) -> None:
        """Are the dimensions of atomcoords 1 x natom x 3?"""
        assert data._ccCollection._parsed_data[0].atomcoords.shape == (
            1,
            data._ccCollection._parsed_data[0].natom,
            3,
        )

    @skipForParser("Jaguar", "Data file does not contain enough information")
    def testdimmocoeffs(self, data) -> None:
        """Are the dimensions of mocoeffs equal to 2 x nmo x nbasis?"""
        if hasattr(data, "mocoeffs"):
            assert isinstance(data.mocoeffs, list)
            assert len(data.mocoeffs) == 2
            assert data._ccCollection._parsed_data[0].mocoeffs[0].shape == (
                data._ccCollection._parsed_data[0].nmo,
                data._ccCollection._parsed_data[0].nbasis,
            )
            assert data._ccCollection._parsed_data[0].mocoeffs[1].shape == (
                data._ccCollection._parsed_data[0].nmo,
                data._ccCollection._parsed_data[0].nbasis,
            )

    @skipForParser("Jaguar", "Data file does not contain enough information")
    @skipForParser("DALTON", "mocoeffs not implemented yet")
    def testfornoormo(self, data) -> None:
        """Do we have NOs or MOs?"""
        assert hasattr(data._ccCollection._parsed_data[0], "nocoeffs") or hasattr(
            data._ccCollection._parsed_data[0], "mocoeffs"
        )

    @skipForParser("Gaussian", "The parser is still being developed for version 2")
    @skipForParser("Psi4", "The parser is still being developed for version 2")
    def testdimnoccnos(self, data) -> None:
        """Is the length of nooccnos equal to nmo?"""
        if hasattr(data._ccCollection._parsed_data[0], "nooccnos"):
            assert isinstance(data._ccCollection._parsed_data[0].nooccnos, numpy.ndarray)
            # FIXME
            assert data._ccCollection._parsed_data[0].nooccnos.shape in [
                (data._ccCollection._parsed_data[0].nmo,),
                (2, data._ccCollection._parsed_data[0].nmo),
            ]

    @skipForParser("Gaussian", "The parser is still being developed for version 2")
    @skipForParser("Psi4", "The parser is still being developed for version 2")
    def testdimnocoeffs(self, data) -> None:
        """Are the dimensions of nocoeffs equal to 2 x nmo x nmo?"""
        if hasattr(data._ccCollection._parsed_data[0], "nocoeffs"):
            assert isinstance(data._ccCollection._parsed_data[0].nocoeffs, numpy.ndarray)
            assert data.nocoeffs.shape == (
                2,
                data._ccCollection._parsed_data[0].nmo,
                data._ccCollection._parsed_data[0].nmo,
            )

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    def testcharge_and_mult(self, data) -> None:
        """Are the charge and multiplicity correct?"""
        assert data._ccCollection._parsed_data[0].charge == 1
        assert data._ccCollection._parsed_data[0].mult == 2

    @skipForParser("Gaussian", "The parser is still being developed for version 2")
    @skipForParser("Psi4", "The parser is still being developed for version 2")
    def testhomos(self, data) -> None:
        """Are the homos correct?"""
        msg = f"{numpy.array_repr(data._ccCollection._parsed_data[0].homos)} != array([34,33],'i')"
        numpy.testing.assert_array_equal(
            data._ccCollection._parsed_data[0].homos, numpy.array([34, 33], "i"), msg
        )

    def testmoenergies(self, data) -> None:
        """Are the dims of the moenergies equals to 2 x nmo?"""
        if hasattr(data, "moenergies"):
            assert len(data._ccCollection._parsed_data[0].moenergies) == 2
            assert (
                len(data._ccCollection._parsed_data[0].moenergies[0])
                == data._ccCollection._parsed_data[0].nmo
            )
            assert (
                len(data._ccCollection._parsed_data[0].moenergies[1])
                == data._ccCollection._parsed_data[0].nmo
            )

    @skipForParser("FChk", "Fchk files do not have a section for symmetry")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Molpro", "?")
    @skipForParser("ORCA", "ORCA has no support for symmetry yet")
    @skipForParser("Psi4", "The parser is still being developed for version 2")
    def testmosyms(self, data) -> None:
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape = (
            len(data._ccCollection._parsed_data[0].mosyms),
            len(data._ccCollection._parsed_data[0].mosyms[0]),
        )
        assert shape == (2, data._ccCollection._parsed_data[0].nmo)


class GenericROSPTest(GenericSPunTest):
    """Customized restricted open-shell single point unittest"""

    @skipForParser("DALTON", "mocoeffs not implemented yet")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testdimmocoeffs(self, data) -> None:
        """Are the dimensions of mocoeffs equal to 1 x nmo x nbasis?"""
        assert isinstance(data._ccCollection._parsed_data[0].mocoeffs, list)
        assert len(data._ccCollection._parsed_data[0].mocoeffs) == 1
        assert data._ccCollection._parsed_data[0].mocoeffs[0].shape == (
            data._ccCollection._parsed_data[0].nmo,
            data._ccCollection._parsed_data[0].nbasis,
        )

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testhomos(self, data) -> None:
        """Are the HOMO indices equal to 34 and 33 (one more alpha electron
        than beta electron)?
        """
        msg = (
            f"{numpy.array_repr(data._ccCollection._parsed_data[0].homos)} != array([34, 33], 'i')"
        )
        numpy.testing.assert_array_equal(
            data._ccCollection._parsed_data[0].homos, numpy.array([34, 33], "i"), msg
        )

    @skipForParser("QChem", "prints 2 sets of different MO energies?")
    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testmoenergies(self, data) -> None:
        """Are the dims of the moenergies equals to 1 x nmo?"""
        assert len(data._ccCollection._parsed_data[0].moenergies) == 1
        assert (
            len(data._ccCollection._parsed_data[0].moenergies[0])
            == data._ccCollection._parsed_data[0].nmo
        )

    @skipForParser("Molcas", "The parser is still being developed so we skip this test")
    @skipForParser("Turbomole", "The parser is still being developed so we skip this test")
    def testmosyms(self, data) -> None:
        """Are the dims of the mosyms equals to 1 x nmo?"""
        shape = (
            len(data._ccCollection._parsed_data[0].mosyms),
            len(data._ccCollection._parsed_data[0].mosyms[0]),
        )
        assert shape == (1, data._ccCollection._parsed_data[0].nmo)


class GamessUK70SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testdimmocoeffs(self, data) -> None:
        """Are the dimensions of mocoeffs equal to 2 x (homos+6) x nbasis?"""

        assert isinstance(data._ccCollection._parsed_data[0].mocoeffs, list)
        assert len(data._ccCollection._parsed_data[0].mocoeffs) == 2

        # This is only an issue in version 7.0 (and before?), since in the version 8.0
        # logfile all eigenvectors are happily printed.
        shape_alpha = (
            data._ccCollection._parsed_data[0].homos[0] + 6,
            data._ccCollection._parsed_data[0].nbasis,
        )
        shape_beta = (
            data._ccCollection._parsed_data[0].homos[1] + 6,
            data._ccCollection._parsed_data[0].nbasis,
        )
        assert data._ccCollection._parsed_data[0].mocoeffs[0].shape == shape_alpha
        assert data._ccCollection._parsed_data[0].mocoeffs[1].shape == shape_beta

    def testnooccnos(self, data) -> None:
        """Are natural orbital occupation numbers the right size?"""
        assert data._ccCollection._parsed_data[0].nooccnos.shape == (
            data._ccCollection._parsed_data[0].nmo,
        )


class GamessUK80SPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testnooccnos(self, data) -> None:
        """Are natural orbital occupation numbers the right size?"""
        assert data.nooccnos.shape == (data.nmo,)


class GaussianSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    @skipForParser("Gaussian", "The parser is still being developed for version 2")
    def testatomspins(self, data) -> None:
        """Are atomic spins from Mulliken population analysis consistent with
        natom and sum to one (doublet)?
        """
        spins = data.atomspins["mulliken"]
        assert len(spins) == data.natom
        assert abs(sum(spins) - 1.0) < 0.001


class JaguarSPunTest(GenericSPunTest):
    """Customized unrestricted single point unittest"""

    def testmoenergies(self, data) -> None:
        """Are the dims of the moenergies equal to 2 x homos+11?"""
        assert len(data.moenergies) == 2
        assert len(data.moenergies[0]) == data.homos[0] + 11
        assert len(data.moenergies[1]) == data.homos[1] + 11

    def testmosyms(self, data) -> None:
        """Are the dims of the mosyms equals to 2 x nmo?"""
        shape0 = (len(data.mosyms), len(data.mosyms[0]))
        shape1 = (len(data.mosyms), len(data.mosyms[1]))
        assert shape0 == (2, data.homos[0] + 11)
        assert shape1 == (2, data.homos[1] + 11)
