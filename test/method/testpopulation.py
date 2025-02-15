# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the various population analyses (MPA, LPA, CSPA, Bickelhaupt) in cclib"""

import logging
import os
import sys
from typing import Type

from cclib.method import CSPA, LPA, MPA, OPA, Bickelhaupt
from cclib.method.calculationmethod import Method, MissingAttributeError
from cclib.parser import Gaussian

import numpy

sys.path.insert(1, "..")

import pytest

from ..test_data import getdatafile


class PopulationTest:
    """Generic population method tests."""

    methods = (CSPA, LPA, MPA, OPA, Bickelhaupt)

    def calculate(self, method_class: Type[Method]) -> None:
        if not hasattr(self, "data"):
            self.data = parse()
        self.analysis = method_class(self.data)
        self.analysis.logger.setLevel(0)
        self.analysis.calculate()

    def testmissingrequiredattributes(self) -> None:
        """Is an error raised when required attributes are missing?"""
        for missing_attribute in MPA.required_attrs:
            self.data = parse()
            delattr(self.data, missing_attribute)
            for method_class in self.methods:
                with pytest.raises(MissingAttributeError):
                    self.calculate(method_class)

    def testmissingoverlaps(self) -> None:
        """Is an error raised when no overlaps are available?"""
        self.data = parse()
        for overlap_attribute in MPA.overlap_attributes:
            if hasattr(self.data, overlap_attribute):
                delattr(self.data, overlap_attribute)
        for method_class in self.methods:
            if method_class.overlap_attributes:
                with pytest.raises(MissingAttributeError):
                    self.calculate(method_class)


class GaussianMPATest:
    """Mulliken Population Analysis test"""

    @classmethod
    def setup_class(cls) -> None:
        cls.data = parse()
        cls.analysis = MPA(cls.data)
        cls.analysis.logger.setLevel(0)
        cls.analysis.calculate()

    def testsumcharges(self) -> None:
        """Do the Mulliken charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        assert abs(totalpopulation - formalcharge) < 1.0e-3

    def testsumspins(self) -> None:
        """Do the Mulliken spins sum up to the total formal spin?"""
        formalspin = self.data.homos[0] - self.data.homos[1]
        totalspin = sum(self.analysis.fragspins)
        assert abs(totalspin - formalspin) < 1.0e-3


class GaussianLPATest:
    """Lowdin Population Analysis test"""

    @classmethod
    def setup_class(cls) -> None:
        cls.data = parse()
        cls.analysis = LPA(cls.data)
        cls.analysis.logger.setLevel(0)
        cls.analysis.calculate()

    def testsumcharges(self) -> None:
        """Do the Lowdin charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        assert abs(totalpopulation - formalcharge) < 0.001

    def testsumspins(self) -> None:
        """Do the Lowdin spins sum up to the total formal spin?"""
        formalspin = self.data.homos[0] - self.data.homos[1]
        totalspin = sum(self.analysis.fragspins)
        assert abs(totalspin - formalspin) < 1.0e-3


class GaussianCSPATest:
    """C-squared Population Analysis test"""

    @classmethod
    def setup_class(cls) -> None:
        cls.data = parse()
        cls.analysis = CSPA(cls.data)
        cls.analysis.logger.setLevel(0)
        cls.analysis.calculate()

    def testsumcharges(self) -> None:
        """Do the CSPA charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        assert abs(totalpopulation - formalcharge) < 1.0e-3

    def testsumspins(self) -> None:
        """Do the CSPA spins sum up to the total formal spin?"""
        formalspin = self.data.homos[0] - self.data.homos[1]
        totalspin = sum(self.analysis.fragspins)
        assert abs(totalspin - formalspin) < 1.0e-3


class GaussianBickelhauptTest:
    """Bickelhaupt Population Analysis test"""

    @classmethod
    def setup_class(cls) -> None:
        cls.data = parse()
        cls.analysis = Bickelhaupt(cls.data)
        cls.analysis.logger.setLevel(0)
        cls.analysis.calculate()

    def testsumcharges(self) -> None:
        """Do the Bickelhaupt charges sum up to the total formal charge?"""
        formalcharge = sum(self.data.atomnos) - self.data.charge
        totalpopulation = sum(self.analysis.fragcharges)
        assert abs(totalpopulation - formalcharge) < 1.0e-3

    def testsumspins(self) -> None:
        """Do the Bickelhaupt spins sum up to the total formal spin?"""
        formalspin = self.data.homos[0] - self.data.homos[1]
        totalspin = sum(self.analysis.fragspins)
        assert abs(totalspin - formalspin) < 1.0e-3

    def test_dvb_sp(self) -> None:
        """Testing Bickelhaupt charges (restricted) against outputs from Multiwfn."""
        data, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_sp.out"])
        bpa = Bickelhaupt(data)
        bpa.logger.setLevel(logging.ERROR)
        bpa.calculate()

        e_bpa = numpy.loadtxt(f"{os.path.dirname(os.path.realpath(__file__))}/dvb_sp.bpa")
        assert numpy.all(bpa.fragcharges >= e_bpa - 0.05)
        assert numpy.all(bpa.fragcharges <= e_bpa + 0.05)

    def test_dvb_un_sp(self) -> None:
        """Testing Bickelhaupt charges (unrestricted) against outputs from Multiwfn."""
        data = parse()
        bpa = Bickelhaupt(data)
        bpa.logger.setLevel(logging.ERROR)
        bpa.calculate()

        e_bpaalpha = numpy.loadtxt(f"{os.path.dirname(os.path.realpath(__file__))}/dvb_un_sp.bpa")
        e_bpaspin = numpy.loadtxt(
            f"{os.path.dirname(os.path.realpath(__file__))}/dvb_un_sp.bpaspin"
        )

        assert numpy.all(bpa.fragcharges >= e_bpaalpha - 0.05)
        assert numpy.all(bpa.fragcharges <= e_bpaalpha + 0.05)
        assert numpy.all(bpa.fragspins >= e_bpaspin - 0.05)
        assert numpy.all(bpa.fragspins <= e_bpaspin + 0.05)


def parse():
    data, _ = getdatafile(Gaussian, "basicGaussian09", ["dvb_un_sp.log"])
    return data
