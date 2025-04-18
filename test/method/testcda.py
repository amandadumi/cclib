# Copyright (c) 2025, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Test the CDA method in cclib"""

import logging
import sys

sys.path.insert(1, "..")

from cclib.method import CDA
from cclib.parser import Gaussian

from ..test_data import getdatafile


def main(log: bool = True) -> CDA:
    data1, logfile1 = getdatafile(Gaussian, "CDA", ["BH3CO-sp.log"])
    data2, logfile2 = getdatafile(Gaussian, "CDA", ["BH3.log"])
    data3, logfile3 = getdatafile(Gaussian, "CDA", ["CO.log"])
    fa = CDA(data1)
    if not log:
        fa.logger.setLevel(logging.ERROR)
    fa.calculate([data2, data3])

    return fa


def printResults() -> None:
    fa = main()

    print("       d       b       r")
    print("---------------------------")

    spin = 0
    for i in range(len(fa.donations[0])):
        print(
            f"{int(i):2}: {fa.donations[spin][i]:7.3f} {fa.bdonations[spin][i]:7.3f} {fa.repulsions[spin][i]:7.3f}"
        )

    print("---------------------------")
    print(
        f"T:  {fa.donations[0].sum():7.3f} {fa.bdonations[0].sum():7.3f} {fa.repulsions[0].sum():7.3f}"
    )
    print("\n\n")


class CDATest:
    def runTest(self) -> None:
        """Testing CDA results against Frenking's code"""
        fa = main(log=False)

        donation = fa.donations[0].sum()
        bdonation = fa.bdonations[0].sum()
        repulsion = fa.repulsions[0].sum()

        assert round(abs(donation - 0.181), 3) == 0
        assert round(abs(bdonation - 0.471), 3) == 0
        assert round(abs(repulsion - -0.334), 3) == 0


if __name__ == "__main__":
    printResults()
