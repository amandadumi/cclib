# -*- coding: utf-8 -*-
#
# Copyright (c) 2020, the cclib development team
#
# This file is part of cclib (http://cclib.github.io) and is distributed under
# the terms of the BSD 3-Clause License.

"""Parser for Formatted Checkpoint files"""


from __future__ import print_function

import re

import numpy

from cclib.parser import data
from cclib.parser import logfileparser
from cclib.parser import utils


class CFOUR(logfileparser.Logfile):
    """ A CFOUR output file
    """

    def __init__(self, *args, **kwargs):
        # Call the __init__ method of the superclass
        super(CFOUR, self).__init__(logname="CFOUR", *args, **kwargs)
        self.start = True

    def __str__(self):
        """Return a string representation of the object."""
        return "CFOUR output file {}".format(self.filename)

    def __repr__(self):
        """Return a representation of the object."""
        return 'CFOUR("{}")'.format(self.filename)

    def normalisesym(self, symlabel):
        """Just return label"""
        return symlabel


    def extract(self, inputfile, line):
        # just opened file, skip first line to get basis
        if 'Version' in line:
            package_version = line.split()[1]
            self.metadata["package_version"] = package_version

        if '@GETXYZ-I' in line:
            self.natom = int(line.split()[1])

        if 'Coordinates (in bohr)' in line:
            # skip labeling lines
            #         Coordinates used in calculation (QCOMP)
            # ----------------------------------------------------------------
            # Z-matrix   Atomic            Coordinates (in bohr)
            #  Symbol    Number           X              Y              Z
            # ----------------------------------------------------------------
            #     O         8         0.00000000     0.00000000     0.00000000
            for i in range(2):
                 line = next(inputfile)
            coord_block = numpy.array(self._parse_block(inputfile, self.natom,str, 'parsing coordinate blocks'))
            self.atomnos = numpy.array(coord_block[:,1],dtype=int)
            self.atomcoords = numpy.array(coord_block[:,2:],dtype=float)

        # find the number of basis functions
        ang_mom_map = {'S':'S', 'X':'P', 'XX':'D', 'XXX':'F'}
        if  'There are' in line and 'basis functions.' in line:
            self.nbasis = line.split()[2]
        if  'GAUSSIAN BASIS INFORMATION' in line:
            for i in range(4):
                 line = next(inputfile)
            self.gbasis = [[]]
            split_line = line.split() # example line:  O #1  1    S
            atom_num = int(split_line[1].strip('#'))
            new_atom_num = 0
            ang_mom = split_line[3]
            basis_done = False
            line = next(inputfile)
            count = 0
            new = True
            while not basis_done:
                # data structure to store right shape of coeffs, exp pairs.
                overall_basis = []
                while ('#' not in line) and (basis_done is False):
                    ## the basis section end is marked by a line with two numbers.
                    if len(line.split()) == 2:
                        basis_done=True
                        break
                    # there are new lines between each basis shell.
                    if line =='\n':
                        line = next(inputfile)
                        new = True
                        continue
                    if ang_mom not in ['S','X','XX','XXX']: # CFOUR outputs each component, though cclib only stores as S, P, D.
                        line = next(inputfile)
                        while ('#' not in line):
                            line = next(inputfile)
                            if len(line) == 2:
                                basis_done = True
                        print(line)
                    else:
                        line = line.strip('+')
                        split_line = line.split()
                        exp = float(split_line[1])
                        coeffs = numpy.array(split_line[2:],dtype=float)
                        print(numpy.shape(overall_basis))
                        print(len(overall_basis))
                        for idx, i in enumerate(coeffs):
                            # if you are starting on a new atom, or starting a new angular momenta shell.
                            print('new is {}'.format(new))
                            if new == True or len(overall_basis) ==0:
                                print('we are here')
                                overall_basis.append([(exp,i)])
                                if idx+1 == len(coeffs):
                                    new = False
                            # you are in the same shell
                            else:
                                print('now we are here!')
                                overall_basis[idx].append((exp,i))
                        # pairs = zip(numpy.ones(len(coeffs))*exp,coeffs)
                        line = next(inputfile)
                        # overall_basis.concatenate(basis_list,axis=1)
                for i in range(len(overall_basis)):
                    self.gbasis[-1].append((ang_mom_map[ang_mom],overall_basis[i]))
                if basis_done:
                    break
                print('line after appending new basis is:')
                print(line)
                split_line = line.split() # example line:  O #1  1    S
                new_atom_num = int(split_line[1].strip('#'))
                ang_mom = split_line[3]
                line = next(inputfile)
                if new_atom_num != atom_num:
                    self.gbasis.append([])
                    atom_num = new_atom_num
                count +=1
                if count == 3:
                    exit
                print(self.gbasis)
        else:
            pass
            # raise Warning('CFOUR job was not run with PRINT=1 option, basis is not detailed enough')




    def _parse_block(self, inputfile, count, type, msg):
        """
        pases through a block
        PARAMETERS
        ---
        inputfile: str
          input file currently being parsed
        count: int
          the number of lines to parse
        type: object type definition
          type of data returned
        msg: str
          message to pass to the update log.
        """
        property = []
        while len(property) < count :
            self.updateprogress(inputfile, msg, self.fupdate)
            line = next(inputfile)
            property.append([ type(x) for x in line.split()])
            print(property)
        print(numpy.shape(property))
        return property