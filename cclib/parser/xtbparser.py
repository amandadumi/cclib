import re
from itertools import groupby
from typing import List, Optional

from cclib.parser import logfileparser

import numpy as np


class XTB(logfileparser.Logfile):
    """An output parser for the xTB code"""

    def __init__(self, *args, **kwargs):
        super().__init__(logname="xTB", *args, **kwargs)

    def __str__(self) -> str:
        """Return a string representation of the object."""
        return f"xTB log file {self.filename}"

    def __repr__(self) -> str:
        """Return a representation of the object."""
        return f'xTB("{self.filename}")'

    def normalisesym(self, label):
        """xTB does not require normalizing symmetry labels."""
        return label

    def before_parsing(self) -> None:
        """Set attributes before parsing"""
        self.atomprop = {}
        self.bondprop = {}

    def after_parsing(self) -> None:
        """Delete empty attributes after parsing"""
        if not self.atomprop:
            delattr(self, "atomprop")
        if not self.bondprop:
            delattr(self, "bondprop")

    def _extract_version(self, line: str) -> Optional[str]:
        """
        Extract xtb version from the following:

              -----------------------------------------------------------
             |                   =====================                   |
             |                           x T B                           |
             |                   =====================                   |
             |                         S. Grimme                         |
             |          Mulliken Center for Theoretical Chemistry        |
             |                    University of Bonn                     |
              -----------------------------------------------------------

           * xtb version 6.6.1 (8d0f1dd) compiled by 'conda@1efc2f54142f' on 2023-08-01
        """

        version_match = re.search(r"xtb version (\d+(\.\d+)+)", line)
        return version_match.group(1) if version_match else None

    def _extract_coord_file(self, line: str) -> Optional[str]:
        """
        Extract the coordinate filename, from which we can strip out the type

           -------------------------------------------------
          |                Calculation Setup                |
           -------------------------------------------------

          program call               : xtb coord.xyz --opt
          coordinate file            : coord.xyz
          omp threads                :                    20
        """
        coordinate_file_match = re.search(r"coordinate file\s+:\s+(\S+)", line)
        return coordinate_file_match.group(1) if coordinate_file_match else None

    def _extract_charge(self, line: str) -> Optional[int]:
        """
        Extract the total charge. It can always be found in the
        summary formatted as a float.

         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         ::                     SUMMARY                     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
         :: total energy              -5.070544440612 Eh    ::
         :: gradient norm              0.000057326562 Eh/a0 ::
         :: HOMO-LUMO gap             14.391809984508 eV    ::
         ::.................................................::
         :: SCC energy                -5.104920280363 Eh    ::
         :: -> isotropic ES            0.031458595179 Eh    ::
         :: -> anisotropic ES          0.000396760551 Eh    ::
         :: -> anisotropic XC         -0.000881430881 Eh    ::
         :: -> dispersion             -0.000141085082 Eh    ::
         :: repulsion energy           0.034375839725 Eh    ::
         :: add. restraining           0.000000000000 Eh    ::
         :: total charge              -0.000000000000 e     ::
         :::::::::::::::::::::::::::::::::::::::::::::::::::::
        """
        charge_match = re.search(r"total charge\s+([-+]?\d+\.\d+)", line)
        return round(float(charge_match.group(1))) if charge_match else None

    def _extract_final_energy(self, line: str) -> Optional[float]:
        """
        Extract the final total energy from the result table.

           -------------------------------------------------
          | TOTAL ENERGY               -5.070544323569 Eh   |
          | GRADIENT NORM               0.000458081396 Eh/α |
          | HOMO-LUMO GAP              14.381252816459 eV   |
           -------------------------------------------------
        """
        total_energy_match = re.search(r"TOTAL ENERGY\s+([-+]?\d+\.\d+)", line)
        return float(total_energy_match.group(1)) if total_energy_match else None

    def _extract_geom_energy(self, line: str) -> Optional[float]:
        """
                Extract the energies for a geometry step.

        ........................................................................
        .............................. CYCLE    1 ..............................
        ........................................................................

                 iter      E             dE          RMSdq      gap      omega  full diag
                   1     -5.1048382 -0.510484E+01  0.417E-06   14.38       0.0  T
                   2     -5.1048382  0.000000E+00  0.234E-06   14.38   24706.5  T
                   3     -5.1048382  0.000000E+00  0.437E-07   14.38  100000.0  T
                     SCC iter.                  ...        0 min,  0.005 sec
                     gradient                   ...        0 min,  0.010 sec
                 * total energy  :    -5.0705443 Eh     change       -0.4369838E-12 Eh
                   gradient norm :     0.0004582 Eh/α   predicted     0.0000000E+00 (-100.00%)
                   displ. norm   :     0.0005728 α      lambda       -0.1688374E-06
                   maximum displ.:     0.0005029 α      in ANC's #2, #3, #1, ...

                   *** GEOMETRY OPTIMIZATION CONVERGED AFTER 1 ITERATIONS ***
        """

        geom_energy_match = re.search(r"\*\s+total energy\s+:\s+([-+]?\d+\.\d+)", line)
        return float(geom_energy_match.group(1)) if geom_energy_match else None

    def _extract_final_structure(self, line: str) -> Optional[str]:
        """
        Extract the final structure for a geometry optimization.

        For XYZ files:

        ================
        final structure:
        ================
        3
        xtb: 6.6.1 (8d0f1dd)
        O            1.07015391331798       -0.01769828395654        0.04981203402603
        H            2.02952514441169       -0.00813780275851        0.03338237327689
        H            0.76845094227033        0.44031608671506       -0.73728440730292

        For SDF/mol files:

        """

    def _is_cycle_line(self, line: str) -> bool:
        """
        Extract if the line indicates it is a geometry optimization.

        ........................................................................
        .............................. CYCLE    1 ..............................
        ........................................................................
        """

        return bool(re.search(r"CYCLE\s+\d+", line))

    def _is_geom_end_line(self, line: str) -> bool:
        """
        Extract if the line indicates the optimization is over.

        *** GEOMETRY OPTIMIZATION CONVERGED AFTER 1 ITERATIONS ***

        or

        *** FAILED TO CONVERGE GEOMETRY OPTIMIZATION ***
        """
        return bool(re.search(r"\s+GEOMETRY OPTIMIZATION\s+", line))

    def _is_finished(self, line: str) -> bool:
        """
        Extract if the job finished.

        ------------------------------------------------------------------------
        * finished run on 2023/10/26 at 13:18:28.705
        ------------------------------------------------------------------------
        """
        return bool(re.search(r"\*\s+finished run on", line))

    def extract(self, inputfile: List[str], line: str) -> None:
        # Get the xTB version
        version = self._extract_version(line)
        if version:
            self.metadata["legacy_package_version"] = version

        # Get the coordinate file type
        coord_filename = self._extract_coord_file(line)
        if coord_filename:
            self.metadata["coord_type"] = coord_filename.split(".")[-1]

        # TODO: How to handle mult...
        # self.set_attribute("mult", mult)

        # Get the net charge
        charge = self._extract_charge(line)
        if charge:
            self.set_attribute("charge", charge)

        # Get if it's a marker indicating that it's an optimization
        is_geom_opt = self._is_cycle_line(line)

        # Cycle through the gemoetry steps to get the total energies,
        # if applicable
        scf_energies = []
        if is_geom_opt:
            while not self._is_geom_end_line(line):
                scf_energy = self._extract_geom_energy(line)
                if scf_energy:
                    scf_energies.append(scf_energy)
                line = next(inputfile)

            # Get the final geometry
            if line.strip()[:15] == "final structure":
                self.skip_line(inputfile, "=")

                if self.metadata["coord_type"] == "xyz":
                    atomnos = []
                    atomcoords = []
                    for line in inputfile:
                        # Ending criteria for xyz is a blank line at the end of the coords block
                        if line == " \n":
                            break
                        if line[0].isupper():
                            atom, x, y, z = line.split()
                            atomnos.append(self.table.number[atom])
                            atomcoords.append([float(x), float(y), float(z)])
                    self.set_attribute("natom", len(atomnos))
                    self.set_attribute("atomnos", atomnos)
                    self.set_attribute("atomcoords", atomcoords)

                elif self.metadata["coord_type"] in ("sdf", "mol"):
                    atomnos = []
                    atomcoords = []
                    # Ending criteria for sdf\mol is the END at the end of the coord block
                    while line.strip()[-3:] != "END":
                        # Atoms block start with 3 blank spaces, bonds block starts with 1
                        if line[:3] == "   ":
                            x, y, z, atom = line.split()[:4]
                            atomnos.append(self.table.number[atom])
                            atomcoords.append([float(x), float(y), float(z)])
                        line = next(inputfile)
                    self.set_attribute("natom", len(atomnos))
                    self.set_attribute("atomnos", atomnos)
                    self.set_attribute("atomcoords", atomcoords)

        # Get Molecular Orbitals energies and HOMO index
        # xTB trunctaes the MO list so we need to take care of that.
        # Unkown energies will be given NaN as a value
        #
        # * Orbital Energies and Occupations

        #      #    Occupation            Energy/Eh            Energy/eV
        #   -------------------------------------------------------------
        #      1        2.0000           -0.7817342             -21.2721
        #    ...           ...                  ...                  ...
        #     21        2.0000           -0.5177364             -14.0883
        #     22        2.0000           -0.5133906             -13.9701
        #     23        2.0000           -0.5119411             -13.9306
        #     24        2.0000           -0.5103339             -13.8869
        #     25        2.0000           -0.5064217             -13.7804
        #     26        2.0000           -0.4793904             -13.0449
        #     27        2.0000           -0.4762317             -12.9589
        #     28        2.0000           -0.4705819             -12.8052
        #     29        2.0000           -0.4558376             -12.4040
        #     30        2.0000           -0.4505134             -12.2591
        #     31        2.0000           -0.4390552             -11.9473
        #     32        2.0000           -0.4371482             -11.8954
        #     33        2.0000           -0.4083272             -11.1111 (HOMO)
        #     34                         -0.2990289              -8.1370 (LUMO)
        #     35                         -0.2703399              -7.3563
        #     36                         -0.2376187              -6.4659
        #     37                         -0.2246900              -6.1141
        #     38                         -0.2213822              -6.0241
        #     39                         -0.2016539              -5.4873
        #     40                         -0.1317437              -3.5849
        #     41                         -0.1173862              -3.1942
        #     42                          0.0207011               0.5633
        #     43                          0.0378419               1.0297
        #     44                          0.0843351               2.2949
        #    ...                                ...                  ...
        #     60                          1.1799189              32.1072
        #   -------------------------------------------------------------
        #               HL-Gap            0.1092983 Eh            2.9742 eV
        #          Fermi-level           -0.3536781 Eh           -9.6241 eV
        #

        if line.strip() == "* Orbital Energies and Occupations":
            # Skip 4 lines to get to the table
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)

            mooccnos = [[]]
            moenergies = [[]]
            monumbers = []

            # Ending criteria is the dashed line at the end of the MO block.
            while line.strip() != "-------------------------------------------------------------":
                line_split = line.split()
                monumbers.append(line_split[0])

                # Parsing the lines before the LUMO line and the HOMO line itself.
                # All MOs are occupied
                if (len(line_split) == 4 and line_split[-1] != "(LUMO)") or (len(line_split) == 5):
                    mooccnos[0].append(line_split[1])
                    moenergies[0].append(line_split[3])

                # For the LUMO line and after we assume 0 electrons
                # Since they are not explicit
                if len(line_split) == 3 or line_split[-1] == "(LUMO)":
                    mooccnos[0].append(0.0)
                    moenergies[0].append(line_split[2])

                # xTB gives the index of the HOMO
                # Occupation can be not an integer which complicates the way
                # to parse the HOMOs for unrestricted and openshell calculations.
                # Keeping it this way for now.
                if line_split[-1] == "(HOMO)" and not hasattr(self, "homos"):
                    self.set_attribute("homos", [int(line_split[0]) - 1])

                line = next(inputfile)

            # Find in index of the "..."
            # And fixing the type of the rest of the values
            fill_in_idx = []
            for idx, monumber in enumerate(monumbers):
                if monumber == "...":
                    fill_in_idx.append(idx)
                else:
                    monumbers[idx] = int(monumbers[idx])
                    moenergies[0][idx] = float(moenergies[0][idx])
                    mooccnos[0][idx] = float(mooccnos[0][idx])

            # Filling in the gaps since xTB truncates the list to include only
            # the first and last MOs, and some number of MOs before and
            # after the HOMO and LUMO.
            #
            # Electron occupency is assumed to be 2 for the missing values
            # before the HOMO, and 0 for the missing values after the LUMO.
            #
            # NaN is put for the missing MO energies since we have no way
            # to extrapolate them.

            missing_mos_num = 0
            for idx in fill_in_idx:
                fixed_idx = idx + missing_mos_num
                first_mo = int(monumbers[fixed_idx - 1])
                last_mo = int(monumbers[fixed_idx + 1])
                missing_mos_num = last_mo - first_mo - 1
                na_list = [np.nan] * missing_mos_num

                if mooccnos[0][fixed_idx - 1] == 2:
                    twos_list = [2.0] * missing_mos_num
                    mooccnos[0] = mooccnos[0][:fixed_idx] + twos_list + mooccnos[0][fixed_idx + 1 :]
                else:
                    zeros_list = [0.0] * missing_mos_num
                    mooccnos[0] = (
                        mooccnos[0][:fixed_idx] + zeros_list + mooccnos[0][fixed_idx + 1 :]
                    )
                monumbers = (
                    monumbers[:fixed_idx]
                    + list(np.arange(first_mo + 1, last_mo))
                    + monumbers[fixed_idx + 1 :]
                )
                moenergies[0] = moenergies[0][:fixed_idx] + na_list + moenergies[0][fixed_idx + 1 :]

                missing_mos_num -= 1

            self.set_attribute("moenergies", moenergies)

        # Grabbing atomic properties: Coordination number CN, Atomic partial charge q, Dispersion coefficient C6, Polarizability alpha:
        #
        #  #   Z          covCN         q      C6AA      α(0)
        #  1   6 C        3.056     0.087    26.011     8.364
        #  2   8 O        1.729    -0.155    16.635     5.507
        # ...
        # 17   1 H        0.927     0.059     2.210     2.325
        # 18   1 H        0.926     0.066     2.131     2.283
        #
        if line.strip()[:5] == "#   Z":
            line = next(inputfile)
            atom_convcn = []
            atom_q = []
            atom_c6aa = []
            atom_alpha = []
            while line.strip() != "":
                line_split = line.strip().split()
                atom_convcn.append(float(line_split[3]))
                atom_q.append(float(line_split[4]))
                atom_c6aa.append(float(line_split[5]))
                atom_alpha.append(float(line_split[6]))

                line = next(inputfile)

            self.atomprop["convcn"] = atom_convcn
            self.atomprop["q"] = atom_q
            self.atomprop["c6aa"] = atom_c6aa
            self.atomprop["alpha"] = atom_alpha

        # Grabbing the Wiberg bond orders (WBO)
        #
        # Wiberg/Mayer (AO) data.
        # largest (>0.10) Wiberg bond orders for each atom
        #
        #  ---------------------------------------------------------------------------
        #      #   Z sym  total        # sym  WBO       # sym  WBO       # sym  WBO
        #  ---------------------------------------------------------------------------
        #      1   6 C    3.893 --     3 C    1.482     2 O    1.174    15 H    0.964
        #                              7 C    0.116
        #      2   8 O    2.471 --     1 C    1.174     7 C    1.096
        #      3   6 C    3.982 --     1 C    1.482     6 C    1.254     4 C    1.103
        #      4   6 C    3.951 --     5 N    2.811     3 C    1.103
        #      5   7 N    3.010 --     4 C    2.811
        #      6   6 C    3.982 --     7 C    1.505     3 C    1.254    16 H    0.970
        #      7   6 C    3.913 --     6 C    1.505     8 C    1.116     2 O    1.096
        #                              1 C    0.116
        #      8   6 C    3.919 --    10 C    1.399     9 O    1.131     7 C    1.116
        #                             14 C    0.115
        #      9   8 O    2.468 --    14 C    1.133     8 C    1.131
        #     10   6 C    3.979 --     8 C    1.399    13 C    1.240    11 C    1.119
        #     11   6 C    3.952 --    12 N    2.789    10 C    1.119
        #     12   7 N    2.997 --    11 C    2.789
        #     13   6 C    3.982 --    14 C    1.611    10 C    1.240    17 H    0.973
        #     14   6 C    3.895 --    13 C    1.611     9 O    1.133    18 H    0.964
        #                              8 C    0.115
        #     15   1 H    0.994 --     1 C    0.964
        #     16   1 H    0.996 --     6 C    0.970
        #     17   1 H    0.996 --    13 C    0.973
        #     18   1 H    0.994 --    14 C    0.964
        #  ---------------------------------------------------------------------------
        #

        if line.strip() == "Wiberg/Mayer (AO) data.":
            # Skip 6 lines to get to the first line of data
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)
            line = next(inputfile)

            wbo = []
            # Iterating over lines until the end dashes line
            while (
                line.strip()
                != "---------------------------------------------------------------------------"
            ):
                if line[5] != " ":
                    line_split = line.strip().split()
                    wbo_total = float(line_split[3])

                    wbo.append([wbo_total])
                    for i in range(0, len(line_split) - 5, 3):
                        wbo_next = line_split[i + 5 : i + 8]
                        wbo_idx = int(wbo_next[0]) - 1
                        wbo_order = float(wbo_next[2])

                        wbo[-1].append([wbo_idx, wbo_order])
                else:
                    line_split = line.strip().split()

                    for i in range(0, len(line_split), 3):
                        wbo_next = line_split[i : i + 3]
                        wbo_idx = int(wbo_next[0]) - 1
                        wbo_order = float(wbo_next[2])
                        wbo[-1].append([wbo_idx, wbo_order])
                line = next(inputfile)

            self.bondprop["wbo"] = wbo

        # Get Fukui Indecis for each atom
        #    #        f(+)     f(-)     f(0)
        #  1N      -0.075   -0.084   -0.080
        #  2C      -0.021   -0.014   -0.018
        #  3N      -0.084   -0.092   -0.088
        # ....................
        # 14H      -0.059   -0.057   -0.058
        # 15H      -0.070   -0.065   -0.068
        # 16H      -0.067   -0.064   -0.065
        # Creates a list of lists the size of (n,3) where n is the number of atoms.
        # Each atom had the following info, in that order:
        #   - f(+)
        #   - f(-)
        #   - f(0)

        if line.strip() == "#        f(+)     f(-)     f(0)":
            line = next(inputfile)
            atom_fp = []
            atom_fn = []
            atom_fz = []
            for _ in range(self.natom):
                try:
                    atom_fp.append(float(line[9:19]))
                except Exception:
                    atom_fp.append(-1000000.0)
                try:
                    atom_fn.append(float(line[19:28]))
                except Exception:
                    atom_fn.append(-1000000.0)
                try:
                    atom_fz.append(float(line[28:]))
                except Exception:
                    atom_fz.append(1000000.0)
                line = next(inputfile)

            self.atomprop["fukui"] = [atom_fp, atom_fn, atom_fz]

        # Get LMO data
        #
        #     LMO Fii/eV  ncent    charge center   contributions...
        #     1 sigma -21.33   1.84  12.23319  -5.99934  -6.73018   13O :  0.64   12C :  0.37
        #     2 sigma -21.05   1.82  11.03031  -5.52370  -7.69423   13O :  0.66   31H :  0.34
        #     3 sigma -20.30   1.94  13.15771  -3.57423  -2.57537    7N :  0.58    6C :  0.42
        #     4 sigma -20.29   1.94  13.85375  -4.46664   4.16798   11N :  0.56   10C :  0.45
        # .......................
        #    37 sigma -17.29   1.97  10.58062   1.21250  -3.64843    4C :  0.52   21H :  0.48
        #    38 sigma -17.22   1.98  12.43345   1.40582  -2.52759    4C :  0.52   22H :  0.48
        #    39 sigma -17.22   1.98  13.80996   3.97178  -5.25134    3C :  0.53   19H :  0.48
        #    40 LP    -15.87   1.05  14.66884  -4.88820   5.98676   11N :  0.98
        #    41 LP    -15.58   1.02  13.59986   8.28526  -4.78561    1N :  0.99
        #
        # Creates a list of lists size of (n,4) where n is the number of atoms.
        # Each atom has the following info, in this order:
        #   - Higest priority LMO type encoding (0 for LP, 1 for pi, 2 for delpi, 3 for sigma) - where 0 (LP) is the highest priority.
        #   - Average contribution of the highest priority LMO
        #   - Average Fii/eV of the highest priority LMO
        #   - Average ncent of the highest priority LMO

        if line.startswith(" LMO Fii/eV"):  # and 'donescf' in self.attributes.keys():
            line = next(inputfile)
            lmo_list = []
            while line[:5].strip().isnumeric():
                split = [x for x in re.split("\s+|:", line) if x != ""]
                lmo_num = split[0]
                lmo_type = split[1]
                try:
                    lmo_fii = float(split[2])
                except Exception:
                    lmo_fii = split[2]
                try:
                    lmo_ncent = float(split[3])
                except Exception:
                    lmo_ncent = split[3]

                lmo_cont = split[7:]
                # if lmo_type in ['pi','LP']:
                for i in range(len(lmo_cont) // 2):
                    if (
                        (lmo_type == "pi" and float(lmo_cont[2 * i - 1]) > 0.3)
                        or (lmo_type == "LP" and float(lmo_cont[2 * i - 1]) > 0.7)
                        or (lmo_type in ["sigma", "delpi"])
                    ):
                        count_atoms = re.findall("(\d+|\D+)", lmo_cont[2 * i - 2])
                        lmo_list.append(
                            {
                                "AtomIdx": int(count_atoms[0]),
                                "Contribution": float(lmo_cont[2 * i - 1]),
                                "LMO Num": lmo_num,
                                "LMO Type": lmo_type,
                                "Fii/eV": lmo_fii,
                                "ncent": lmo_ncent,
                            }
                        )
                line = next(inputfile)

            LMO_ORDER = {"LP": 0, "pi": 1, "delpi": 2, "sigma": 3}
            lmo_list = sorted(lmo_list, key=lambda x: (x["AtomIdx"], LMO_ORDER[x["LMO Type"]]))
            lmo_list_cleaned = []
            keys_list = []
            for (key, _), group in groupby(
                lmo_list, key=lambda x: (x["AtomIdx"], LMO_ORDER[x["LMO Type"]])
            ):
                if key not in keys_list:
                    keys_list.append(key)
                    temp_list = list(group)
                    atom_cont = sum(d["Contribution"] for d in temp_list) / len(temp_list)
                    atom_fii = sum(d["Fii/eV"] for d in temp_list) / len(temp_list)
                    atom_ncent = sum(d["ncent"] for d in temp_list) / len(temp_list)
                    lmo_code = LMO_ORDER[temp_list[0]["LMO Type"]]
                    lmo_list_cleaned.append([lmo_code, atom_cont, atom_fii, atom_ncent])

            self.atomprop["lmo"] = lmo_list_cleaned

        # Get the final total energy
        final_energy = self._extract_final_energy(line)

        # Patch the final total energy to be the last SCF energy
        # since it is higher precision and also always available
        if final_energy:
            scf_energies = scf_energies[-1] if scf_energies else [final_energy]

        if scf_energies:
            self.set_attribute("scfenergies", scf_energies)

        # find if job ended successfuly
        is_finished = self._is_finished(line)
        if is_finished:
            self.metadata["success"] = True