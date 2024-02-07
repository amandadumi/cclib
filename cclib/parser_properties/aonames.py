from cclib.parser_properties import utils
from cclib.parser_properties.base_parser import base_parser

import numpy as np


class aonames(base_parser):
    """
    Docstring? Units?
    """

    known_codes = ["gaussian"]
    @staticmethod
    def gaussian(file_handler, ccdata) -> list | None:
        dependency_list = ["nmo", "nbasis"]
        line = file_handler.last_line
        if (
            line[5:35] == "Molecular Orbital Coefficients"
            or line[5:41] == "Alpha Molecular Orbital Coefficients"
        ):
            if not base_parser.check_dependencies(dependency_list, ccdata, "atombasis"):
                return None
            aonames = []
            colmNames = file_handler.virtual_next()
            for base in range(0, ccdata.nmo, 5):
                symmetries = file_handler.virtual_next()
                eigenvalues = file_handler.virtual_next()
                for i in range(ccdata.nbasis):
                    line = file_handler.virtual_next()
                    if i == 0:
                        # Find location of the start of the basis function name
                        start_of_basis_fn_name = line.find(line.split()[3]) - 1
                    if base == 0:  # Just do this the first time 'round
                        parts = line[:start_of_basis_fn_name].split()
                        if len(parts) > 1:  # New atom
                            atomname = f"{parts[2]}{parts[1]}"
                        orbital = line[start_of_basis_fn_name:20].strip()
                        aonames.append(f"{atomname}_{orbital}")
            return aonames
        return None

    @staticmethod
    def parse(file_handler, program: str, ccdata) -> list | None:
        constructed_data = None
        if program in aonames.known_codes:
            file_handler.virtual_set()
            program_parser = getattr(aonames, program)
            constructed_data = program_parser(file_handler, ccdata)
            file_handler.virtual_reset()
        return constructed_data
