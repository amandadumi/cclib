from dataclasses import dataclass

import cclib.parser_properties as cprops


@dataclass
class combinator:
    name: str
    job_list: list


@dataclass
class sp_combinator(combinator):
    def __init__(self):
        self.name = "single_point"
        self.job_list = [[cprops.scfenergies,cprops.nmo,cprops.atommasses,cprops.charge,cprops.nbasis,cprops.atommasses,cprops.atomcoords]]