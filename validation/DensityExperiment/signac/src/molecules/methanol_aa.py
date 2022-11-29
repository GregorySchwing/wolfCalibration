"""Create atomistic representation of methanol."""
import os

import mbuild as mb

from reproducibility_project.src import molecules


class MethanolAA(mb.Compound):
    """Create a single atomistic methanol compound."""

    def __init__(self):
        super(MethanolAA, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/methanol_aa.mol2"), label="MTO")


def main():
    """Create a MethanolAA compound and print basic properties."""
    methanol = MethanolAA()
    print(methanol)
    print(methanol.name)
    print(methanol.labels)
    print(methanol["MTO"])


if __name__ == "__main__":
    main()
