"""Create atomistic representation of ethanol."""
import os

import mbuild as mb

from reproducibility_project.src import molecules


class mspce(mb.Compound):
    """Create a single atomistic ethanol compound."""

    def __init__(self):
        super(mspce, self).__init__()
        abs_path = os.path.dirname(os.path.abspath(molecules.__file__))
        self.add(mb.load(f"{abs_path}/mspce.mol2"), label="mspce")


def main():
    """Create a MSPCE compound and print basic properties."""
    mspce = mspce()
    print(mspce)
    print(mspce.name)
    print(mspce.labels)
    print(mspce["mspce"])


if __name__ == "__main__":
    main()
