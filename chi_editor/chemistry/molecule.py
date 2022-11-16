from atom import Atom
from bond import Bond


class Molecule:
    atoms: list[Atom]
    bonds: list[Bond]

    def __init__(self, atoms: list[Atom], bonds: list[Bond]) -> None:
        self.atoms = atoms
        self.bonds = bonds
