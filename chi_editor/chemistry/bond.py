from .atom import Atom


class Bond:
    atoms: list[Atom]
    order: int

    def __init__(self, atoms: list[Atom], order: int) -> None:
        self.atoms = atoms
        self.order = order
