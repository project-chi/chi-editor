from typing import TYPE_CHECKING
from weakref import WeakSet

from chi_editor.bases.molecule.molecule_anchor import MoleculeAnchor

if TYPE_CHECKING:
    from chi_editor.bases.alpha_atom import AlphaAtom


class Molecule:
    atoms: "WeakSet[AlphaAtom]"
    anchor: "MoleculeAnchor"

    def __init__(self, atom: "AlphaAtom") -> None:
        self.atoms = WeakSet()
        self.atoms.add(atom)
        self.anchor = MoleculeAnchor(self.atoms)

    def add_atom(self, atom: "AlphaAtom") -> None:
        self.atoms.add(atom)
        atom.molecule = self

    def remove_atom(self, atom: "AlphaAtom") -> None:
        self.atoms.remove(atom)
        if len(self.atoms) == 0:
            self.anchor.remove()

    def destroy(self):
        atoms_to_remove: "WeakSet[AlphaAtom]" = self.atoms.copy()
        for atom in atoms_to_remove:
            atom.remove()

    def update_atoms(self) -> None:
        queue: "list[AlphaAtom]" = [atom for atom in self.atoms]
        while queue:
            current_atom: "AlphaAtom" = queue.pop()
            if current_atom.molecule != self:
                current_atom.molecule.remove_atom(current_atom)
            self.add_atom(current_atom)
            queue.extend(
                [
                    x
                    for x in current_atom.get_adjacent_atoms()
                    if x not in self.atoms and x not in queue
                ]
            )
        self.anchor.update_position([atom for atom in self.atoms])
        return [atom for atom in self.atoms]
