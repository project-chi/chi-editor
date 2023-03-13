import weakref
from typing import Dict

from rdkit import Chem

from chi_editor.bases.alpha_atom import AlphaAtom
from chi_editor.bases.molecule.molecule import Molecule


def mol_from_graphs(molecule: Molecule):
    # create empty editable mol object
    mol = Chem.RWMol()

    atom_to_index: Dict[AlphaAtom, int] = {}
    # add atoms to mol and keep track of index
    for atom in molecule.atoms:
        atom_to_index[atom] = mol.AddAtom(Chem.Atom(atom.text))

    # add bonds between adjacent atoms
    for atom in atom_to_index:
        for line in atom.lines:
            atom1: int = atom_to_index[line.vertex1]
            atom2: int = atom_to_index[line.vertex2]
            if mol.GetBondBetweenAtoms(atom1, atom2) is None:
                bond_type: Chem.rdchem.BondType = {
                    line.multiplicity == 1: Chem.BondType.SINGLE,
                    line.multiplicity == 2: Chem.BondType.DOUBLE,
                    line.multiplicity == 3: Chem.BondType.TRIPLE,
                }[True]
                mol.AddBond(atom1, atom2, bond_type)

    # Convert RWMol to Mol object
    mol = mol.GetMol()

    return mol


def is_line_between(atom1: AlphaAtom, atom2: AlphaAtom):
    for i in atom1.lines:
        if i.vertex1 == atom2 or i.vertex2 == atom2:
            return i.multiplicity
    return 0


def matrix_from_item(atom: AlphaAtom):
    alpha_atoms: list[AlphaAtom] = [atom for atom in atom.molecule.atoms]
    adjacency: list[list[int]]

    adjacency = list(
        list(0 for _ in range(len(alpha_atoms))) for _ in range(len(alpha_atoms))
    )
    for i in range(len(alpha_atoms)):
        for j in range(len(alpha_atoms)):
            adjacency[i][j] = is_line_between(alpha_atoms[i], alpha_atoms[j])

    return [atom.text for atom in alpha_atoms], adjacency, alpha_atoms
