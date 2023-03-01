from rdkit import Chem

from chi_editor.bases.alpha_atom import AlphaAtom


def mol_from_graphs(node_list, adjacency_matrix):

    # create empty editable mol object
    mol = Chem.RWMol()

    # add atoms to mol and keep track of index
    node_to_idx = {}
    for i in range(len(node_list)):
        a = Chem.Atom(node_list[i])
        molIdx = mol.AddAtom(a)
        node_to_idx[i] = molIdx

    # add bonds between adjacent atoms
    for ix, row in enumerate(adjacency_matrix):
        for iy, bond in enumerate(row):

            # only traverse half the matrix
            if iy <= ix:
                continue

            # add relevant bond type (there are many more of these)
            if bond == 0:
                continue
            elif bond == 1:
                bond_type = Chem.rdchem.BondType.SINGLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
            elif bond == 2:
                bond_type = Chem.rdchem.BondType.DOUBLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)
            elif bond == 3:
                bond_type = Chem.rdchem.BondType.TRIPLE
                mol.AddBond(node_to_idx[ix], node_to_idx[iy], bond_type)

    # Convert RWMol to Mol object
    mol = mol.GetMol()

    return mol


def is_line_between(atom1: AlphaAtom, atom2: AlphaAtom):
    for i in atom1.lines:
        if i.vertex1 == atom2 or i.vertex2 == atom2:
            return i.multiplicity
    return 0


def matrix_from_item(atom: AlphaAtom):
    alpha_atoms: list[AlphaAtom]
    atoms: list[str]
    adjacency: list[list[int]]

    alpha_atoms = []
    queue = [atom]
    while queue:
        current_atom = queue.pop(0)
        if current_atom not in alpha_atoms:
            alpha_atoms.append(current_atom)
            queue += (list(map(lambda x: x.vertex2 if x.vertex1 == current_atom else x.vertex1, current_atom.lines)))
    adjacency = list(list(0 for i in range(len(alpha_atoms))) for i in range(len(alpha_atoms)))
    for i in range(len(alpha_atoms)):
        for j in range(len(alpha_atoms)):
            adjacency[i][j] = is_line_between(alpha_atoms[i], alpha_atoms[j])

    return (list(map(lambda x: x.text, alpha_atoms)), adjacency, alpha_atoms)
