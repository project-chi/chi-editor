import rdkit.Chem.rdDepictor
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem, QGraphicsPixmapItem
from PyQt6.QtCore import Qt, QPointF
from rdkit import Chem
from datamol import incorrect_valence

from ...bases.alpha_atom import AlphaAtom
from ...bases.tool import Tool
from ...chem_bonds.double_bond import DoubleBond
from ...chem_bonds.single_bond import SingleBond
from ...chem_bonds.triple_bond import TripleBond
from ...playground import mol_from_graphs, matrix_from_item


def create_molecule(atom: AlphaAtom) -> (Chem.Mol, list):
    molecule_matrix = matrix_from_item(atom)
    nodes = molecule_matrix[0]
    adjacent = molecule_matrix[1]
    molecule_smiles = Chem.MolToSmiles(mol_from_graphs(nodes, adjacent))
    molecule_dm = Chem.MolFromSmiles(molecule_smiles)
    Chem.Kekulize(molecule_dm)
    return molecule_dm, molecule_matrix


def create_atoms(molecule: Chem.Mol, position) -> list:
    Chem.rdDepictor.Compute2DCoords(molecule)
    atoms = [None for _ in range(molecule.GetNumAtoms())]

    molecule_center: QPointF = get_geometrical_center(
        [
            QPointF(
                molecule.GetConformer().GetAtomPosition(i).x,
                molecule.GetConformer().GetAtomPosition(i).y
            )
            for i in range(molecule.GetNumAtoms())
        ]
    )

    for i, atom in enumerate(molecule.GetAtoms()):
        positions = molecule.GetConformer().GetAtomPosition(i)
        new_atom = AlphaAtom(atom.GetSymbol())

        new_atom.setPos(position.x() + (positions.x - molecule_center.x()) * 100,
                        position.y() + (positions.y - molecule_center.y()) * 100)

        new_atom.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        new_atom.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
        atoms[i] = new_atom
    return atoms


def create_bonds(molecule: Chem.Mol, atoms: list) -> list:
    bonds = []
    for bond in molecule.GetBonds():
        start_position = bond.GetBeginAtomIdx()
        end_position = bond.GetEndAtomIdx()
        bond_type = bond.GetBondTypeAsDouble()
        if bond_type == 1:
            new_bond = SingleBond(atoms[start_position], atoms[end_position])
        elif bond_type == 2:
            new_bond = DoubleBond(atoms[start_position], atoms[end_position])
        elif bond_type == 3:
            new_bond = TripleBond(atoms[start_position], atoms[end_position])
        atoms[start_position].add_line(new_bond)
        atoms[end_position].add_line(new_bond)
        bonds.append(bond)
    return bonds


def get_geometrical_center(points: list[QPointF]) -> QPointF:
    x: float = sum(point.x() for point in points)
    y: float = sum(point.y() for point in points)
    atoms_count: int = len(points)
    return QPointF(x/atoms_count, y/atoms_count)


class Structure(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            items: list[AlphaAtom] = self.canvas.items(event.scenePos(), Qt.ItemSelectionMode.IntersectsItemShape)
            if items == [] or not isinstance(items[0], AlphaAtom):
                return super(Structure, self).mouse_press_event(event)
            molecule, molecule_matrix = create_molecule(items[0])

            if not self.check_correctness(molecule, molecule_matrix):
                return

            print(items[0].get_molecule_atoms())

            atoms = create_atoms(
                molecule,
                get_geometrical_center(
                    [atom.pos() for atom in items[0].get_molecule_atoms()]
                )
            )

            for atom in atoms:
                self.canvas.addItem(atom)
            self.put_bonds(molecule, atoms)
            self.remove_obsolete(molecule_matrix)
        else:
            atoms = []
            old_atoms = []
            for item in filter(lambda x: isinstance(x, AlphaAtom), self.canvas.items()):
                if item not in old_atoms:
                    molecule, molecule_matrix = create_molecule(item)
                    if not self.check_correctness(molecule, molecule_matrix):
                        return
                    cur_atoms = molecule_matrix[2]
                    old_atoms.extend(cur_atoms)
                    new_atoms = atoms = create_atoms(
                        molecule,
                        get_geometrical_center(
                            [atom.pos() for atom in item.get_molecule_atoms()]
                        )
                    )
                    atoms.extend(new_atoms)
                    self.put_bonds(molecule, new_atoms)
                    self.remove_obsolete(molecule_matrix)
            for atom in atoms:
                self.canvas.addItem(atom)

    def check_correctness(self, molecule: Chem.Mol, molecule_matrix: list) -> bool:
        if molecule is None or incorrect_valence(molecule):
            image = QImage('resources//stathem.jpg')
            molecule = QGraphicsPixmapItem(QPixmap.fromImage(image))
            for i in molecule_matrix[2]:
                for j in i.lines:
                    self.canvas.removeItem(j)
                self.canvas.removeItem(i)
            self.canvas.addItem(molecule)
            return False
        return True

    def put_bonds(self, molecule: Chem.Mol, atoms: list):
        for bond in molecule.GetBonds():
            start_position = bond.GetBeginAtomIdx()
            end_position = bond.GetEndAtomIdx()
            bond_type = bond.GetBondTypeAsDouble()
            if bond_type == 1:
                new_bond = SingleBond(atoms[start_position], atoms[end_position])
            elif bond_type == 2:
                new_bond = DoubleBond(atoms[start_position], atoms[end_position])
            elif bond_type == 3:
                new_bond = TripleBond(atoms[start_position], atoms[end_position])

            atoms[start_position].add_line(new_bond)
            atoms[end_position].add_line(new_bond)
            self.canvas.addItem(new_bond)

    def remove_obsolete(self, molecule_matrix):
        for i in molecule_matrix[2]:
            for j in i.lines:
                self.canvas.removeItem(j)
            self.canvas.removeItem(i)
        self.canvas.selectedItems()

    @property
    def asset(self) -> str:
        return 'structure'
