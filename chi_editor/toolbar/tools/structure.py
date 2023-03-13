import rdkit.Chem.rdDepictor
from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsPixmapItem
from datamol import incorrect_valence
from rdkit import Chem

from ...bases.alpha_atom import AlphaAtom
from ...bases.tool import Tool
from ...chem_bonds.double_bond import DoubleBond
from ...chem_bonds.single_bond import SingleBond
from ...chem_bonds.triple_bond import TripleBond
from ...playground import mol_from_graphs, matrix_from_item


def create_molecule(atom: AlphaAtom) -> (Chem.Mol, list):
    molecule_matrix = matrix_from_item(atom)
    molecule_smiles = Chem.MolToSmiles(mol_from_graphs(atom.molecule))
    molecule_dm = Chem.MolFromSmiles(molecule_smiles)
    Chem.Kekulize(molecule_dm)
    return molecule_dm, molecule_matrix


def create_atoms(molecule: Chem.Mol, position) -> list[AlphaAtom]:
    Chem.rdDepictor.Compute2DCoords(molecule)
    atoms: list[AlphaAtom] = []

    molecule_center: QPointF = get_geometrical_center(
        [
            QPointF(
                molecule.GetConformer().GetAtomPosition(i).x,
                molecule.GetConformer().GetAtomPosition(i).y,
            )
            for i in range(molecule.GetNumAtoms())
        ]
    )

    for atom in molecule.GetAtoms():
        positions = molecule.GetConformer().GetAtomPosition(atom.GetIdx())
        new_atom = AlphaAtom(atom.GetSymbol())

        new_atom.setPos(
            position + (QPointF(positions.x, positions.y) - molecule_center) * 100,
        )

        atoms.append(new_atom)

    return atoms


def get_geometrical_center(points: list[QPointF]) -> QPointF:
    x: float = sum(point.x() for point in points)
    y: float = sum(point.y() for point in points)
    atoms_count: int = len(points)
    return QPointF(x / atoms_count, y / atoms_count)


class Structure(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            items: list[AlphaAtom] = self.canvas.items(
                event.scenePos(), Qt.ItemSelectionMode.IntersectsItemShape
            )
            if items == [] or not isinstance(items[0], AlphaAtom):
                return super(Structure, self).mouse_press_event(event)
            molecule, molecule_matrix = create_molecule(items[0])

            if not self.check_correctness(molecule, molecule_matrix):
                return

            atoms = create_atoms(molecule, items[0].molecule.molecule_drawer.pos())

            for atom in atoms:
                atom.add_to_canvas(self.canvas)
            self.put_bonds(molecule, atoms)
            items[0].molecule.destroy()
        else:
            atoms = []
            old_atoms = []
            for item in filter(lambda x: isinstance(x, AlphaAtom), self.canvas.items()):
                if item not in old_atoms:
                    molecule, molecule_matrix = create_molecule(item)
                    if not self.check_correctness(molecule, molecule_matrix):
                        return
                    old_atoms.extend(item.molecule.atoms)
                    new_atoms = create_atoms(
                        molecule, item.molecule.molecule_drawer.pos()
                    )
                    atoms.extend(new_atoms)
                    self.put_bonds(molecule, new_atoms)
                    item.molecule.destroy()
            for atom in atoms:
                atom.add_to_canvas(self.canvas)

    def check_correctness(self, molecule: Chem.Mol, molecule_matrix: list) -> bool:
        if molecule is None or incorrect_valence(molecule):
            image = QImage("resources//stathem.jpg")
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
            new_bond = {
                bond_type == 1: SingleBond(atoms[start_position], atoms[end_position]),
                bond_type == 2: DoubleBond(atoms[start_position], atoms[end_position]),
                bond_type == 3: TripleBond(atoms[start_position], atoms[end_position]),
            }[True]

            atoms[start_position].add_line(new_bond)
            atoms[end_position].add_line(new_bond)

            self.canvas.addItem(new_bond)

    @property
    def asset(self) -> str:
        return "structure"
