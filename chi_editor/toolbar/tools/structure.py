from __future__ import annotations

from datamol import incorrect_valence
from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtWidgets import QGraphicsPixmapItem, QGraphicsSceneMouseEvent
from rdkit import Chem
from rdkit.Chem import Mol

from ...bases.alpha_atom import AlphaAtom
from ...bases.line import Line
from ...bases.tool import Tool
from ...chem_bonds.double_bond import DoubleBond
from ...chem_bonds.single_bond import SingleBond
from ...chem_bonds.triple_bond import TripleBond
from ...chem_utils import mol_from_graphs


def create_atoms(molecule: Chem.Mol, position: QPointF) -> list[AlphaAtom]:
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


def create_bonds(molecule: Chem.Mol, atoms: list[AlphaAtom]) -> list[Line]:
    result: list[Line] = []
    for bond in molecule.GetBonds():
        start_position: int = bond.GetBeginAtomIdx()
        end_position: int = bond.GetEndAtomIdx()
        bond_type: float = bond.GetBondTypeAsDouble()
        new_bond: Line = {
            bond_type == 1: SingleBond(atoms[start_position], atoms[end_position]),
            bond_type == 2: DoubleBond(atoms[start_position], atoms[end_position]),
            bond_type == 3: TripleBond(atoms[start_position], atoms[end_position]),
        }[True]

        atoms[start_position].add_line(new_bond)
        atoms[end_position].add_line(new_bond)
        atoms[start_position].molecule.update_atoms()

        result.append(new_bond)
    return result


def create_molecule(canvas, atom: AlphaAtom) -> Mol | None:
    # molecule_smiles: str = Chem.MolToSmiles(mol_from_graphs(atom.molecule))
    # molecule_dm: Chem.Mol = Chem.MolFromSmiles(molecule_smiles)
    molecule_dm = mol_from_graphs(atom.molecule)
    if not check_correctness(canvas, molecule_dm):
        return None
    Chem.Kekulize(molecule_dm)
    return molecule_dm


def check_correctness(canvas, mol: Chem.Mol) -> bool:
    if mol is None or incorrect_valence(mol):
        image = QImage("resources//stathem.jpg")
        mol = QGraphicsPixmapItem(QPixmap.fromImage(image))
        for item in canvas.items():
            if isinstance(item, AlphaAtom):
                item.molecule.destroy()
        canvas.addItem(mol)
        return False
    return True


def put_bonds(canvas, molecule: Chem.Mol, atoms: list[AlphaAtom]) -> None:
    bonds = create_bonds(molecule, atoms)
    for bond in bonds:
        canvas.addItem(bond)


def put_molecule(canvas, molecule: Chem.Mol, position: QPointF) -> None:
    if molecule is None:
        return
    atoms: list[AlphaAtom] = create_atoms(molecule, position)
    for atom in atoms:
        atom.add_to_canvas(canvas)
    put_bonds(canvas, molecule, atoms)


class Structure(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            items: list = self.canvas.items(
                event.scenePos(), Qt.ItemSelectionMode.IntersectsItemShape
            )
            if items == [] or not isinstance(items[0], AlphaAtom):
                return super(Structure, self).mouse_press_event(event)

            current_atom: AlphaAtom = items[0]
            molecule: Chem.Mol = create_molecule(self.canvas, current_atom)

            if molecule is None:
                return

            atoms: list[AlphaAtom] = create_atoms(
                molecule, current_atom.molecule.anchor.pos()
            )

            for atom in atoms:
                atom.add_to_canvas(self.canvas)
            put_bonds(self.canvas, molecule, atoms)
            current_atom.molecule.destroy()
        else:
            atoms: list[AlphaAtom] = []
            old_atoms: list[AlphaAtom] = []
            for item in self.canvas.items():
                if isinstance(item, AlphaAtom) and item not in old_atoms:
                    molecule = create_molecule(self.canvas, item)
                    if molecule is None:
                        return
                    old_atoms.extend(item.molecule.atoms)
                    new_atoms = create_atoms(molecule, item.molecule.anchor.pos())
                    atoms.extend(new_atoms)
                    put_bonds(self.canvas, molecule, new_atoms)
                    item.molecule.destroy()
            for atom in atoms:
                atom.add_to_canvas(self.canvas)

    @property
    def picture(self) -> str:
        return "structure"
