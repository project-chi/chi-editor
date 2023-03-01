import rdkit.Chem.rdDepictor
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem, QGraphicsPixmapItem
from PyQt6.QtCore import Qt
from rdkit import Chem
from datamol import incorrect_valence

from ...bases.alpha_atom import AlphaAtom
from ...bases.tool import Tool
from ...chem_bonds.double_bond import DoubleBond
from ...chem_bonds.single_bond import SingleBond
from ...chem_bonds.triple_bond import TripleBond
from ...playground import mol_from_graphs, matrix_from_item


class Structure(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        items = self.canvas.items(event.scenePos(), Qt.ItemSelectionMode.IntersectsItemShape)
        if items == [] or not isinstance(items[0], AlphaAtom):
            return super(Structure, self).mouse_press_event(event)

        molecule_matrix = matrix_from_item(items[0])
        nodes = molecule_matrix[0]
        adjacent = molecule_matrix[1]
        molecule_smiles = Chem.MolToSmiles(mol_from_graphs(nodes, adjacent))
        molecule_dm = Chem.MolFromSmiles(molecule_smiles)
        Chem.Kekulize(molecule_dm)
        if molecule_dm is None or incorrect_valence(molecule_dm):
            image = QImage('resources//stathem.jpg')
            molecule = QGraphicsPixmapItem(QPixmap.fromImage(image))
            for i in molecule_matrix[2]:
                for j in i.lines:
                    self.canvas.removeItem(j)
                self.canvas.removeItem(i)
            self.canvas.addItem(molecule)
            return

        Chem.rdDepictor.Compute2DCoords(molecule_dm)
        zeros = [0, 0]
        atoms = [None for _ in range(molecule_dm.GetNumAtoms())]
        for i, atom in enumerate(molecule_dm.GetAtoms()):
            positions = molecule_dm.GetConformer().GetAtomPosition(i)
            print(atom.GetSymbol(), positions.x, positions.y)
            new_atom = AlphaAtom(atom.GetSymbol())
            if i == 0:
                new_atom.setPos(event.scenePos())
                zeros[0] = positions.x
                zeros[1] = positions.y
            else:
                new_atom.setPos((positions.x - zeros[0]) * 100 + event.scenePos().x(),
                                (positions.y - zeros[1]) * 100 + event.scenePos().y())
            new_atom.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
            new_atom.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)
            self.canvas.addItem(new_atom)
            atoms[i] = new_atom

        for bond in molecule_dm.GetBonds():
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

        for i in molecule_matrix[2]:
            for j in i.lines:
                self.canvas.removeItem(j)
            self.canvas.removeItem(i)
        self.canvas.selectedItems()


    @property
    def asset(self) -> str:
        return 'structure'
