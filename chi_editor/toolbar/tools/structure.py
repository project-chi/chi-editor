import datamol.viz
import rdkit.Chem.rdDepictor
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtSvgWidgets import QGraphicsSvgItem
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem, QGraphicsTextItem, QGraphicsPixmapItem
from PyQt6.QtCore import Qt
from rdkit import Chem
from datamol import viz, to_mol, incorrect_valence

from ... import bases
from ...bases.alpha_atom import AlphaAtom
from ...bases.line import Line
from ...bases.tool import Tool
from ...constants import RESOURCES
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
        molecule_dm = to_mol(mol=molecule_smiles)
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
            bonds = atom.GetBonds()
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
            bond_type = bond.GetBondType()
            new_bond = Line(atoms[start_position], atoms[end_position])
            atoms[start_position].add_line(new_bond)
            atoms[end_position].add_line(new_bond)
            self.canvas.addItem(new_bond)

        # viz.to_image(mols=molecule_dm, use_svg=True, outfile=RESOURCES / 'molecule.svg')
        # molecule = QGraphicsSvgItem('resources//molecule.svg')
        # molecule.setPos(event.scenePos())
        # molecule.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)

        for i in molecule_matrix[2]:
            for j in i.lines:
                self.canvas.removeItem(j)
            self.canvas.removeItem(i)
        #self.canvas.addItem(molecule)
        self.canvas.selectedItems()


    @property
    def asset(self) -> str:
        return 'structure'
