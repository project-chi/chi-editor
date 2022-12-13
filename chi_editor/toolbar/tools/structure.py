import datamol.viz
from PyQt6.QtGui import QImage, QPixmap
from PyQt6.QtSvgWidgets import QGraphicsSvgItem
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem, QGraphicsTextItem, QGraphicsPixmapItem
from PyQt6.QtCore import Qt
from rdkit import Chem
from datamol import viz, to_mol, incorrect_valence

from ...bases.alpha_atom import AlphaAtom
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
        viz.to_image(mols=molecule_dm, use_svg=True, outfile=RESOURCES / 'molecule.svg')
        molecule = QGraphicsSvgItem('resources//molecule.svg')
        molecule.setPos(event.scenePos())
        molecule.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)

        for i in molecule_matrix[2]:
            for j in i.lines:
                self.canvas.removeItem(j)
            self.canvas.removeItem(i)
        self.canvas.addItem(molecule)
        self.canvas.selectedItems()


    @property
    def asset(self) -> str:
        return 'structure'
