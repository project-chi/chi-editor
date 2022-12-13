import datamol.viz
from PyQt6.QtSvgWidgets import QGraphicsSvgItem
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem, QGraphicsTextItem
from rdkit import Chem
from datamol import viz, to_mol

from ...bases.alpha_atom import AlphaAtom
from ...bases.tool import Tool
from ...constants import RESOURCES
from ...playground import mol_from_graphs, matrix_from_item


class Structure(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        selected = self.canvas.selectedItems()
        molecule_matrix = matrix_from_item(selected[0])
        nodes = molecule_matrix[0]
        adjacent = molecule_matrix[1]
        molecule_smiles = Chem.MolToSmiles(mol_from_graphs(nodes, adjacent))
        molecule_dm = to_mol(mol=molecule_smiles)
        viz.to_image(mols=molecule_dm, use_svg=True, outfile=RESOURCES / 'molecule.svg')
        molecule = QGraphicsSvgItem('resources//molecule.svg')
        molecule.setPos(event.scenePos())
        molecule.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)


        self.canvas.addItem(molecule)
        self.canvas.selectedItems()


    @property
    def asset(self) -> str:
        return 'structure'
