import datamol.viz
from PyQt6.QtSvgWidgets import QGraphicsSvgItem
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem, QGraphicsTextItem
from rdkit import Chem
from datamol import viz, to_mol

from ...bases.tool import Tool
from ...constants import RESOURCES
from ...playground import MolFromGraphs


class Structure(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        nodes = ["C", "C", "N", "C", "C", "C", "C", "O"]
        adjacent = [[0., 1., 0., 0., 0., 0., 0., 0.],
                            [1., 0., 2., 0., 0., 0., 0., 1.],
                            [0., 2., 0., 1., 0., 0., 0., 0.],
                            [0., 0., 1., 0., 1., 0., 0., 0.],
                            [0., 0., 0., 1., 0., 1., 0., 0.],
                            [0., 0., 0., 0., 1., 0., 1., 1.],
                            [0., 0., 0., 0., 0., 1., 0., 0.],
                            [0., 1., 0., 0., 0., 1., 0., 0.]]
        molecule_smiles = Chem.MolToSmiles(MolFromGraphs(nodes, adjacent))
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
