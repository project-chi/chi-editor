import datamol.viz
from PyQt6.QtSvgWidgets import QGraphicsSvgItem
from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem, QGraphicsTextItem
from rdkit import Chem
from datamol import viz, to_mol

from ...bases.tool import Tool
from ...constants import RESOURCES
from ...playground import MolFromGraphs


class Eraser(Tool):
    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        self.canvas.clear()


    @property
    def asset(self) -> str:
        return 'eraser'
