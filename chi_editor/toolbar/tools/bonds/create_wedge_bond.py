from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QGraphicsItem

from ....bases.line import Line
from ....chem_bonds.wedge_bond import WedgeBond
from ..bond import Bond


class CreateWedgeBond(Bond):

    def get_line(self, start_atom: QGraphicsItem, mouse_pos: QPointF) -> Line:
        return WedgeBond(start_atom, mouse_pos)

    @property
    def picture(self) -> str:
        return "bond1"
