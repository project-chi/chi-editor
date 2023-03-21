from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QGraphicsItem

from ....bases.line import Line
from ....chem_bonds.double_bond import DoubleBond
from ..bond import Bond


class CreateDoubleBond(Bond):

    def get_line(self, start_atom: QGraphicsItem, mouse_pos: QPointF) -> Line:
        return DoubleBond(start_atom, mouse_pos)

    @property
    def asset(self) -> str:
        return "bond2"
