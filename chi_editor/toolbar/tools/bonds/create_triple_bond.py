from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QGraphicsItem

from ....bases.line import Line
from ....chem_bonds.triple_bond import TripleBond
from ..bond import Bond


class CreateTripleBond(Bond):

    def get_line(self, start_atom: QGraphicsItem, mouse_pos: QPointF) -> Line:
        return TripleBond(start_atom, mouse_pos)

    @property
    def asset(self) -> str:
        return "bond3"
