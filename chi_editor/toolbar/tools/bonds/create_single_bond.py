from PyQt6.QtWidgets import QGraphicsItem
from PyQt6.QtCore import QPointF

from ..bond import Bond
from ....bases.line import Line
from ....chem_bonds.single_bond import SingleBond


class CreateSingleBond(Bond):

    def get_line(self, start_atom: QGraphicsItem, mouse_pos: QPointF) -> Line:
        return SingleBond(start_atom, mouse_pos)

    @property
    def asset(self) -> str:
        return 'bond1'
