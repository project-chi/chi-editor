from PyQt6.QtCore import QPointF
from PyQt6.QtWidgets import QGraphicsItem

from ....bases.line import Line
from ....chem_bonds.wedge_bond import WedgeBond
from chi_editor.toolbar.tools.bonds.create_bond import CreateBond


class CreateWedgeBond(CreateBond):

    def get_line(self, start_atom: QGraphicsItem, mouse_pos: QPointF) -> Line:
        return WedgeBond(start_atom, mouse_pos)

    @property
    def picture(self) -> str:
        return "bond1"
