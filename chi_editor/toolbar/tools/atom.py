from PyQt6.QtWidgets import QGraphicsSceneMouseEvent, QGraphicsItem

from chi_editor.bases.molecule.molecule import MoleculeDrawer
from ...bases.tool import Tool
from ...bases.alpha_atom import AlphaAtom


class Atom(Tool):
    _element: str

    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        new_atom = AlphaAtom(self._element)
        new_atom.setPos(event.scenePos() - new_atom.sceneBoundingRect().center())
        new_atom.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        new_atom.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)

        new_molecule = new_atom.molecule.molecule_drawer
        new_molecule.setPos(new_atom.pos())

        self.canvas.addItem(new_atom)
        self.canvas.addItem(new_molecule)

    @property
    def asset(self) -> str:
        return "carbon"
