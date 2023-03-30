from __future__ import annotations

from PyQt6.QtCore import QPointF, Qt
from PyQt6.QtWidgets import (
    QGraphicsEllipseItem,
    QGraphicsItem,
    QGraphicsSceneMouseEvent,
)

from ...bases.alpha_atom import AlphaAtom
from ...bases.line import Line
from ...bases.tool import Tool


class Bond(Tool):
    startItem: AlphaAtom = None
    bond: Line = None

    def atom_at(self, pos: QPointF) -> AlphaAtom | None:
        for item in self.canvas.items(pos, Qt.ItemSelectionMode.IntersectsItemShape):
            if isinstance(item, AlphaAtom):
                return item
        return None

    def mouse_press_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            atom = self.atom_at(event.scenePos())
            if atom is not None:
                self.startItem = atom
                self.bond = self.get_line(atom, event.scenePos())
                self.canvas.addItem(self.bond)

    def mouse_move_event(self, event: QGraphicsSceneMouseEvent) -> None:
        if self.bond is not None:
            atom = self.atom_at(event.scenePos())
            if atom is not None and atom != self.startItem:
                new_end = atom.sceneBoundingRect().center()  # lock on atom
            else:
                new_end = event.scenePos()

            mouse_item = QGraphicsEllipseItem(0, 0, 0, 0)
            mouse_item.setPos(new_end)
            self.bond.update_pixmap(mouse_item, following_mouse=True)

    def mouse_release_event(self, event) -> None:
        if self.bond is None:
            return

        end_atom = self.atom_at(event.scenePos())
        if end_atom is None or end_atom == self.startItem:
            self.canvas.removeItem(self.bond)
            return

        self.bond.set_v2(end_atom)
        only_one_line_between = self.startItem.add_line(self.bond)
        if not only_one_line_between:
            self.canvas.removeItem(self.bond)
        else:  # if line didn't exist before, we add it
            end_atom.add_line(self.bond)
        self.startItem.molecule.update_atoms()

    # should be @property
    def get_line(self, start_atom: QGraphicsItem, mouse_pos: QPointF) -> Line:
        pass

    @property
    def picture(self) -> str:
        return "bond"
