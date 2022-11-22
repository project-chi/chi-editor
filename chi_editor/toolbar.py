from enum import Enum

from PyQt6.QtGui import QIcon, QActionGroup
from PyQt6.QtWidgets import QToolBar

from .toolbar_buttons import *


class Tools(str, Enum):
    Arrow = 'arrow'
    Hand = 'hand'
    Bond = 'bond'
    Structure = 'structure'
    Atom = 'atom'
    Block = 'block'
    Reaction = 'reaction'
    Text = 'text'
    Eraser = 'eraser'


class ToolBar(QToolBar):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.arrow_button = ArrowButton()
        self.hand_button = HandButton()
        self.bond_button = BondButton()
        self.structure_button = StructureButton()
        self.atom_button = AtomButton()
        self.block_button = BlockButton()
        self.reaction_button = ReactionButton()
        self.text_button = TextButton()
        self.eraser_button = EraserButton()

        group_buttons = QActionGroup(self)
        group_buttons.setExclusive(True)
        for tool in Tools:
            btn = getattr(self, '%s_button' % tool)
            group_buttons.addAction(btn)
            self.addAction(btn)
            btn.setIcon(QIcon('../resources/%s.png' % tool))

        self.setMovable(False)
        self.setStyleSheet("""QToolBar { background-color: rgb(212, 204, 234); }""")
