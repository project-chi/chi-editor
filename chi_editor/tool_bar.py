from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QToolBar
from arrow_button import ArrowButton
from bond_button import BondButton
from structure_button import StructureButton
from hand_button import HandButton
from atom_button import AtomButton
from block_button import BlockButton
from reaction_button import ReactionButton
from text_button import TextButton
from eraser_button import EraserButton


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

        self.addAction(self.arrow_button)
        self.addAction(self.hand_button)
        self.addAction(self.bond_button)
        self.addAction(self.structure_button)
        self.addAction(self.atom_button)
        self.addAction(self.block_button)
        self.addAction(self.reaction_button)
        self.addAction(self.text_button)
        self.addAction(self.eraser_button)
        self.setMovable(False)
        self.setStyleSheet("""QToolBar { background-color: rgb(212, 204, 234); }""")
