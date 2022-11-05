from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QToolBar
from tool_bar_buttons.arrow_button import ArrowButton
from tool_bar_buttons.bond_button import BondButton
from tool_bar_buttons.structure_button import StructureButton
from tool_bar_buttons.hand_button import HandButton
from tool_bar_buttons.atom_button import AtomButton
from tool_bar_buttons.block_button import BlockButton
from tool_bar_buttons.reaction_button import ReactionButton
from tool_bar_buttons.text_button import TextButton
from tool_bar_buttons.eraser_button import EraserButton


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

        self._set_actions()
        self._set_icons()
        self.setMovable(False)
        self.setStyleSheet("""QToolBar { background-color: rgb(212, 204, 234); }""")

    def _set_actions(self):
        self.addAction(self.arrow_button)
        self.addAction(self.hand_button)
        self.addAction(self.bond_button)
        self.addAction(self.structure_button)
        self.addAction(self.atom_button)
        self.addAction(self.block_button)
        self.addAction(self.reaction_button)
        self.addAction(self.text_button)
        self.addAction(self.eraser_button)

    def _set_icons(self):
        self.arrow_button.setIcon(QIcon("../resources/arrow.png"))
        self.hand_button.setIcon(QIcon("../resources/hand.png"))
        self.bond_button.setIcon(QIcon("../resources/bond.png"))
        self.structure_button.setIcon(QIcon("../resources/structure.png"))
        self.atom_button.setIcon(QIcon("../resources/atom.png"))
        self.block_button.setIcon(QIcon("../resources/block.png"))
        self.reaction_button.setIcon(QIcon("../resources/reaction.png"))
        self.text_button.setIcon(QIcon("../resources/text.png"))
        self.eraser_button.setIcon(QIcon("../resources/eraser.png"))
