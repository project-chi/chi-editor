from PyQt6.QtGui import QIcon, QActionGroup
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

TOOLS = [
    'arrow', 'hand', 'bond', 'structure',
    'atom', 'block', 'reaction', 'text', 'eraser'
]


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
        for tool in TOOLS:
            btn = getattr(self, '%s_button' % tool)
            group_buttons.addAction(btn)
            self.addAction(btn)
            btn.setIcon(QIcon('../resources/%s.png' % tool))

        self.setMovable(False)
        self.setStyleSheet("""QToolBar { background-color: rgb(212, 204, 234); }""")
