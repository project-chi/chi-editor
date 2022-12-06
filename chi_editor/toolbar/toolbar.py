from PyQt6.QtGui import QActionGroup
from PyQt6.QtWidgets import QToolBar

from .button_factory import Tools, create_buttons


def create_toolbar() -> QToolBar:
    toolbar = QToolBar()
    toolbar.setStyleSheet("""QToolBar { background-color: rgb(212, 204, 234); }""")
    toolbar.setMovable(False)

    action_group = QActionGroup(toolbar)
    action_group.setExclusive(True)

    for tool_button in create_buttons(Tools):
        toolbar.addAction(tool_button)
        action_group.addAction(tool_button)
        action_group.setExclusive(True)

    return toolbar
