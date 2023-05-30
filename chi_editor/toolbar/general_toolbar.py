from typing import TYPE_CHECKING, cast

from PyQt6.QtCore import QSize
from PyQt6.QtGui import QActionGroup
from PyQt6.QtWidgets import QToolBar, QToolButton

from chi_editor.toolbar.tools import tools, menus
from PyQt6.QtGui import QAction

if TYPE_CHECKING:
    from chi_editor.canvas import Canvas


class GeneralToolBar(QToolBar):
    _canvas: "Canvas"
    _action_group: "QActionGroup"

    def __init__(self, *args, canvas: "Canvas", **kwargs) -> "None":
        super().__init__(*args, **kwargs)
        self._canvas = canvas
        self.setStyleSheet("""QToolBar { background-color: rgb(212, 204, 234); }""")
        self.setIconSize(QSize(20, 20))
        self.setMovable(False)

        self.action_group = QActionGroup(self)
        self.action_group.setExclusive(True)

        for Tool in tools:
            tool = Tool(canvas)
            self.addAction(tool)
            self.action_group.addAction(tool)
            self.widgetForAction(tool).setStyleSheet("padding: 5px")

        for menu in menus:
            menu_tools = menu[0]
            menu_button = menu[1]()
            self.addWidget(menu_button)
            menu_button.setStyleSheet("padding: 5px")
            menu_button.setPopupMode(QToolButton.ToolButtonPopupMode.MenuButtonPopup)
            for Tool in menu_tools:
                tool = Tool(canvas)
                menu_button.addAction(tool)
                self.action_group.addAction(tool)
                tool.toggled.connect(menu_button.setChecked)  # focus on menu button when one of its items is chosen
            menu_button.triggered.connect(self.changeAction)

        self.actionTriggered.connect(self.changeAction)

    def changeAction(self, action: "QAction") -> "None":
        self._canvas.current_action = action

    def change_canvas(self, canvas: "Canvas"):
        self._canvas = canvas
        for Action in self.actions():
            Action.change_canvas(canvas)
