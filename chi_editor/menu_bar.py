from PyQt6.QtWidgets import QMenuBar

from menu_buttons.file_menu import FileMenu
from menu_buttons.edit_menu import EditMenu
from menu_buttons.view_menu import ViewMenu
from menu_buttons.help_menu import HelpMenu


class MenuBar(QMenuBar):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._create_menus()
        self.setStyleSheet("""QMenuBar { background-color: rgb(212, 204, 234); }""")

    def _create_menus(self):
        self.edit_menu = EditMenu()
        self.view_menu = ViewMenu()
        self.help_menu = HelpMenu()
        self.addMenu(FileMenu())
        self.addMenu(self.edit_menu)
        self.addMenu(self.view_menu)
        self.addMenu(self.help_menu)

