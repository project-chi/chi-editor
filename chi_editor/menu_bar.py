from PyQt6.QtWidgets import QMenuBar


class MenuBar(QMenuBar):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._create_menus()
        self.setStyleSheet("""QMenuBar { background-color: rgb(212, 204, 234); }""")

    def _create_menus(self):
        file_menu = self.addMenu("File")
        edit_menu = self.addMenu("Edit")
        view_menu = self.addMenu("View")
        help_menu = self.addMenu("Help")
