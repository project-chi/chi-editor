import sys

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication, QMainWindow

from .canvas import Canvas
from .menu_bar import MenuBar
from .toolbar import create_toolbar


class Editor(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        self.menu_bar = MenuBar()
        self.tool_bar = create_toolbar()
        self.canvas = Canvas()
        self.setWindowTitle("Project Chi")
        # to be changed to relative dimensions or whatever
        self.resize(400, 200)
        self.setWindowIcon(QIcon("..\\resources\\ProjectChi.png"))
        self.setMenuBar(self.menu_bar)
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.tool_bar)
        self.setCentralWidget(self.canvas)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Editor()
    win.show()
    sys.exit(app.exec())
