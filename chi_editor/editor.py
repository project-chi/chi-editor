import sys

from PyQt6.QtCore import Qt

from canvas import Canvas
from menu_bar import MenuBar
from tool_bar import ToolBar
from PyQt6.QtWidgets import QApplication, QMainWindow
from PyQt6.QtGui import QIcon


class Editor(QMainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)
        menu_bar = MenuBar()
        tool_bar = ToolBar()
        canvas = Canvas()
        self.setWindowTitle("Project Chi")
        # to be changed to relative dimensions or whatever
        self.resize(400, 200)
        self.setWindowIcon(QIcon("..\\resources\\ProjectChi.png"))
        self.setMenuBar(menu_bar)
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, tool_bar)
        self.setCentralWidget(canvas)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Editor()
    win.show()
    sys.exit(app.exec())
