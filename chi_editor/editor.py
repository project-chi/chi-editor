import sys

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QIcon
from PyQt6.QtWidgets import QApplication, QMainWindow, QGraphicsScene, QGraphicsView, QGraphicsItem, QButtonGroup

from .canvas import Canvas
from .menu_bar import MenuBar
from .toolbar import create_toolbar
from .experiment_ui import Ui_Dialog


class Editor(QMainWindow, Ui_Dialog):
    def __init__(self, parent=None):
        super(Editor, self).__init__(parent)
        self.menu_bar = MenuBar()
        self.tool_bar = create_toolbar()

        self.setupUi(self)
        self.scene = Canvas(self)
        self.graphicsView.setScene(self.scene)
        self.graphicsView.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)

        self.setWindowTitle("Project Chi")
        # to be changed to relative dimensions or whatever
        self.resize(400, 200)
        self.setWindowIcon(QIcon("..\\resources\\ProjectChi.png"))
        self.setMenuBar(self.menu_bar)
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.tool_bar)
        for action in self.tool_bar.actions():
            if action.text() == "Text":
                action.triggered.connect(self.scene.setText)
            else:
                action.triggered.connect(self.scene.setArrow)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Editor()
    win.show()
    sys.exit(app.exec())
