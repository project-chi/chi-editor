import sys

from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QIcon, QPainter
from PyQt6.QtWidgets import QApplication, QMainWindow, QGraphicsScene, QGraphicsView, QGraphicsItem, QButtonGroup, QSizePolicy

from .canvas import Canvas
from .menu_bar import MenuBar
from .toolbar import create_toolbar
from .experiment_ui import Ui_Dialog


class Editor(QMainWindow, Ui_Dialog):
    def __init__(self, parent=None):
        super(Editor, self).__init__(parent)
        self.menu_bar = MenuBar()
        self.tool_bar = create_toolbar()

        self.scene = Canvas(self)
        self.scene.setSceneRect(QRectF(self.geometry()))

        self.graphicsView = QGraphicsView(self)
        self.graphicsView.setScene(self.scene)
        self.graphicsView.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.graphicsView.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        for action in self.tool_bar.actions():
            if action.text() == "Text":
                action.triggered.connect(self.scene.setText)
            else:
                action.triggered.connect(self.scene.setArrow)


        self.setWindowTitle("Project Chi")
        # to be changed to relative dimensions or whatever
        self.resize(400, 200)
        self.setWindowIcon(QIcon("../resources/project-chi.png"))
        self.setMenuBar(self.menu_bar)
        self.addToolBar(Qt.ToolBarArea.LeftToolBarArea, self.tool_bar)
        self.setCentralWidget(self.graphicsView)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Editor()
    win.show()
    sys.exit(app.exec())
