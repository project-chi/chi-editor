from PyQt6 import QtGui
from PyQt6.QtWidgets import QApplication, QMainWindow, QGraphicsScene, QGraphicsView, QGraphicsItem
from PyQt6.QtGui import QPen, QBrush, QColor, QFont

import sys


class Window(QMainWindow):
    def __init__(self):
        super().__init__()
        self.title = "PyQt5 QGraphicView"
        self.top = 200
        self.left = 500
        self.width = 600
        self.height = 500
        self.InitWindow()

    def InitWindow(self):
        self.setWindowIcon(QtGui.QIcon("icon.png"))
        self.setWindowTitle(self.title)
        self.setGeometry(self.left, self.top, self.width, self.height)
        self.createGraphicView()
        self.show()

    def createGraphicView(self):
        self.scene = QGraphicsScene()
        self.greenBrush = QBrush(QColor("green"))
        self.grayBrush = QBrush(QColor("grey"))
        self.pen = QPen(QColor("red"))
        graphicView = QGraphicsView(self.scene, self)
        graphicView.setGeometry(0, 0, 600, 500)
        self.shapes()

    def shapes(self):
        ellipse = self.scene.addEllipse(20, 20, 200, 200, self.pen, self.greenBrush)
        rect = self.scene.addRect(-100, -100, 200, 200, self.pen, self.grayBrush)
        text = self.scene.addText("sadasdasd", QFont("Helvetica"))
        text.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        ellipse.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        rect.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsMovable)
        ellipse.setFlag(QGraphicsItem.GraphicsItemFlag.ItemIsSelectable)

App = QApplication(sys.argv)
window = Window()
sys.exit(App.exec())
