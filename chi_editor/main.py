import sys

from PyQt6.QtWidgets import QApplication

from .editor import Editor


def main():
    app = QApplication(sys.argv)
    win = Editor()
    win.show()
    sys.exit(app.exec())
