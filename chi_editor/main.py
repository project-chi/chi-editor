import sys

from PyQt6.QtWidgets import QApplication

from .editor import Editor


def main():
    application = QApplication(sys.argv)
    window = Editor()
    window.show()
    sys.exit(application.exec())
