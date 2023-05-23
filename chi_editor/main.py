import sys

from PyQt6.QtWidgets import QApplication

from chi_editor.editor import Editor


def main():
    application = QApplication(sys.argv)
    window = Editor()
    window.show()
    sys.exit(application.exec())
