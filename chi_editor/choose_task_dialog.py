from PyQt6.QtWidgets import QDialog, QTreeView, QSizePolicy, QVBoxLayout, QHBoxLayout, QPushButton
from PyQt6.QtCore import Qt


class ChooseTaskDialog(QDialog):
    # View that holds all the tasks links
    view: QTreeView

    # View's layout to hold buttons on top of tasks
    view_layout: QHBoxLayout

    # Buttons to manipulate items
    accept_button: QPushButton

    # Layout that holds view to make it expandable
    layout: QVBoxLayout

    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.setWindowTitle("Choose a task")

        # View init
        self.view = QTreeView(self)
        self.view.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.MinimumExpanding)

        # Main layout
        self.layout = QVBoxLayout(self)
        self.layout.setContentsMargins(0, 2, 0, 0)
        self.layout.addWidget(self.view)

        # View layout
        self.view_layout = QHBoxLayout(self.view)
        self.view_layout.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignBottom)
        self.view_layout.setContentsMargins(0, 0, 2, 2)

        # Buttons
        self.accept_button = QPushButton("Choose task")
        self.accept_button.setFixedSize(self.accept_button.sizeHint())  # sizeHint() is minimal size to fit the text
        self.accept_button.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

        self.view_layout.addWidget(self.accept_button)
