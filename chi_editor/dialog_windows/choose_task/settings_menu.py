from PyQt6.QtWidgets import QMenu


class SettingsMenu(QMenu):
    def __init__(self, *args, **kwargs) -> None:
        super().__init__(*args, **kwargs)
        self.addAction("privet")
