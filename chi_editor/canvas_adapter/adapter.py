from PyQt6.QtWidgets import QWidget
from entity import Entity


class Adapter:
    def __init__(self, base: QWidget) -> None:
        self.base = base
        self.entities = {}

    def add(self, new_entity: Entity) -> None:
        self.entities[new_entity] = {"scaling": 0.0}
