from enum import Enum, auto
from uuid import UUID

from pydantic import BaseModel


class Kind(Enum):
    Chain = auto()
    Molecule = auto()
    Reaction = auto()


class Task(BaseModel):
    identifier: UUID
    name: str
    kind: Kind
    problem: str
    initial: str
    solution: str
