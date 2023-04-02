from enum import Enum, auto
from uuid import UUID

from pydantic import BaseModel


class Kind(Enum):
    Chain = auto()
    Molecule = auto()
    Reaction = auto()


class Task(BaseModel):
    identifier: UUID | None
    name: str
    kind: Kind
    problem: str
    initial: str | None
    solution: str
