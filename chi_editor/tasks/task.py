from rdkit import Chem

from .task_type import TaskType


class Task:
    # Task's title
    _title: str

    # Answer represented as one SMILES string
    _correct_answer: str

    # Task's type
    _type: TaskType

    def __init__(self, title: str, correct_answer: str, task_type: TaskType):
        self._title = title
        self._correct_answer = Chem.CanonSmiles(correct_answer)
        self._type = task_type

    def checkAnswer(self, user_answer: str) -> bool:
        canon_user_answer = Chem.CanonSmiles(user_answer)
        return self._correct_answer == canon_user_answer
