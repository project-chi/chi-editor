from chi_editor.tasks.task import Task
from chi_editor.tasks.task_type import TaskType

alkanes_list: tuple[list[Task], TaskType] = (
    [
        Task("Methane", "Draw a methane", "C", TaskType.ALKANES),
        Task("Ethane", "Draw an ethane", "CC", TaskType.ALKANES)
    ],
    TaskType.ALKANES
)

alcohols_list: tuple[list[Task], TaskType] = (
    [
        Task("Ethanol", "Draw an ethanol", "CCO", TaskType.ALCOHOLS),
        Task("1-Butanol", "Draw a 1-butanol", "CCCCO", TaskType.ALCOHOLS),
    ],
    TaskType.ALCOHOLS
)

aromatic_list: tuple[list[Task], TaskType] = (
    [
        Task("Benzene", "Draw a benzene", "c1ccccc1", TaskType.AROMATIC),
        Task("Toluene", "Draw a toluene", "Cc1ccccc1", TaskType.AROMATIC),
    ],
    TaskType.AROMATIC
)

tasks_list: list[tuple[list[Task], TaskType]] = [
    alkanes_list,
    alcohols_list,
    aromatic_list,
]
