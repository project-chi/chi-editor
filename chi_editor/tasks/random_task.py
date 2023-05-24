from random import choice

from chi_editor.api.server import Server, default_url
from chi_editor.tasks import Task


class RandomTask(Task):
    def __init__(self, *args, **kwargs) -> None:
        server = Server(default_url)
        task_ids = server.get_tasks_identifiers()
        task = server.get_task(choice(task_ids))
        super().__init__(task.name, task.problem, task.solution, task.kind)
