from chi_editor.tasks import Task, TaskType


class RandomTask(Task):
    def __init__(self, *args, **kwargs):
        server = Server("http://kapkekes.site:8000")
        super().__init__("Benzene", "Draw a benzene", "c1ccccc1", TaskType.Random)

