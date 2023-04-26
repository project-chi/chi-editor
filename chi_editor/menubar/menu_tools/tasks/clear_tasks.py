from chi_editor.api.server import default_url, Server
from chi_editor.bases.menu_tool import MenuTool


class ClearTasks(MenuTool):
    def exec(self) -> None:
        server = Server(default_url)
        tasks = server.get_tasks()
        for task in tasks:
            server.delete_task(task)
        pass

    @property
    def asset(self) -> str:
        return ""

    @property
    def text_asset(self) -> str:
        return "Clear tasks"
