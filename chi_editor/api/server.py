from uuid import UUID

from pydantic import ValidationError
from requests import get, post

from chi_editor.api.task import Kind, Task


class Server:
    _url: str

    def __init__(self, url: str) -> None:
        response = get(f"{url}/docs")
        if not response.ok:
            raise ValueError("can't find server at this URL")

        self._url = url

    def get_tasks_identifiers(self) -> list[UUID]:
        url = f"{self._url}/tasks"
        response = get(url)
        return [UUID(identifier) for identifier in response.json()]

    def get_task(self, identifier: str | UUID) -> Task | None:
        if isinstance(identifier, str):
            identifier = UUID(identifier)

        if not isinstance(identifier, UUID):
            raise TypeError(
                f"identifier should be a str or UUID instance, not {type(identifier)}"
            )

        url = f"{self._url}/task/{identifier}"
        data = get(url).json()
        if data is None:
            return None

        return Task.parse_obj(data)

    def create_task(
        self,
        name: str,
        kind: Kind,
        problem: str,
        solution: str,
        initial: str = "",
    ) -> Task | None:
        try:
            task = Task(
                name=name,
                kind=kind,
                problem=problem,
                initial=initial,
                solution=solution,
            )
        except ValidationError:
            raise TypeError("check types of args")

        url = f"{self._url}/task"
        data = task.dict(exclude={"identifier"})
        data["kind"] = kind.value
        response = post(
            url, json=data,
            headers={"accept": "application/json", "Content-Type": "application/json"}
        )
        if not response.ok:
            print(response.json())
            return None

        return Task.parse_obj(response.json())

