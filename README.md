# chi-editor

The molecule editor in the Project Chi software package, currently in early development.

## Preparations

Firstly, ensure that you have an installed Python 3.10. If you don't, get the one, [Python 3.10.8](https://www.python.org/downloads/release/python-3108/) for example. Notice: you need to remember the installation path, so we can pass it to Poetry later!

After that, you need to set up Poetry: you can find all instructions [here](https://python-poetry.org/docs/#installation).

## Setup

Clone the repository and open the terminal inside it, execute:

```shell
poetry env use {PATH-TO-PYTHON-3.10-EXECUTABLE}
poetry install
```

After that, Poetry should create a virtual environment for this project, so you can try to play with it!

## Usage

To start the application, execute:

```shell
poetry run python -m chi_editor
```
