import os

_ROOT = os.path.abspath(os.path.dirname(__file__))


def get_data(path):
    return os.path.join(_ROOT, 'data', path)


def get_julia(path):
    """
    get the julia code directory.
    """
    return os.path.join(_ROOT, 'julia', path)


def get_requirements():
    return os.path.join(_ROOT, "..", "requirements.txt")


__version__ = "0.1.0"
