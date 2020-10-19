import pkg_resources
from pkg_resources import DistributionNotFound, VersionConflict
from . import _ROOT
from os.path import join

if __name__ == "__main__":
    with open(join(_ROOT, "..", "requirements.txt")) as f:
        dependencies = f.read().splitlines()

    # here, if a dependency is not met, a DistributionNotFound or VersionConflict
    # exception is thrown.
    pkg_resources.require(dependencies)
    print("All your dependencies have been installed.")
