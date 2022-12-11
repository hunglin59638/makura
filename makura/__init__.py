import importlib.metadata
from pathlib import Path

_DISTRIBUTION_METADATA = importlib.metadata.metadata("makura")
__version__ = _DISTRIBUTION_METADATA["version"]
PKG_PATH = Path(__file__).parent
