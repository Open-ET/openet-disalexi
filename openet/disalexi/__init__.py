from .disalexi import Image
from .disalexi import ta_mosaic_interpolate
from .disalexi import ta_mosaic_min_bias
from .collection import Collection
from .landsat import Landsat_C02_L2

MODEL_NAME = 'DISALEXI'

from importlib import metadata

__version__ = metadata.version(__package__.replace('.', '-') or __name__.replace('.', '-'))
# __version__ = metadata.version('openet-ssebop')
