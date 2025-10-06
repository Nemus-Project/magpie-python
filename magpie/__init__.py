"""
magpie
======

Version: |release|
------------------

Provides
  1. A function for exploring rectangular plates with generalised boundary conditions :code:`magpie`
  2. A function for generating a biharmonic sparse matrix :code:`bhmat`
  3. A function for generating impulse responses of a given plate :code:`modal_time_integration`
  4. A function for generating estimating the young's modulus of a given plate :code:`youngcalc`

Credits
-------
:Authors:
    Michele Ducceschi
    Matthew Hamilton
    Alexis Mousseau
"""
# from . import magpie
from .magpie import magpie
# from . import bhmat
from .bhmat import bhmat
# from . import youngcalc
from .youngcalc import youngcalc
# from . import modal_time_integration
from .modal_time_integration import modal_time_integration

__version__ = "0.0.5"

__all__ = [
    "magpie",
    "bhmat",
    "youngcalc",
    "modal_time_integration",
]
