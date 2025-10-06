"""
magpie
======

Provides
  1. A function for exploring rectangular plates with generalised boundary conditions :code:`magpie`
  2. A function for generating a biharmonic sparse matrix :code:`bhmat`
  3. A function for generating impulse responses of a given plate :code:`modal_time_integration`
  4. A function for generating estimating the young's modulus of a given plate :code:`youngcalc`

How to use the documentation
----------------------------


Utilities
---------

__version__
    magpie version string

Viewing documentation using IPython
-----------------------------------


Credits
-------
:Authors:
    Michele Ducceschi
    Matthew Hamilton
    Alexis Mousseau

:Version: 1.0 of 2001/01/01
:Dedication:

"""
# from . import magpie
from .magpie import magpie
# from . import bhmat
from .bhmat import bhmat
# from . import youngcalc
from .youngcalc import youngcalc
# from . import modal_time_integration
from .modal_time_integration import modal_time_integration


def test_sphinx():
    """
    This is a placeholder function that allows testing of the sphinx documentation format.
    This function doesn't do anything, so don't worry.

    The folder tree for setting up looks as follows (using attached attribute
    names rather than paths):

    Formatting code
    ---------------

    .. code-block:: python
        :linenos:

        def test_sphinx():
            pass

        if __name__ == '__main__':
            test_sphinx()

    For more info, see the `Sphinx Docs <https://www.sphinx-doc.org/en/master/usage/restructuredtext/directives.html#directive-code-block>`_

    inline code can be written like :code:`val = number`

    The documentation for roles can be found `here <https://www.sphinx-doc.org/en/master/usage/restructuredtext/roles.html>`_


    Title
    =====

    Titles are underlined (or over- and underlined) with
    a nonalphanumeric character at least as long as the
    text.

    A lone top-level section is lifted up to be the
    document's title.

    Any non-alphanumeric character can be used, but
    Python convention is:

    * ``#`` with overline, for parts
    * ``*`` with overline, for chapters
    * ``=``, for sections
    * ``-``, for subsections
    * ``^``, for subsubsections
    * ``"``, for paragraphs

    =====  =====  ======
       Inputs     Output
    ------------  ------
      A      B    A or B
    =====  =====  ======
    False  False  False
    True   False  True
    False  True   True
    True   True   True
    =====  =====  ======

    [Summary]

    :param [ParamName]: [ParamDescription], defaults to [DefaultParamVal]
    :type [ParamName]: [ParamType](, optional)
    ...
    :raises [ErrorType]: [ErrorDescription]
    ...
    :return: [ReturnDescription]
    :rtype: [ReturnType]

    Credits
    ==========
    :Authors:
        Michele Ducceschi
        Matthew Hamilton
        Alexis Mousseau

    :Version: 1.0 of 2001/01/01
    :Dedication:

    This is a conceptual class representation of a simple BLE device
    (GATT Server). It is essentially an extended combination of the
    :class:`bluepy.btle.Peripheral` and :class:`bluepy.btle.ScanEntry` classes

    :param client: A handle to the :class:`simpleble.SimpleBleClient` client
        object that detected the device
    :type client: class:`simpleble.SimpleBleClient`
    :param addr: Device MAC address, defaults to None
    :type addr: str, optional
    :param addrType: Device address type - one of ADDR_TYPE_PUBLIC or
        ADDR_TYPE_RANDOM, defaults to ADDR_TYPE_PUBLIC
    :type addrType: str, optional
    :param iface: Bluetooth interface number (0 = /dev/hci0) used for the
        connection, defaults to 0
    :type iface: int, optional
    :param data: A list of tuples (adtype, description, value) containing the
        AD type code, human-readable description and value for all available
        advertising data items, defaults to None
    :type data: list, optional
    :param rssi: Received Signal Strength Indication for the last received
        broadcast from the device. This is an integer value measured in dB,
        where 0 dB is the maximum (theoretical) signal strength, and more
        negative numbers indicate a weaker signal, defaults to 0
    :type rssi: int, optional
    :param connectable: `True` if the device supports connections, and `False`
        otherwise (typically used for advertising ‘beacons’).,
        defaults to `False`
    :type connectable: bool, optional
    :param updateCount: Integer count of the number of advertising packets
        received from the device so far, defaults to 0
    :type updateCount: int, optional
    """
    pass


__all__ = [
    "magpie",
    "bhmat",
    "youngcalc",
    "modal_time_integration",
]
