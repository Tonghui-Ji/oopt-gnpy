import os,sys
from pathlib import Path

class GlobalControl(object):
    _instance = None

    def __new__(cls, *args, **kw):
        if cls._instance is None:
            cls._instance = object.__new__(cls, *args, **kw)
        return cls._instance

    def __init__(self, bool_fig=True, log_level='info'):
        self.bool_fig = bool_fig
        self.log_level = log_level

        self.opticommpy_path = Path(os.path.dirname(os.path.dirname(os.path.abspath(__file__ if '__file__' in globals() else sys.executable))))

        pass
