# coding=utf-8

__doc__ = """\
Module for opening and dealing with files (ixo stands for IO), including:
- Fast and efficient parsing of binary files from the Foldometer software
- Parsing of txt files from the old setup
- Data preparation and conversion to force and extension from raw data
"""

import os
import glob

SOURCE_FILES = glob.glob(os.path.dirname(__file__) + "/*.py")
__all__ = [os.path.basename(f)[: -3] for f in SOURCE_FILES]



from .binary import read_file
from .data_conversion import process_file, process_data
from .old_setup import read_file_old_setup