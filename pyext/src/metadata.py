"""@namespace IMP.pmi.metadata
Classes for attaching metadata to PMI objects.
"""

from __future__ import print_function, division
from IMP.pmi.tools import OrderedDict
import os

class Metadata(object):
    """Base class for all metadata"""
    pass


class RootMetadata(Metadata):
    """Metadata that only makes sense for the top-level PMI object."""
    pass


class PythonScript(RootMetadata):
    """A Python script used as part of the modeling."""
    def __init__(self, location):
        self.location = location


class ChimeraXCommandScript(RootMetadata):
    """A ChimeraX command script used to visualize the model."""
    def __init__(self, location):
        self.location = location
