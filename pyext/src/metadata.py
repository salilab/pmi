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


class Software(RootMetadata):
    """Software (other than IMP) used as part of the modeling protocol."""
    def __init__(self, name, classification, description, url, type='program',
                 version=None):
        self.name = name
        self.classification = classification
        self.description = description
        self.url = url
        self.type = type
        self.version = version


class Citation(RootMetadata):
    """A publication that describes the modeling."""
    def __init__(self, pmid, title, journal, volume, page_range, year, authors,
                 doi):
        self.title, self.journal, self.volume = title, journal, volume
        self.page_range, self.year = page_range, year
        self.pmid, self.authors, self.doi = pmid, authors, doi


class PythonScript(RootMetadata):
    """A Python script used as part of the modeling."""
    def __init__(self, location):
        self.location = location


class ChimeraXCommandScript(RootMetadata):
    """A ChimeraX command script used to visualize the model."""
    def __init__(self, location):
        self.location = location


class Dataset(Metadata):
    """A set of input data, for example, a crystal structure or EM map."""

    _eq_keys = ['location']

    # Datasets compare equal iff they are the same class and have the
    # same attributes
    def _eq_vals(self):
        return tuple([self.__class__]
                     + [getattr(self, x) for x in self._eq_keys])
    def __eq__(self, other):
        return self._eq_vals() == other._eq_vals()
    def __hash__(self):
        return hash(self._eq_vals())

    _data_type = 'unspecified'
    def __init__(self, location):
        self.location = location
        self._parents = OrderedDict()

    def add_parent(self, dataset):
        """Add another Dataset from which this one was derived.
           For example, a 3D EM map may be derived from a set of 2D images."""
        self._parents[dataset] = None

    def add_primary(self, dataset):
        """Add another Dataset from which the ultimate parent of this one
           was derived."""
        if len(self._parents) == 0:
            self.add_parent(dataset)
        elif len(self._parents) == 1:
            list(self._parents.keys())[0].add_parent(dataset)
        else:
            raise ValueError("This dataset has multiple parents - don't "
                             "know which one to add to")

class CXMSDataset(Dataset):
    """Processed crosslinks from a CX-MS experiment"""
    _data_type = 'CX-MS data'

class MassSpecDataset(Dataset):
    """Raw mass spectrometry files such as peaklists"""
    _data_type = 'Mass Spectrometry data'

class EMDensityDataset(Dataset):
    """A 3D electron microscopy dataset"""
    _data_type = '3DEM volume'

class PDBDataset(Dataset):
    """An experimentally-determined 3D structure as a set of a coordinates,
       usually in a PDB file"""
    _data_type = 'Experimental model'

class ComparativeModelDataset(Dataset):
    """A 3D structure determined by comparative modeling"""
    _data_type = 'Comparative model'

class IntegrativeModelDataset(Dataset):
    """A 3D structure determined by integrative modeling"""
    _data_type = 'Integrative model'

class EMMicrographsDataset(Dataset):
    """Raw 2D electron micrographs"""
    _eq_keys = Dataset._eq_keys + ['number']

    _data_type = 'EM raw micrographs'
    def __init__(self, location, number):
        super(EMMicrographsDataset, self).__init__(location)
        self.number = number

class EM2DClassDataset(Dataset):
    """2DEM class average"""
    _data_type = '2DEM class average'

class SASDataset(Dataset):
    """SAS data"""
    _data_type = 'SAS data'
