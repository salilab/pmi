"""@namespace IMP.pmi.metadata
Classes for attaching metadata to PMI objects.
"""

from __future__ import print_function, division
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
    def __init__(self, pmid, title, journal, volume, page_range, year, authors):
        self.title, self.journal, self.volume = title, journal, volume
        self.page_range, self.year = page_range, year
        self.pmid, self.authors = pmid, authors


class Repository(Metadata):
    """A repository containing modeling files.
       This can be used if the PMI script is part of a repository, which
       has been archived somewhere with a DOI.
       This can be used to construct permanent references to files
       used in this modeling, even if they haven't been uploaded to
       a database such as PDB or EMDB."""

    def __init__(self, doi, root):
        """Constructor.
           @param doi the Digital Object Identifer for the repository.
           @param root the relative path to the top-level directory
                  of the repository from the working directory of the script.
        """
        self.doi, self._root = doi, root

    def get_path(self, fname):
        """Return a path relative to the top of the repository"""
        return os.path.relpath(fname, self._root)
