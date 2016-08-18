"""@namespace IMP.pmi.mmcif
   @brief Support for the mmCIF file format.
"""

from __future__ import print_function
import IMP.atom
import IMP.pmi.representation
import IMP.pmi.tools
from IMP.pmi.tools import OrderedDict
import IMP.pmi.output
import IMP.pmi.metadata
import re
import sys
import os
import operator
import textwrap

class _LineWriter(object):
    def __init__(self, writer, line_len=80, multi_line_len=70):
        self.writer = writer
        self.line_len = line_len
        self.multi_line_len = multi_line_len
        self.column = 0
    def write(self, val):
        if isinstance(val, str) and len(val) > self.multi_line_len:
            self.writer.fh.write("\n;")
            for i in range(0, len(val), self.multi_line_len):
                self.writer.fh.write(val[i:i+self.multi_line_len])
                self.writer.fh.write("\n")
            self.writer.fh.write(";\n")
            self.column = 0
            return
        val = self.writer._repr(val)
        if self.column > 0:
            if self.column + len(val) + 1 > self.line_len:
                self.writer.fh.write("\n")
                self.column = 0
            else:
                self.writer.fh.write(" ")
                self.column += 1
        self.writer.fh.write(val)
        self.column += len(val)


class CifCategoryWriter(object):
    def __init__(self, writer, category):
        self.writer = writer
        self.category = category
    def write(self, **kwargs):
        self.writer._write(self.category, kwargs)
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        pass


class CifLoopWriter(object):
    def __init__(self, writer, category, keys):
        self.writer = writer
        self.category = category
        self.keys = keys
        self._empty_loop = True
    def write(self, **kwargs):
        if self._empty_loop:
            f = self.writer.fh
            f.write("#\nloop_\n")
            for k in self.keys:
                f.write("%s.%s\n" % (self.category, k))
            self._empty_loop = False
        l = _LineWriter(self.writer)
        for k in self.keys:
            l.write(kwargs.get(k, self.writer.omitted))
        self.writer.fh.write("\n")
    def __enter__(self):
        return self
    def __exit__(self, exc_type, exc_value, traceback):
        if not self._empty_loop:
            self.writer.fh.write("#\n")


class CifWriter(object):
    omitted = '.'
    unknown = '?'
    _boolmap = {False: 'NO', True: 'YES'}

    def __init__(self, fh):
        self.fh = fh
    def category(self, category):
        return CifCategoryWriter(self, category)
    def loop(self, category, keys):
        return CifLoopWriter(self, category, keys)
    def write_comment(self, comment):
        for line in textwrap.wrap(comment, 78):
            self.fh.write('# ' + line + '\n')
    def _write(self, category, kwargs):
        for key in kwargs:
            self.fh.write("%s.%s %s\n" % (category, key,
                                          self._repr(kwargs[key])))
    def _repr(self, obj):
        if isinstance(obj, str) and '"' not in obj \
           and "'" not in obj and " " not in obj:
            return obj
        elif isinstance(obj, float):
            return "%.3f" % obj
        elif isinstance(obj, bool):
            return self._boolmap[obj]
        else:
            return repr(obj)

def _get_by_residue(p):
    """Determine whether the given particle represents a specific residue
       or a more coarse-grained object."""
    return IMP.atom.Residue.get_is_setup(p) or IMP.atom.Atom.get_is_setup(p)

class AsymIDMapper(object):
    """Map a Particle to an asym_id (chain ID)"""
    def __init__(self, prot):
        self.o = IMP.pmi.output.Output()
        self.prot = prot
        self.name = 'cif-output'
        self.o.dictionary_pdbs[self.name] = self.prot
        self.o._init_dictchain(self.name, self.prot)

    def __getitem__(self, p):
        protname, is_a_bead = self.o.get_prot_name_from_particle(self.name, p)
        return self.o.dictchain[self.name][protname]

class Dumper(object):
    """Base class for helpers to dump output to mmCIF"""
    def __init__(self, simo):
        self.simo = simo

    def finalize(self):
        pass


class EntryDumper(Dumper):
    def dump(self, writer):
        entry_id = 'imp_model'
        # Write CIF header (so this dumper should always be first)
        writer.fh.write("data_%s\n" % entry_id)
        with writer.category("_entry") as l:
            l.write(id=entry_id)


class SoftwareDumper(Dumper):
    software = [
       IMP.pmi.metadata.Software(
             name="Integrative Modeling Platform (IMP)",
             version=IMP.__version__,
             classification="integrative model building",
             description="integrative model building",
             url='https://integrativemodeling.org'),
       IMP.pmi.metadata.Software(
            name="IMP PMI module",
            version=IMP.pmi.__version__,
            classification="integrative model building",
            description="integrative model building",
            url='https://integrativemodeling.org')
       ]

    def __init__(self, simo):
        super(SoftwareDumper, self).__init__(simo)
        self.modeller_used = False

    def set_modeller_used(self, version, date):
        if self.modeller_used:
            return
        self.modeller_used = True
        self.software.append(IMP.pmi.metadata.Software(
                    name='MODELLER', classification='comparative modeling',
                    description='Comparative modeling by satisfaction '
                                'of spatial restraints, build ' + date,
                    url='https://salilab.org/modeller/',
                    version=version))

    def dump(self, writer):
        ordinal = 1
        with writer.loop("_software",
                         ["pdbx_ordinal", "name", "classification", "version",
                          "type", "location"]) as l:
            for m in self.software + self.simo._metadata:
                if isinstance(m, IMP.pmi.metadata.Software):
                    l.write(pdbx_ordinal=ordinal, name=m.name,
                            classification=m.classification, version=m.version,
                            type=m.type, location=m.url)
                    ordinal += 1

class AuditAuthorDumper(Dumper):
    """Populate the mmCIF audit_author category (a list of the people that
       authored this mmCIF file; here we assume that's just the authors of
       any associated publications)"""
    def dump(self, writer):
        citations = [m for m in self.simo._metadata
                     if isinstance(m, IMP.pmi.metadata.Citation)]
        seen_authors = {}
        with writer.loop("_audit_author",
                         ["name", "pdbx_ordinal"]) as l:
            ordinal = 1
            for n, c in enumerate(citations):
                for a in c.authors:
                    if a not in seen_authors:
                        seen_authors[a] = None
                        l.write(name=a, pdbx_ordinal=ordinal)
                        ordinal += 1


class CitationDumper(Dumper):
    def dump(self, writer):
        citations = [m for m in self.simo._metadata
                     if isinstance(m, IMP.pmi.metadata.Citation)]
        with writer.loop("_citation",
                         ["id", "title", "journal_abbrev", "journal_volume",
                          "page_first", "page_last", "year",
                          "pdbx_database_id_PubMed",
                          "pdbx_database_id_DOI"]) as l:
            for n, c in enumerate(citations):
                l.write(id=n+1, title=c.title, journal_abbrev=c.journal,
                        journal_volume=c.volume, page_first=c.page_range[0],
                        page_last=c.page_range[1], year=c.year,
                        pdbx_database_id_PubMed=c.pmid,
                        pdbx_database_id_DOI=c.doi)

        with writer.loop("_citation_author",
                         ["citation_id", "name", "ordinal"]) as l:
            ordinal = 1
            for n, c in enumerate(citations):
                for a in c.authors:
                    l.write(citation_id=n+1, name=a, ordinal=ordinal)
                    ordinal += 1


class EntityDumper(Dumper):
    # todo: we currently only support amino acid sequences here
    def dump(self, writer):
        all_entities = [x for x in sorted(self.simo.entities.items(),
                                          key=operator.itemgetter(1))]
        with writer.loop("_entity",
                         ["id", "type", "src_method", "pdbx_description",
                          "formula_weight", "pdbx_number_of_molecules",
                          "details"]) as l:
            for name, entity_id in all_entities:
                l.write(id=entity_id, type='polymer', src_method='man',
                        pdbx_description=name, formula_weight=writer.unknown,
                        pdbx_number_of_molecules=1, details=writer.unknown)


class EntityPolyDumper(Dumper):
    # todo: we currently only support amino acid sequences here
    def __init__(self, simo):
        super(EntityPolyDumper, self).__init__(simo)
        self.output = IMP.pmi.output.Output()

    def dump(self, writer):
        all_entities = [x for x in sorted(self.simo.entities.items(),
                                          key=operator.itemgetter(1))]
        with writer.loop("_entity_poly",
                         ["entity_id", "type", "nstd_linkage",
                          "nstd_monomer", "pdbx_strand_id",
                          "pdbx_seq_one_letter_code",
                          "pdbx_seq_one_letter_code_can"]) as l:
            for name, entity_id in all_entities:
                seq = self.simo.sequence_dict[name]
                first_copy = self.simo.copies[name][0]
                chain = self.simo.chains[(name, first_copy)]
                l.write(entity_id=entity_id, type='polypeptide(L)',
                        nstd_linkage='no', nstd_monomer='no',
                        pdbx_strand_id=self.output.chainids[chain],
                        pdbx_seq_one_letter_code=seq,
                        pdbx_seq_one_letter_code_can=seq)

class EntityPolySeqDumper(Dumper):
    def dump(self, writer):
        all_entities = [x for x in sorted(self.simo.entities.items(),
                                          key=operator.itemgetter(1))]
        with writer.loop("_entity_poly_seq",
                         ["entity_id", "num", "mon_id", "hetero"]) as l:
            for name, entity_id in all_entities:
                seq = self.simo.sequence_dict[name]
                for num, one_letter_code in enumerate(seq):
                    restyp = IMP.atom.get_residue_type(one_letter_code)
                    l.write(entity_id=entity_id, num=num + 1,
                            mon_id=restyp.get_string(),
                            hetero=CifWriter.omitted)

class StructAsymDumper(Dumper):
    def __init__(self, simo):
        super(StructAsymDumper, self).__init__(simo)
        self.output = IMP.pmi.output.Output()

    def dump(self, writer):
        with writer.loop("_struct_asym",
                         ["id", "entity_id", "details"]) as l:
            for comp in self.simo.all_components:
                entity_id = self.simo.entities[comp]
                for copy in self.simo.copies[comp]:
                    chain = self.simo.chains[(comp, copy)]
                    l.write(id=self.output.chainids[chain],
                            entity_id=entity_id,
                            details=copy)

class _PDBFragment(object):
    """Record details about part of a PDB file used as input
       for a component."""
    primitive = 'sphere'
    granularity = 'by-residue'
    num = CifWriter.omitted
    def __init__(self, m, component, start, end, offset, pdbname, chain, hier):
        self.component, self.start, self.end, self.offset, self.pdbname \
              = component, start, end, offset, pdbname
        self.chain, self.hier = chain, hier
        sel = IMP.atom.NonWaterNonHydrogenPDBSelector() \
              & IMP.atom.ChainPDBSelector(chain)
        self.starting_hier = IMP.atom.read_pdb(pdbname, m, sel)

    def combine(self, other):
        pass

class _BeadsFragment(object):
    """Record details about beads used to represent part of a component."""
    primitive = 'sphere'
    granularity = 'by-feature'
    chain = None
    def __init__(self, m, component, start, end, num, hier):
        self.component, self.start, self.end, self.num, self.hier \
              = component, start, end, num, hier

    def combine(self, other):
        # todo: don't combine if one fragment is rigid and the other flexible
        if type(other) == type(self) and other.start == self.end + 1:
            self.end = other.end
            self.num += other.num
            return True

class _StartingModel(object):
    """Record details about an input model (e.g. comparative modeling
       template) used for a component."""

    source = CifWriter.unknown
    db_name = CifWriter.unknown
    db_code = CifWriter.unknown
    sequence_identity = CifWriter.unknown

    def __init__(self, fragment):
        self.fragments = [fragment]

class ModelRepresentationDumper(Dumper):
    def __init__(self, simo):
        super(ModelRepresentationDumper, self).__init__(simo)
        # dict of fragments, ordered by component name
        self.fragments = OrderedDict()

    def add_fragment(self, fragment):
        """Add a model fragment."""
        comp = fragment.component
        if comp not in self.fragments:
            self.fragments[comp] = []
        fragments = self.fragments[comp]
        if len(fragments) == 0 or not fragments[-1].combine(fragment):
            fragments.append(fragment)

    def get_model_mode(self, fragment):
        """Determine the model_mode for a given fragment ('rigid' or
           'flexible')"""
        leaves = IMP.atom.get_leaves(fragment.hier)
        # Assume all leaves are set up as rigid/flexible in the same way
        if IMP.core.RigidMember.get_is_setup(leaves[0]):
            return 'rigid'
        else:
            return 'flexible'

    def dump(self, writer):
        segment_id = 1
        with writer.loop("_ihm_model_representation",
                         ["segment_id", "entity_id", "entity_description",
                          "seq_id_begin", "seq_id_end",
                          "model_object_primitive", "starting_model_id",
                          "model_mode", "model_granularity",
                          "model_object_count"]) as l:
            for comp, fragments in self.fragments.items():
                for f in fragments:
                    starting_model_id = CifWriter.omitted
                    if hasattr(f, 'pdbname'):
                        starting_model_id = self.starting_model[f.pdbname].name
                    l.write(segment_id=segment_id,
                            entity_id=self.simo.entities[f.component],
                            entity_description=f.component,
                            seq_id_begin=f.start,
                            seq_id_end=f.end,
                            model_object_primitive=f.primitive,
                            starting_model_id=starting_model_id,
                            model_mode=self.get_model_mode(f),
                            model_granularity=f.granularity,
                            model_object_count=f.num)
                    segment_id += 1

class PDBSource(object):
    """An experimental PDB file used as part of a starting model"""
    source = 'experimental model'
    db_name = 'PDB'
    sequence_identity = 100.0

    def __init__(self, model, db_code, chain_id, metadata):
        self.db_code = db_code
        self.chain_id = chain_id
        self.metadata = metadata

    def get_seq_id_range(self, model):
        # Assume the structure covers the entire sequence
        return (model.seq_id_begin, model.seq_id_end)

class TemplateSource(object):
    """A PDB file used as a template for a comparative starting model"""
    source = 'comparative model'
    db_name = 'PDB'

    def __init__(self, code, seq_id_begin, seq_id_end, seq_id,
                 model):
        # Assume a code of 1abcX refers to a real PDB structure
        if len(code) == 5:
            self.db_code = code[:4].upper()
            self.chain_id = code[4]
        else:
            self.db_code = self.chain_id = CifWriter.unknown
        self.sequence_identity = seq_id
        self._seq_id_begin, self._seq_id_end = seq_id_begin, seq_id_end

    def get_seq_id_range(self, model):
        # The template may cover more than the current starting model
        seq_id_begin = max(model.seq_id_begin, self._seq_id_begin)
        seq_id_end = min(model.seq_id_end, self._seq_id_end)
        return (seq_id_begin, seq_id_end)

class UnknownSource(object):
    """Part of a starting model from an unknown source"""
    source = CifWriter.unknown
    db_code = CifWriter.unknown
    db_name = CifWriter.unknown
    chain_id = CifWriter.unknown
    sequence_identity = CifWriter.unknown

    def __init__(self, model):
        pass

    def get_seq_id_range(self, model):
        return (model.seq_id_begin, model.seq_id_end)

class DatasetLocation(object):
    """External location of a dataset"""
    pass

class RepoDatasetLocation(DatasetLocation):
    """Pointer to a dataset stored in a repository"""
    doi = content_filename = CifWriter.unknown

    def __init__(self, loc):
        if loc:
            self.doi = loc.doi
            self.content_filename = loc.path

class DBDatasetLocation(DatasetLocation):
    """Pointer to a dataset stored in an official database (e.g. PDB)"""
    def __init__(self, db_name, db_code, version, details):
        self.db_name = db_name
        self.access_code = db_code
        self.version, self.details = version, details

class Dataset(object):
    location = None

    def set_location(self, loc):
        self.location = RepoDatasetLocation(loc)

class CXMSDataset(Dataset):
    data_type = 'CX-MS data'

class EMMicrographsDataset(Dataset):
    data_type = 'EM raw micrographs'
    def __init__(self, number):
        self.number = number

class EM2DClassDataset(Dataset):
    data_type = '2DEM class average'

class CompModelDataset(Dataset):
    """A comparative model dataset.
       Currently it is assumed that models are stored in a repository."""
    data_type = 'Comparative model'
    def __init__(self, location):
        self.set_location(location)

class PDBDataset(Dataset):
    """An experimental PDB structure dataset."""
    data_type = 'experimental model'
    def __init__(self, db_code, version, details):
        self.location = DBDatasetLocation('PDB', db_code, version, details)

class DatasetGroup(object):
    """A group of datasets"""
    def __init__(self, datasets):
        self.datasets = datasets[:]

class DatasetDumper(Dumper):
    def __init__(self, simo):
        super(DatasetDumper, self).__init__(simo)
        self.datasets = []
        self.dataset_groups = {}

    def get_all_group(self):
        """Get a DatasetGroup encompassing all datasets so far"""
        num_datasets = len(self.datasets)
        if num_datasets not in self.dataset_groups:
            g = DatasetGroup(self.datasets)
            self.dataset_groups[num_datasets] = g
            g.id = len(self.dataset_groups)
        return self.dataset_groups[num_datasets]

    def add(self, dataset):
        self.datasets.append(dataset)
        dataset.id = len(self.datasets)

    def dump(self, writer):
        ordinal = 1
        groups = sorted(self.dataset_groups.values(), key=lambda x: x.id)
        with writer.loop("_ihm_dataset_list",
                         ["ordinal_id", "id", "group_id", "data_type",
                          "database_hosted"]) as l:
            for g in groups:
                for d in g.datasets:
                    l.write(ordinal_id=ordinal, id=d.id, group_id=g.id,
                            data_type=d.data_type,
                            database_hosted=not isinstance(d.location,
                                                           RepoDatasetLocation))
                    ordinal += 1
        self.dump_other((d for d in self.datasets
                         if isinstance(d.location, RepoDatasetLocation)),
                        writer)
        self.dump_rel_dbs((d for d in self.datasets
                           if isinstance(d.location, DBDatasetLocation)),
                          writer)

    def dump_rel_dbs(self, datasets, writer):
        ordinal = 1
        with writer.loop("_ihm_dataset_related_db_reference",
                         ["id", "dataset_list_id", "db_name",
                          "access_code", "version", "data_type",
                          "details"]) as l:
            for d in datasets:
                l.write(id=ordinal, dataset_list_id=d.id,
                        db_name=d.location.db_name,
                        access_code=d.location.access_code,
                        version=d.location.version,
                        data_type=d.data_type, details=d.location.details)
                ordinal += 1

    def dump_other(self, datasets, writer):
        ordinal = 1
        with writer.loop("_ihm_dataset_other",
                         ["id", "dataset_list_id", "data_type",
                          "doi", "content_filename"]) as l:
            for d in datasets:
                l.write(id=ordinal, dataset_list_id=d.id,
                        data_type=d.data_type, doi=d.location.doi,
                        content_filename=d.location.content_filename)
                ordinal += 1


class ExperimentalCrossLink(object):
    def __init__(self, r1, c1, r2, c2, label, length, dataset):
        self.r1, self.c1, self.r2, self.c2, self.label = r1, c1, r2, c2, label
        self.length, self.dataset = length, dataset

class CrossLink(object):
    def __init__(self, ex_xl, p1, p2, sigma1, sigma2, psi):
        self.ex_xl, self.sigma1, self.sigma2 = ex_xl, sigma1, sigma2
        self.p1, self.p2 = p1, p2
        self.psi = psi

class CrossLinkDumper(Dumper):
    def __init__(self, simo):
        super(CrossLinkDumper, self).__init__(simo)
        self.cross_links = []
        self.exp_cross_links = []

    def add_experimental(self, cross_link):
        self.exp_cross_links.append(cross_link)
        cross_link.id = len(self.exp_cross_links)

    def add(self, cross_link):
        self.cross_links.append(cross_link)
        cross_link.id = len(self.cross_links)

    def dump(self, writer):
        self.dump_list(writer)
        self.dump_restraint(writer)

    def dump_list(self, writer):
        with writer.loop("_ihm_cross_link_list",
                         ["id", "group_id", "entity_description_1",
                          "entity_id_1", "seq_id_1", "comp_id_1",
                          "entity_description_2",
                          "entity_id_2", "seq_id_2", "comp_id_2", "type",
                          "dataset_list_id"]) as l:
            for xl in self.exp_cross_links:
                seq1 = self.simo.sequence_dict[xl.c1]
                seq2 = self.simo.sequence_dict[xl.c2]
                rt1 = IMP.atom.get_residue_type(seq1[xl.r1-1])
                rt2 = IMP.atom.get_residue_type(seq2[xl.r2-1])
                l.write(id=xl.id,
                        entity_description_1=xl.c1,
                        entity_id_1=self.simo.entities[xl.c1],
                        seq_id_1=xl.r1,
                        comp_id_1=rt1.get_string(),
                        entity_description_2=xl.c2,
                        entity_id_2=self.simo.entities[xl.c2],
                        seq_id_2=xl.r2,
                        comp_id_2=rt2.get_string(),
                        type=xl.label,
                        dataset_list_id=xl.dataset.id)

    def _granularity(self, xl):
        """Determine the granularity of a cross link"""
        if _get_by_residue(xl.p1) and _get_by_residue(xl.p2):
            return 'by-residue'
        else:
            return 'by-feature'

    def dump_restraint(self, writer):
        asym = AsymIDMapper(self.simo.prot)
        with writer.loop("_ihm_cross_link_restraint",
                         ["id", "group_id", "entity_id_1", "asym_id_1",
                          "seq_id_1", "comp_id_1",
                          "entity_id_2", "asym_id_2", "seq_id_2", "comp_id_2",
                          "type", "conditional_crosslink_flag",
                          "model_granularity", "distance_threshold",
                          "psi", "sigma_1", "sigma_2"]) as l:
            for xl in self.cross_links:
                seq1 = self.simo.sequence_dict[xl.ex_xl.c1]
                seq2 = self.simo.sequence_dict[xl.ex_xl.c2]
                rt1 = IMP.atom.get_residue_type(seq1[xl.ex_xl.r1-1])
                rt2 = IMP.atom.get_residue_type(seq2[xl.ex_xl.r2-1])
                l.write(id=xl.id,
                        group_id=xl.ex_xl.id,
                        entity_id_1=self.simo.entities[xl.ex_xl.c1],
                        asym_id_1=asym[xl.p1],
                        seq_id_1=xl.ex_xl.r1,
                        comp_id_1=rt1.get_string(),
                        entity_id_2=self.simo.entities[xl.ex_xl.c2],
                        asym_id_2=asym[xl.p2],
                        seq_id_2=xl.ex_xl.r2,
                        comp_id_2=rt2.get_string(),
                        type=xl.ex_xl.label,
                        # todo: any circumstances where this could be ANY?
                        conditional_crosslink_flag="ALL",
                        model_granularity=self._granularity(xl),
                        distance_threshold=xl.ex_xl.length,
                        psi=xl.psi, sigma_1=xl.sigma1, sigma_2=xl.sigma2)

class EM2DRestraint(object):
    def __init__(self, dataset, resolution, pixel_size,
                 image_resolution, projection_number, micrographs_dataset):
        self.dataset, self.resolution = dataset, resolution
        self.pixel_size, self.image_resolution = pixel_size, image_resolution
        self.projection_number = projection_number
        self.micrographs_dataset = micrographs_dataset

class EM2DDumper(Dumper):
    def __init__(self, simo):
        super(EM2DDumper, self).__init__(simo)
        self.restraints = []

    def add(self, rsr):
        self.restraints.append(rsr)
        rsr.id = len(self.restraints)

    def dump(self, writer):
        with writer.loop("_ihm_2dem_class_average_restraint",
                         ["id", "dataset_list_id", "number_raw_micrographs",
                          "raw_micrographs_dataset_list_id",
                          "pixel_size_width", "pixel_size_height",
                          "image_resolution", "image_segment_flag",
                          "number_of_projections", "struct_assembly_id",
                          "details"]) as l:
            for r in self.restraints:
                mgd = r.micrographs_dataset
                unk = CifWriter.unknown
                l.write(id=r.id, dataset_list_id=r.dataset.id,
                        number_raw_micrographs=mgd.number if mgd else unk,
                        raw_micrographs_dataset_list_id=mgd.id if mgd else unk,
                        pixel_size_width=r.pixel_size,
                        pixel_size_height=r.pixel_size,
                        image_resolution=r.image_resolution,
                        number_of_projections=r.projection_number,
                        struct_assembly_id=self.simo.default_assembly.id,
                        image_segment_flag=False)

class Assembly(list):
    """A collection of components. Currently simply implemented as a list of
       the component names."""
    pass

class AssemblyDumper(Dumper):
    def __init__(self, simo):
        super(AssemblyDumper, self).__init__(simo)
        self.assemblies = []
        self.output = IMP.pmi.output.Output()

    def add(self, a):
        self.assemblies.append(a)
        a.id = len(self.assemblies)

    def dump(self, writer):
        ordinal = 1
        with writer.loop("_ihm_struct_assembly",
                         ["ordinal_id", "assembly_id", "entity_description",
                          "entity_id", "asym_id", "seq_id_begin",
                          "seq_id_end"]) as l:
            for a in self.assemblies:
                for comp in a:
                    for copy in self.simo.copies[comp]:
                        seq = self.simo.sequence_dict[comp]
                        chain = self.simo.chains[(comp, copy)]
                        l.write(ordinal_id=ordinal, assembly_id=a.id,
                                entity_description=comp,
                                entity_id=self.simo.entities[comp],
                                asym_id=self.output.chainids[chain],
                                seq_id_begin=1,
                                seq_id_end=len(seq))
                        ordinal += 1


class ReplicaExchangeProtocol(object):
    def __init__(self, rex):
        if rex.monte_carlo_sample_objects is not None:
            self.step_method = 'Replica exchange monte carlo'
        else:
            self.step_method = 'Replica exchange molecular dynamics'
        self.num_models_end = rex.vars["number_of_frames"]

class Model(object):
    def __init__(self, prot, simo):
        o = IMP.pmi.output.Output()
        name = 'cif-output'
        o.dictionary_pdbs[name] = prot
        o._init_dictchain(name, prot)
        (particle_infos_for_pdb,
         self.geometric_center) = o.get_particle_infos_for_pdb_writing(name)
        self.entity_for_chain = {}
        for protname, chain_id in o.dictchain[name].items():
            self.entity_for_chain[chain_id] = simo.entities[protname]
        self.spheres = [t for t in particle_infos_for_pdb if t[1] is None]
        self.atoms = [t for t in particle_infos_for_pdb if t[1] is not None]

class ModelDumper(Dumper):
    def __init__(self, simo):
        super(ModelDumper, self).__init__(simo)
        self.models = []

    def add(self, prot):
        m = Model(prot, self.simo)
        self.models.append(m)
        m.id = len(self.models)
        return m.id

    def dump(self, writer):
        num_atoms = sum(len(m.atoms) for m in self.models)
        num_spheres = sum(len(m.spheres) for m in self.models)
        self.dump_atoms(writer)
        self.dump_spheres(writer)

    def dump_atoms(self, writer):
        ordinal = 1
        with writer.loop("_atom_site",
                         ["id", "label_atom_id", "label_comp_id",
                          "label_seq_id",
                          "label_asym_id", "Cartn_x",
                          "Cartn_y", "Cartn_z", "label_entity_id",
                          "model_id"]) as l:
            for model in self.models:
                for atom in model.atoms:
                    (xyz, atom_type, residue_type, chain_id, residue_index,
                     all_indexes, radius) = atom
                    l.write(id=ordinal, label_atom_id=atom_type.get_string(),
                            label_comp_id=residue_type.get_string(),
                            label_asym_id=chain_id,
                            label_entity_id=model.entity_for_chain[chain_id],
                            label_seq_id=residue_index,
                            Cartn_x=xyz[0] - model.geometric_center[0],
                            Cartn_y=xyz[1] - model.geometric_center[1],
                            Cartn_z=xyz[2] - model.geometric_center[2],
                            model_id=model.id)
                    ordinal += 1

    def dump_spheres(self, writer):
        ordinal = 1
        with writer.loop("_ihm_sphere_obj_site",
                         ["ordinal_id", "entity_id", "seq_id_begin",
                          "seq_id_end", "asym_id", "Cartn_x",
                          "Cartn_y", "Cartn_z", "object_radius",
                          "model_id"]) as l:
            for model in self.models:
                for sphere in model.spheres:
                    (xyz, atom_type, residue_type, chain_id, residue_index,
                     all_indexes, radius) = sphere
                    if all_indexes is None:
                        all_indexes = (residue_index,)
                    l.write(ordinal_id=ordinal,
                            entity_id=model.entity_for_chain[chain_id],
                            seq_id_begin = all_indexes[0],
                            seq_id_end = all_indexes[-1],
                            asym_id=chain_id,
                            Cartn_x=xyz[0] - model.geometric_center[0],
                            Cartn_y=xyz[1] - model.geometric_center[1],
                            Cartn_z=xyz[2] - model.geometric_center[2],
                            object_radius=radius, model_id=model.id)
                    ordinal += 1


class ModelProtocolDumper(Dumper):
    def __init__(self, simo):
        super(ModelProtocolDumper, self).__init__(simo)
        self.protocols = []

    def add(self, protocol):
        self.protocols.append(protocol)
        protocol.id = len(self.protocols)
        # Assume that protocol uses all currently-defined datasets
        protocol.dataset_group = self.simo.dataset_dump.get_all_group()

    def dump(self, writer):
        ordinal = 1
        with writer.loop("_ihm_modeling_protocol",
                         ["ordinal_id", "protocol_id", "step_id",
                          "struct_assembly_id", "dataset_group_id",
                          "assembly_description", "protocol_name",
                          "step_name", "step_method", "num_models_begin",
                          "num_models_end", "multi_scale_flag",
                          "multi_state_flag", "time_ordered_flag"]) as l:
            # todo: handle multiple protocols (e.g. sampling then refinement)
            num_models_begin = 0
            for p in self.protocols:
                l.write(ordinal_id=ordinal, protocol_id=1,
                        step_id=p.id, step_method=p.step_method,
                        step_name='Sampling',
                        struct_assembly_id=self.simo.default_assembly.id,
                        dataset_group_id=p.dataset_group.id,
                        num_models_begin=num_models_begin,
                        num_models_end=p.num_models_end)
                num_models_begin = p.num_models_end
                ordinal += 1


class PDBHelix(object):
    """Represent a HELIX record from a PDB file."""
    def __init__(self, line):
        self.helix_id = line[11:14].strip()
        self.start_resnam = line[14:18].strip()
        self.start_asym = line[19]
        self.start_resnum = int(line[21:25])
        self.end_resnam = line[27:30].strip()
        self.end_asym = line[31]
        self.end_resnum = int(line[33:37])
        self.helix_class = int(line[38:40])
        self.length = int(line[71:76])


class StartingModelDumper(Dumper):
    def __init__(self, simo):
        super(StartingModelDumper, self).__init__(simo)
        # dict of PDB fragments, ordered by component name
        self.fragments = OrderedDict()
        # dict of starting models (entire PDB files), collected from fragments
        self.models = OrderedDict()
        # mapping from pdbname to starting model
        self.starting_model = {}

    def add_pdb_fragment(self, fragment):
        """Add a starting model PDB fragment."""
        comp = fragment.component
        if comp not in self.fragments:
            self.fragments[comp] = []
            self.models[comp] = []
        self.fragments[comp].append(fragment)
        models = self.models[comp]
        if len(models) == 0 \
           or models[-1].fragments[0].pdbname != fragment.pdbname:
            model = _StartingModel(fragment)
            models.append(model)
            model.sources = self.get_sources(model, fragment.pdbname,
                                             fragment.chain)
        else:
            models[-1].fragments.append(fragment)

    def get_templates(self, pdbname, model):
        templates = []
        tmpre = re.compile('REMARK   6 TEMPLATE: (\S+) .* '
                           'MODELS (\S+):\S+ \- (\S+):\S+ AT (\S+)%')

        with open(pdbname) as fh:
            for line in fh:
                if line.startswith('ATOM'): # Read only the header
                    break
                m = tmpre.match(line)
                if m:
                    templates.append(TemplateSource(m.group(1),
                                                    int(m.group(2)),
                                                    int(m.group(3)),
                                                    m.group(4), model))
        # Sort by starting residue, then ending residue
        return sorted(templates, key=lambda x: (x._seq_id_begin, x._seq_id_end))

    def _parse_pdb(self, fh, first_line):
        """Extract information from an official PDB"""
        metadata = []
        details = ''
        for line in fh:
            if line.startswith('TITLE'):
                details += line[10:].rstrip()
            elif line.startswith('HELIX'):
                metadata.append(PDBHelix(line))
        return (first_line[50:59].strip(),
                details if details else CifWriter.unknown, metadata)

    def get_sources(self, model, pdbname, chain):
        # Attempt to identity PDB file vs. comparative model
        fh = open(pdbname)
        first_line = fh.readline()
        if first_line.startswith('HEADER'):
            version, details, metadata = self._parse_pdb(fh, first_line)
            source = PDBSource(model, first_line[62:66].strip(), chain,
                               metadata)
            model.dataset = PDBDataset(source.db_code, version, details)
            self.simo.dataset_dump.add(model.dataset)
            return [source]
        elif first_line.startswith('EXPDTA    THEORETICAL MODEL, MODELLER'):
            self.simo.software_dump.set_modeller_used(
                                        *first_line[38:].split(' ', 1))
            model.dataset = CompModelDataset(self.simo._get_location(pdbname))
            self.simo.dataset_dump.add(model.dataset)
            templates = self.get_templates(pdbname, model)
            if templates:
                return templates
        return [UnknownSource(model)]

    def assign_model_details(self):
        for comp, models in self.models.items():
            for i, model in enumerate(models):
                model.name = "%s-m%d" % (comp, i+1)
                self.starting_model[model.fragments[0].pdbname] = model
                model.seq_id_begin = min(x.start + x.offset
                                         for x in model.fragments)
                model.seq_id_end = max(x.end + x.offset
                                       for x in model.fragments)

    def all_models(self):
        for comp, models in self.models.items():
            for model in models:
                yield model

    def finalize(self):
        self.assign_model_details()

    def dump(self, writer):
        self.dump_details(writer)
        self.dump_coords(writer)

    def dump_details(self, writer):
        writer.write_comment("""IMP will attempt to identify which input models
are crystal structures and which are comparative models, but does not always
have sufficient information to deduce all of the templates used for comparative
modeling. These may need to be added manually below.""")
        with writer.loop("_ihm_starting_model_details",
                     ["id", "entity_id", "entity_description", "seq_id_begin",
                      "seq_id_end", "starting_model_source",
                      "starting_model_db_name", "starting_model_db_code",
                      "starting_model_db_pdb_auth_seq_id",
                      "starting_model_sequence_identity",
                      "starting_model_id",
                      "dataset_list_id"]) as l:
            ordinal = 1
            for model in self.all_models():
                f = model.fragments[0]
                for source in model.sources:
                    seq_id_begin, seq_id_end = source.get_seq_id_range(model)
                    l.write(id=ordinal,
                      entity_id=self.simo.entities[f.component],
                      entity_description=f.component,
                      seq_id_begin=seq_id_begin,
                      seq_id_end=seq_id_end,
                      starting_model_db_pdb_auth_seq_id=source.chain_id,
                      starting_model_id=model.name,
                      starting_model_source=source.source,
                      starting_model_db_name=source.db_name,
                      starting_model_db_code=source.db_code,
                      starting_model_sequence_identity=source.sequence_identity,
                      dataset_list_id=model.dataset.id)
                    ordinal += 1

    def dump_coords(self, writer):
        ordinal = 1
        with writer.loop("_ihm_starting_model_coord",
                     ["starting_model_id", "group_PDB", "id", "type_symbol",
                      "atom_id", "comp_id", "entity_id", "seq_id", "Cartn_x",
                      "Cartn_y", "Cartn_z", "B_iso_or_equiv",
                      "ordinal_id"]) as l:
            for model in self.all_models():
                for f in model.fragments:
                    for a in IMP.atom.get_leaves(f.starting_hier):
                        coord = IMP.core.XYZ(a).get_coordinates()
                        atom = IMP.atom.Atom(a)
                        element = atom.get_element()
                        element = IMP.atom.get_element_table().get_name(element)
                        atom_name = atom.get_atom_type().get_string()
                        group_pdb = 'ATOM'
                        if atom_name.startswith('HET:'):
                            group_pdb = 'HETATM'
                            atom_name = atom_name[4:]
                        res = IMP.atom.get_residue(atom)
                        res_name = res.get_residue_type().get_string()
                        chain = IMP.atom.get_chain(res)
                        l.write(starting_model_id=model.name,
                                group_PDB=group_pdb,
                                id=atom.get_input_index(), type_symbol=element,
                                atom_id=atom_name, comp_id=res_name,
                                entity_id=self.simo.entities[f.component],
                                seq_id=res.get_index(), Cartn_x=coord[0],
                                Cartn_y=coord[1], Cartn_z=coord[2],
                                B_iso_or_equiv=atom.get_temperature_factor(),
                                ordinal_id=ordinal)
                        ordinal += 1

class StructConfDumper(Dumper):
    def all_rigid_fragments(self):
        """Yield all rigid model representation fragments"""
        asym = AsymIDMapper(self.simo.prot)
        model_repr = self.simo.model_repr_dump
        for comp, fragments in model_repr.fragments.items():
            for f in fragments:
                if hasattr(f, 'pdbname') \
                   and model_repr.get_model_mode(f) == 'rigid':
                    yield (f, model_repr.starting_model[f.pdbname],
                           asym[f.hier])

    def all_helices(self):
        """Yield all helices that overlap with rigid model fragments"""
        for f, starting_model, asym_id in self.all_rigid_fragments():
            if len(starting_model.sources) == 1 \
               and isinstance(starting_model.sources[0], PDBSource):
                pdb = starting_model.sources[0]
                for m in pdb.metadata:
                    if isinstance(m, PDBHelix) \
                       and m.start_asym == pdb.chain_id \
                       and m.end_asym == pdb.chain_id \
                       and m.start_resnum >= f.start and m.end_resnum <= f.end:
                        yield (m, max(f.start, m.start_resnum),
                               min(f.end, m.end_resnum), asym_id)

    def dump(self, writer):
        with writer.category("_struct_conf_type") as l:
            l.write(id='HELX_P', criteria=CifWriter.unknown,
                    reference=CifWriter.unknown)
        # Dump helix information for the model. For any model fragment that
        # is rigid and uses an experimental PDB structure as the starting
        # model, inherit any helix information from that PDB file.
        # Note that we can't use the helix id from the original PDB, since
        # it has to be unique and this might not be the case if we inherit
        # from multiple PDBs.
        ordinal = 0
        with writer.loop("_struct_conf",
                     ["id", "conf_type_id", "beg_label_comp_id",
                      "beg_label_asym_id", "beg_label_seq_id",
                      "end_label_comp_id", "end_label_asym_id",
                      "end_label_seq_id",]) as l:
            for h, begin, end, asym_id in self.all_helices():
                ordinal += 1
                l.write(id='HELX_P%d' % ordinal, conf_type_id='HELX_P',
                        beg_label_comp_id=h.start_resnam,
                        beg_label_asym_id=asym_id,
                        beg_label_seq_id=begin,
                        end_label_comp_id=h.end_resnam,
                        end_label_asym_id=asym_id,
                        end_label_seq_id=end)


class CifEntities(dict):
    """Handle mapping from IMP components to CIF entity IDs.
       An entity is a chain with a unique sequence. Thus, multiple
       components may map to the same entity if they share sequence."""
    def __init__(self):
        super(CifEntities, self).__init__()
        self._sequence_dict = {}
        self._first_component = {}

    def add(self, component_name, sequence):
        if sequence not in self._sequence_dict:
            entity_id = len(self._sequence_dict) + 1
            self._first_component[entity_id] = component_name
            self._sequence_dict[sequence] = entity_id
        self[component_name] = self._sequence_dict[sequence]

    def get_all(self):
        """Yield all (entity, component) pairs"""
        return self._first_component.items()


class ProtocolOutput(IMP.pmi.output.ProtocolOutput):
    def __init__(self, fh):
        self._cif_writer = CifWriter(fh)
        self.entities = CifEntities()
        self.chains = {}
        self.copies = {}
        self._default_copy = {}
        self.sequence_dict = {}
        self.all_components = []
        self.model_repr_dump = ModelRepresentationDumper(self)
        self.cross_link_dump = CrossLinkDumper(self)
        self.em2d_dump = EM2DDumper(self)
        self.model_prot_dump = ModelProtocolDumper(self)
        self.dataset_dump = DatasetDumper(self)
        self.starting_model_dump = StartingModelDumper(self)
        self.assembly_dump = AssemblyDumper(self)
        self.default_assembly = Assembly()
        self.assembly_dump.add(self.default_assembly)
        self.model_dump = ModelDumper(self)
        self.model_repr_dump.starting_model \
                    = self.starting_model_dump.starting_model
        self.software_dump = SoftwareDumper(self)
        self._dumpers = [EntryDumper(self), # should always be first
                         AuditAuthorDumper(self),
                         self.software_dump, CitationDumper(self),
                         EntityDumper(self),
                         EntityPolyDumper(self), EntityPolySeqDumper(self),
                         StructAsymDumper(self),
                         self.assembly_dump,
                         self.model_repr_dump, self.dataset_dump,
                         self.cross_link_dump,
                         self.em2d_dump,
                         self.starting_model_dump,
                         StructConfDumper(self),
                         self.model_prot_dump, self.model_dump]

    def create_component(self, name):
        self.all_components.append(name)
        self.default_assembly.append(name)
        # Set the component to have only a single copy to start with.
        self.copies[name] = [name]
        self._default_copy[name] = True
        self.chains[(name, name)] = len(self.chains)

    def add_component_sequence(self, name, seq):
        self.sequence_dict[name] = seq
        self.entities.add(name, seq)

    def add_copy(self, component_name, copy_name):
        # Remove any default-constructed copy
        if self._default_copy[component_name]:
            self.copies[component_name] = []
            self._default_copy[component_name] = False
            del self.chains[(component_name, component_name)]
        self.copies[component_name].append(copy_name)
        self.chains[(component_name, copy_name)] = len(self.chains)

    def flush(self):
        for dumper in self._dumpers:
            dumper.finalize()
        for dumper in self._dumpers:
            dumper.dump(self._cif_writer)

    def add_pdb_element(self, name, start, end, offset, pdbname, chain, hier):
        p = _PDBFragment(self.m, name, start, end, offset, pdbname, chain, hier)
        self.model_repr_dump.add_fragment(p)
        self.starting_model_dump.add_pdb_fragment(p)

    def add_bead_element(self, name, start, end, num, hier):
        b = _BeadsFragment(self.m, name, start, end, num, hier)
        self.model_repr_dump.add_fragment(b)

    def get_cross_link_dataset(self, fname):
        d = CXMSDataset()
        d.set_location(self._get_location(fname))
        self.dataset_dump.add(d)
        return d

    def add_experimental_cross_link(self, r1, c1, r2, c2, label, length,
                                    dataset):
        xl = ExperimentalCrossLink(r1, c1, r2, c2, label, length, dataset)
        self.cross_link_dump.add_experimental(xl)
        return xl

    def add_cross_link(self, ex_xl, p1, p2, sigma1, sigma2, psi):
        self.cross_link_dump.add(CrossLink(ex_xl, p1, p2, sigma1, sigma2, psi))

    def add_replica_exchange(self, rex):
        self.model_prot_dump.add(ReplicaExchangeProtocol(rex))

    def add_em2d_restraint(self, images, resolution, pixel_size,
                           image_resolution, projection_number, micrographs):
        mgd = None
        if micrographs:
            mgd = EMMicrographsDataset(micrographs.number)
            mgd.set_location(self._get_location(None, micrographs.metadata))
            self.dataset_dump.add(mgd)
        for image in images:
            d = EM2DClassDataset()
            d.set_location(self._get_location(image))
            self.dataset_dump.add(d)
            self.em2d_dump.add(EM2DRestraint(d, resolution, pixel_size,
                                      image_resolution, projection_number, mgd))

    def add_model(self):
        return self.model_dump.add(self.prot)

    def _get_location(self, path, metadata=[]):
        """Get the location where the given file is deposited, or None.
           If a RepositoryFile object is available, return that.
           Otherwise, if a Repository object is available, construct a
           RepositoryFile object from that and return it."""
        for m in metadata + self._metadata:
            if isinstance(m, IMP.pmi.metadata.RepositoryFile):
                return m
            if isinstance(m, IMP.pmi.metadata.Repository):
                return m.get_path(path)
