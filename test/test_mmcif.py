from __future__ import print_function
import IMP.test
import IMP.pmi.metadata
import IMP.pmi.representation
import IMP.pmi.mmcif
import sys
import io
if sys.version_info[0] >= 3:
    from io import StringIO
else:
    from io import BytesIO as StringIO

class EmptyObject(object):
    pass

class Tests(IMP.test.TestCase):

    def test_software(self):
        """Test SoftwareDumper"""
        s = IMP.pmi.metadata.Software(name='test', classification='test code',
                                      description='Some test program',
                                      version=1, url='http://salilab.org')
        m = IMP.Model()
        r = IMP.pmi.representation.Representation(m)
        r.add_metadata(s)
        d = IMP.pmi.mmcif.SoftwareDumper(r)
        fh = StringIO()
        w = IMP.pmi.mmcif.CifWriter(fh)
        d.dump(w)
        out = fh.getvalue().split('\n')
        self.assertEqual(out[-3],
                         "3 test 'test code' 1 program http://salilab.org")

    def test_assembly(self):
        """Test AssemblyDumper"""
        class DummyPO(IMP.pmi.mmcif.ProtocolOutput):
            def flush(self):
                pass
        po = DummyPO(EmptyObject())
        d = IMP.pmi.mmcif.AssemblyDumper(po)
        for c, seq in (("foo", "AAA"), ("bar", "AAA"), ("baz", "AA")):
            po.create_component(c)
            po.add_component_sequence(c, seq)
        d.add(IMP.pmi.mmcif.Assembly(["foo", "bar"]))
        d.add(IMP.pmi.mmcif.Assembly(["bar", "baz"]))

        fh = StringIO()
        w = IMP.pmi.mmcif.CifWriter(fh)
        d.dump(w)
        out = fh.getvalue()
        self.assertEqual(out, """#
loop_
_ihm_struct_assembly.ordinal_id
_ihm_struct_assembly.assembly_id
_ihm_struct_assembly.entity_description
_ihm_struct_assembly.entity_id
_ihm_struct_assembly.asym_id
_ihm_struct_assembly.seq_id_begin
_ihm_struct_assembly.seq_id_end
1 1 foo 1 A 1 3
2 1 foo 1 B 1 3
3 2 foo 1 B 1 3
4 2 baz 2 C 1 2
#
""")

    def test_struct_asym(self):
        """Test StructAsymDumper"""
        class DummyPO(IMP.pmi.mmcif.ProtocolOutput):
            def flush(self):
                pass
        po = DummyPO(EmptyObject())
        d = IMP.pmi.mmcif.StructAsymDumper(po)
        for c, seq in (("foo", "AAA"), ("bar", "AAA"), ("baz", "AA")):
            po.create_component(c)
            po.add_component_sequence(c, seq)

        fh = StringIO()
        w = IMP.pmi.mmcif.CifWriter(fh)
        d.dump(w)
        out = fh.getvalue()
        self.assertEqual(out, """#
loop_
_struct_asym.id
_struct_asym.entity_id
_struct_asym.details
A 1 foo
B 1 bar
C 2 baz
#
""")

    def test_citation(self):
        """Test CitationDumper"""
        s = IMP.pmi.metadata.Citation(
              pmid='25161197',
              title="Structural characterization by cross-linking reveals the "
                    "detailed architecture of a coatomer-related heptameric "
                    "module from the nuclear pore complex.",
              journal="Mol Cell Proteomics", volume=13, page_range=(2927,2943),
              year=2014,
              authors=['Shi Y', 'Fernandez-Martinez J', 'Tjioe E', 'Pellarin R',
                       'Kim SJ', 'Williams R', 'Schneidman-Duhovny D', 'Sali A',
                       'Rout MP', 'Chait BT'],
              doi='10.1074/mcp.M114.041673')

        m = IMP.Model()
        r = IMP.pmi.representation.Representation(m)
        r.add_metadata(s)
        d = IMP.pmi.mmcif.CitationDumper(r)
        fh = StringIO()
        w = IMP.pmi.mmcif.CifWriter(fh)
        d.dump(w)
        out = fh.getvalue()
        expected = """#
loop_
_citation.id
_citation.title
_citation.journal_abbrev
_citation.journal_volume
_citation.page_first
_citation.page_last
_citation.year
_citation.pdbx_database_id_PubMed
_citation.pdbx_database_id_DOI
1
;Structural characterization by cross-linking reveals the detailed arch
itecture of a coatomer-related heptameric module from the nuclear pore
 complex.
;
'Mol Cell Proteomics' 13 2927 2943 2014 25161197 10.1074/mcp.M114.041673
#
#
loop_
_citation_author.citation_id
_citation_author.name
_citation_author.ordinal
1 'Shi Y' 1
1 'Fernandez-Martinez J' 2
1 'Tjioe E' 3
1 'Pellarin R' 4
1 'Kim SJ' 5
1 'Williams R' 6
1 'Schneidman-Duhovny D' 7
1 'Sali A' 8
1 'Rout MP' 9
1 'Chait BT' 10
#
"""
        self.assertEqual(out, expected)

    def test_pdb_helix(self):
        """Test PDBHelix class"""
        p = IMP.pmi.mmcif.PDBHelix("HELIX   10  10 ASP A  607  GLU A  624  1                                  18   ")
        self.assertEqual(p.helix_id, '10')
        self.assertEqual(p.start_asym, 'A')
        self.assertEqual(p.start_resnum, 607)
        self.assertEqual(p.end_asym, 'A')
        self.assertEqual(p.end_resnum, 624)
        self.assertEqual(p.helix_class, 1)
        self.assertEqual(p.length, 18)

    def test_asym_id_mapper(self):
        """Test AsymIDMapper class"""
        m = IMP.Model()
        simo = IMP.pmi.representation.Representation(m)
        simo.create_component("Nup84")
        simo.add_component_sequence("Nup84",
                                    self.get_input_file_name("test.fasta"))
        simo.create_component("Nup85")
        simo.add_component_sequence("Nup85",
                                    self.get_input_file_name("test.fasta"))
        h1 = simo.add_component_beads("Nup84", [(1,2), (3,4)])
        h2 = simo.add_component_beads("Nup85", [(1,2), (3,4)])
        mapper = IMP.pmi.mmcif.AsymIDMapper(simo.prot)
        self.assertEqual(mapper[h1[0]], 'A')
        self.assertEqual(mapper[h1[1]], 'A')
        self.assertEqual(mapper[h2[0]], 'B')
        self.assertEqual(mapper[h2[1]], 'B')

    def test_cif_entities(self):
        """Test _EntityMapper class"""
        c = IMP.pmi.mmcif._EntityMapper()
        c.add('foo', 'MELS')
        c.add('bar', 'SELM')
        c.add('foo_2', 'MELS')
        self.assertEqual(c['foo'].id, 1)
        self.assertEqual(c['foo_2'].id, 1)
        self.assertEqual(c['bar'].id, 2)
        a = c.get_all()
        self.assertEqual(len(a), 2)
        self.assertEqual(a[0].id, 1)
        self.assertEqual(a[0].first_component, 'foo')
        self.assertEqual(a[0].description, 'foo')
        self.assertEqual(a[0].sequence, 'MELS')
        self.assertEqual(a[1].id, 2)
        self.assertEqual(a[1].first_component, 'bar')
        self.assertEqual(a[1].description, 'bar')
        self.assertEqual(a[1].sequence, 'SELM')

    def test_dataset_dumper_all_group(self):
        """Test DatasetDumper.get_all_group()"""
        dump = IMP.pmi.mmcif.DatasetDumper(EmptyObject())
        ds1 = IMP.pmi.mmcif.EM2DClassDataset()
        ds2 = IMP.pmi.mmcif.CXMSDataset()
        ds3 = IMP.pmi.mmcif.PDBDataset('1abc', '1.0', 'test details')

        g = dump.get_all_group()
        self.assertEqual(g.id, 1)
        self.assertEqual(g.datasets, [])

        dump.add(ds1)
        dump.add(ds2)
        g = dump.get_all_group()
        self.assertEqual(g.id, 2)
        self.assertEqual(g.datasets, [ds1, ds2])
        g = dump.get_all_group()
        self.assertEqual(g.id, 2)

        dump.add(ds3)
        g = dump.get_all_group()
        self.assertEqual(g.id, 3)
        self.assertEqual(g.datasets, [ds1, ds2, ds3])

    def test_dataset_dumper_duplicates(self):
        """Check that DatasetDumper ignores duplicate datasets"""
        dump = IMP.pmi.mmcif.DatasetDumper(EmptyObject())
        ds1 = dump.add(IMP.pmi.mmcif.PDBDataset('1abc', '1.0', 'test details'))
        self.assertEqual(ds1.id, 1)
        # A duplicate dataset should be ignored even if details differ
        ds2 = dump.add(IMP.pmi.mmcif.PDBDataset('1abc', '1.0', 'other details'))
        self.assertEqual(ds2.id, 1)
        self.assertEqual(id(ds1), id(ds2))

        loc1 = IMP.pmi.mmcif.DBDatasetLocation("mydb", "abc", "1.0", "")
        loc2 = IMP.pmi.mmcif.DBDatasetLocation("mydb", "xyz", "1.0", "")

        # Identical datasets in the same location aren't duplicated
        cx1 = IMP.pmi.mmcif.CXMSDataset()
        cx1.location = loc1
        cx2 = IMP.pmi.mmcif.CXMSDataset()
        cx2.location = loc1
        ds3 = dump.add(cx1)
        ds4 = dump.add(cx1)
        self.assertEqual(ds3.id, 2)
        self.assertEqual(ds4.id, 2)

        # Datasets in different locations are OK
        cx3 = IMP.pmi.mmcif.CXMSDataset()
        cx3.location = loc2
        ds5 = dump.add(cx3)
        self.assertEqual(ds5.id, 3)

        # Different datasets in same location are OK (but odd)
        em2d = IMP.pmi.mmcif.EM2DClassDataset()
        em2d.location = loc2
        ds6 = dump.add(em2d)
        self.assertEqual(ds6.id, 4)

    def test_dataset_dumper_dump(self):
        """Test DatasetDumper.dump()"""
        dump = IMP.pmi.mmcif.DatasetDumper(EmptyObject())
        ds = dump.add(IMP.pmi.mmcif.PDBDataset('1abc', '1.0', 'test details'))
        self.assertEqual(ds.location.access_code, '1abc')

        fh = StringIO()
        w = IMP.pmi.mmcif.CifWriter(fh)
        dump.dump(w)
        out = fh.getvalue()
        self.assertEqual(out, """#
loop_
_ihm_dataset_list.ordinal_id
_ihm_dataset_list.id
_ihm_dataset_list.group_id
_ihm_dataset_list.data_type
_ihm_dataset_list.database_hosted
1 1 1 'Experimental model' YES
#
#
loop_
_ihm_dataset_related_db_reference.id
_ihm_dataset_related_db_reference.dataset_list_id
_ihm_dataset_related_db_reference.db_name
_ihm_dataset_related_db_reference.access_code
_ihm_dataset_related_db_reference.version
_ihm_dataset_related_db_reference.data_type
_ihm_dataset_related_db_reference.details
1 1 PDB 1abc 1.0 'Experimental model' 'test details'
#
""")

    def test_model_dumper_sphere(self):
        """Test ModelDumper sphere_obj output"""
        class DummyPO(IMP.pmi.mmcif.ProtocolOutput):
            def flush(self):
                pass

        m = IMP.Model()
        simo = IMP.pmi.representation.Representation(m)
        po = DummyPO(None)
        simo.add_protocol_output(po)
        simo.create_component("Nup84")
        simo.add_component_sequence("Nup84",
                                    self.get_input_file_name("test.fasta"))
        nup84 = simo.autobuild_model("Nup84",
                                     self.get_input_file_name("test.nup84.pdb"),
                                     "A")

        d = IMP.pmi.mmcif.ModelDumper(po)
        self.assertEqual(d.add(simo.prot), 1)
        fh = StringIO()
        w = IMP.pmi.mmcif.CifWriter(fh)
        d.dump(w)
        out = fh.getvalue()
        self.assertEqual(out, """#
loop_
_ihm_sphere_obj_site.ordinal_id
_ihm_sphere_obj_site.entity_id
_ihm_sphere_obj_site.seq_id_begin
_ihm_sphere_obj_site.seq_id_end
_ihm_sphere_obj_site.asym_id
_ihm_sphere_obj_site.Cartn_x
_ihm_sphere_obj_site.Cartn_y
_ihm_sphere_obj_site.Cartn_z
_ihm_sphere_obj_site.object_radius
_ihm_sphere_obj_site.model_id
1 1 1 1 A 0.000 0.000 0.000 3.068 1
2 1 2 2 A 0.000 0.000 0.000 2.997 1
3 1 3 4 A 0.000 0.000 0.000 3.504 1
#
""")

    def test_chem_comp_dumper(self):
        """Test ChemCompDumper"""
        class DummyPO(IMP.pmi.mmcif.ProtocolOutput):
            def flush(self):
                pass

        po = DummyPO(None)
        po.create_component("Nup84")
        po.add_component_sequence("Nup84", "MELS")
        po.create_component("Nup85")
        po.add_component_sequence("Nup85", "MC")

        d = IMP.pmi.mmcif.ChemCompDumper(po)

        fh = StringIO()
        w = IMP.pmi.mmcif.CifWriter(fh)
        d.dump(w)
        out = fh.getvalue()
        self.assertEqual(out, """#
loop_
_chem_comp.id
MET
GLU
LEU
SER
CYS
#
""")

if __name__ == '__main__':
    IMP.test.main()
