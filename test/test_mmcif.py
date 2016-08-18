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
        po = DummyPO(None)
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
2 1 bar 1 B 1 3
3 2 bar 1 B 1 3
4 2 baz 2 C 1 2
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

if __name__ == '__main__':
    IMP.test.main()
