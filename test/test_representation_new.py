#A new representation module. It helps to construct the hierarchy
#and deal with multi-state, multi-scale, multi-copies

#Usage Example:

import IMP
import IMP.atom
import IMP.pmi
import IMP.pmi.representation_new as r
import IMP.test


class RepresentationNewTest(IMP.test.TestCase):
    
    def test_system_base(self):
        class TestSystemBase(r._SystemBase):
           def __init__(self):
             self.model=IMP.Model()
             
        sb=TestSystemBase()
        root=sb._create_hierarchy()
        child=sb._create_child(root)
        self.assertEqual(child.get_parent(),root)
        sb.build()
    
    def test_create_states(self):
        s=r.System()
        for i in range(10):
            # create a new state
            self.assertEqual(s.get_number_of_states(),i)
            st=s.create_state()
            self.assertEqual(st.state_hierarchy.get_parent(),s.system_hierarchy)
            self.assertEqual(st.model,s.model)
        self.assertEqual(s.get_number_of_states(),10) 
 

    def test_create_molecules(self):
        s=r.System()
        sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st=s.create_state()
        p1=st.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st.create_molecule("Prot2",sequence=sequence_prot2)            
        self.assertEqual(p1.molecule_hierarchies[0].get_parent(),st.state_hierarchy)
        self.assertEqual(p2.molecule_hierarchies[0].get_parent(),st.state_hierarchy) 
        self.assertEqual(p1.model,st.model)                       
        self.assertEqual(p2.model,st.model)                  
        self.assertEqual(p1.molecule_name,"Prot1")                       
        self.assertEqual(p2.molecule_name,"Prot2")              
        self.assertEqual(len(p1.molecule_hierarchies),1)
        self.assertEqual(len(p2.molecule_hierarchies),1) 
        self.assertEqual(len(st.state_hierarchy.get_children()),2)        
        st=s.create_state()
        p3=st.create_molecule("Prot3",sequence=sequence_prot3)         
        self.assertEqual(p3.molecule_hierarchies[0].get_parent(),st.state_hierarchy)
        self.assertEqual(p3.model,st.model)                
        self.assertEqual(p3.molecule_name,"Prot3")            
        self.assertEqual(len(p3.molecule_hierarchies),1) 
        self.assertEqual(len(st.state_hierarchy.get_children()),1)  

    def test_create_copies(self):
        s=r.System()
        sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st1.create_molecule("Prot2",sequence=sequence_prot2) 
        st2=s.create_state()        
        p3=st2.create_molecule("Prot3",sequence=sequence_prot3)  
        p1.add_copy()
        p1.add_copy()
        p1.add_copy()        
        p2.add_copy()

        self.assertEqual(p1.get_number_of_copies(),len(p1.molecule_hierarchies))        
        self.assertEqual(p2.get_number_of_copies(),len(p2.molecule_hierarchies))
        self.assertEqual(p3.get_number_of_copies(),len(p3.molecule_hierarchies))
        self.assertEqual(p1.get_number_of_copies(),len(p1.molecule_hierarchies))        
        self.assertEqual(p2.get_number_of_copies(),len(p2.molecule_hierarchies))
        self.assertEqual(p3.get_number_of_copies(),len(p3.molecule_hierarchies))
        self.assertEqual(p1.get_number_of_copies(),4)   
        self.assertEqual(p2.get_number_of_copies(),2) 
        self.assertEqual(p3.get_number_of_copies(),1)
        self.assertEqual(len(st1.state_hierarchy.get_children()),p1.get_number_of_copies()+p2.get_number_of_copies())        
        self.assertEqual(len(st2.state_hierarchy.get_children()),p3.get_number_of_copies()) 

    def test_add_structure(self):
        # incomplete    
        s=r.System()
        sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st1.create_molecule("Prot2",sequence=sequence_prot2) 
        st2=s.create_state()        
        p3=st2.create_molecule("Prot3",sequence=sequence_prot3)  
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        
    def test_set_representation(self):
        # incomplete    
        s=r.System()
        sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st1.create_molecule("Prot2",sequence=sequence_prot2) 
        st2=s.create_state()        
        p3=st2.create_molecule("Prot3",sequence=sequence_prot3)  
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        p1.set_representation()
        p2.set_representation()
        p3.set_representation()

    def test_build_system(self):
        # incomplete    
        s=r.System()
        sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st1.create_molecule("Prot2",sequence=sequence_prot2) 
        st2=s.create_state()        
        p3=st2.create_molecule("Prot3",sequence=sequence_prot3)  
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        p1.set_representation()
        p2.set_representation()
        p3.set_representation()
        s.build()

    def test_sequence(self):
        # incomplete
        sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        

    def test_build_partial1(self):
        # incomplete    
        s=r.System()
        sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st1.create_molecule("Prot2",sequence=sequence_prot2) 
        st2=s.create_state()        
        p3=st2.create_molecule("Prot3",sequence=sequence_prot3)  
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        p1.set_representation()
        p2.set_representation()
        p3.set_representation()
        st1.build()
        st2.build()

    def test_build_partial2(self):
        # incomplete    
        s=r.System()
        sequence_prot1=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot1")
        sequence_prot2=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot2")
        sequence_prot3=r.Sequence(fastafile="my_fasta_file.fasta",id_fastafile="Prot3")
        # create a new state
        st1=s.create_state()
        p1=st1.create_molecule("Prot1",sequence=sequence_prot1)
        p2=st1.create_molecule("Prot2",sequence=sequence_prot2) 
        st2=s.create_state()        
        p3=st2.create_molecule("Prot3",sequence=sequence_prot3)  
        p1.add_copy()
        p1.add_structure()
        p2.add_structure()
        p3.add_structure()
        p1.set_representation()
        p2.set_representation()
        p3.set_representation()
        p1.build()
        p2.build()
        p3.build()

if __name__ == '__main__':
    IMP.test.main()
