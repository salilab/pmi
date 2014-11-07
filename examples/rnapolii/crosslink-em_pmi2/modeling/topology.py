## \example setting up the system

import shared_functions

read_em=True
em_model_resolution=20.0

datadirectory="../data/"



# compname         hier_name     color    fastafile                      fastaid  pdbname                   chain resrange    read                          beadsize  density_rigid_body super_rigid_body emnum_components, emtxtfilename                               emmrcfilename
domains=[("Rpb1",  "Rpb1_1",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "A",  (1,1140),       5,       0,                 [0,1,2]),
         ("Rpb1",  "Rpb1_2",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "A",  (1141,1274),    5,       1,                 [0,1,2]),
         ("Rpb1",  "Rpb1_3",       0.0,     datadirectory+"1WCM.fasta.txt", "1WCM:A|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "A",  (1275,-1),      5,       0,                 [0,1,2]),
         ("Rpb2",  "Rpb2_1",       1.0,     datadirectory+"1WCM.fasta.txt", "1WCM:B|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "B",  (1,1102),       5,       2,                 [0,1,3]),
         ("Rpb2",  "Rpb2_2",       1.0,     datadirectory+"1WCM.fasta.txt", "1WCM:B|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "B",  (1103,-1),      5,       0,                 [0,1,3]),
         ("Rpb3",  "Rpb3",         0.2,     datadirectory+"1WCM.fasta.txt", "1WCM:C|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "C",  (1,-1),         5,       3,                 [0,4]),
         ("Rpb4",  "Rpb4",         0.3,     datadirectory+"1WCM.fasta.txt", "1WCM:D|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "D",  (1,-1),         5,       4,                 [0,5]),
         ("Rpb5",  "Rpb5",         0.4,     datadirectory+"1WCM.fasta.txt", "1WCM:E|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "E",  (1,-1),         5,       5,                 [0,6]),
         ("Rpb6",  "Rpb6",         0.5,     datadirectory+"1WCM.fasta.txt", "1WCM:F|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "F",  (1,-1),         5,       6,                 [0,7]),
         ("Rpb7",  "Rpb7",         0.6,     datadirectory+"1WCM.fasta.txt", "1WCM:G|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "G",  (1,-1),         5,       7,                 [0,8]),
         ("Rpb8",  "Rpb8",         0.7,     datadirectory+"1WCM.fasta.txt", "1WCM:H|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "H",  (1,-1),         5,       8,                 [0,9]),
         ("Rpb9",  "Rpb9",         0.8,     datadirectory+"1WCM.fasta.txt", "1WCM:I|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "I",  (1,-1),         5,       9,                 [0,10]),
         ("Rpb10", "Rpb10",        0.9,     datadirectory+"1WCM.fasta.txt", "1WCM:J|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "J",  (1,-1),         5,       10,                [0,11]),
         ("Rpb11", "Rpb11",        0.1,     datadirectory+"1WCM.fasta.txt", "1WCM:K|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "K",  (1,-1),         5,       11,                [0,12]),
         ("Rpb12", "Rpb12",        0.35,    datadirectory+"1WCM.fasta.txt", "1WCM:L|PDBID|CHAIN|SEQUENCE",  datadirectory+"1WCM.pdb", "L",  (1,-1),         5,       12,                [0,13])]



resdensities=[]
super_rigid_bodies={}
densities={}

for d in domains:
    comp_name        =d[0]
    hier_name        =d[1]
    color            =d[2]
    fasta_file       =d[3]
    fasta_id         =d[4]
    pdb_name         =d[5]
    chain_id         =d[6]
    res_range        =d[7]
    bead_size        =d[8]
    dens             =d[9]
    super_rb         =d[10]

    if comp_name not in simo.get_component_names():
        simo.create_component(comp_name,color=color)
        simo.add_component_sequence(comp_name,fasta_file,fasta_id)
        outhier=shared_functions.autobuild(simo,comp_name,pdb_name,chain_id,res_range,read=read_em,beadsize=bead_size,color=color)


    #print len(outhier),dens

    if not str(dens) in densities:
        #print "creating a new entry %s %s " % (str(dens),str(len(outhier)))
        densities[str(dens)]=[outhier,set(super_rb)]
    else:
        #print "adding a new entry %s %s " % (str(dens),str(len(outhier)))
        densities[str(dens)][0]+=outhier
        for k in super_rb:
            densities[str(dens)][1].add(k)




for f in densities:
    hier=densities[f][0]

    dens_hier=shared_functions.create_density(simo,str(f)+"_dens",hier,
                             datadirectory+"/densities/"+str(f)+".txt",
                             datadirectory+"/densities/"+str(f)+".mrc",
                             read_em)
    resdensities+=dens_hier
    simo.set_rigid_body_from_hierarchies(hier+dens_hier)

    super_rb=list(densities[f][1])
    for k in super_rb:
        if not k in super_rigid_bodies:
            super_rigid_bodies[k]=hier
        else:
            super_rigid_bodies[k]+=hier
        super_rigid_bodies[k]+=dens_hier


for c in simo.get_component_names():
    simo.setup_component_sequence_connectivity(c,resolution=1.0,scale=3.0)
    simo.setup_component_geometry(c)

for k in super_rigid_bodies:
    print k
    simo.set_super_rigid_body_from_hierarchies(super_rigid_bodies[k])

simo.set_floppy_bodies()
