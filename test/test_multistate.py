# -*- coding: utf-8 -*-
# <nbformat>3.0</nbformat>

# <codecell>

import IMP
import IMP.pmi.representation
import IMP.pmi.restraints.basic

# <codecell>

representations=[]

# <codecell>

m=IMP.Model()
r=IMP.pmi.representation.Representation(m)

r.create_component("particle1",color=0.1)
p11=r.add_component_beads("particle1",[(1,10)])
r.create_component("particle2",color=0.5)
p21=r.add_component_beads("particle2",[(1,10)])
r.create_component("particle3",color=0.9)
p31=r.add_component_beads("particle3",[(1,10)])

representations.append(r)

r=IMP.pmi.representation.Representation(m)

r.create_component("particle1",color=0.1)
p12=r.add_component_beads("particle1",[(1,10)])
r.create_component("particle2",color=0.5)
p22=r.add_component_beads("particle2",[(1,10)])
r.create_component("particle3",color=0.9)
p32=r.add_component_beads("particle3",[(1,10)])

representations.append(r)

representations[0].floppy_bodies.pop(0)
representations[0].floppy_bodies.pop(0)
representations[1].floppy_bodies.pop(0)
representations[1].floppy_bodies.pop(0)

print representations[0].floppy_bodies
print representations[1].floppy_bodies

# <codecell>

import IMP.core 

pp1=IMP.atom.get_leaves(p11[0])[0]
pp2=IMP.atom.get_leaves(p21[0])[0]
pp3=IMP.atom.get_leaves(p31[0])[0]
xyz11=IMP.core.XYZ(pp1.get_particle())
xyz21=IMP.core.XYZ(pp2.get_particle())
xyz31=IMP.core.XYZ(pp3.get_particle())
xyz11.set_coordinates((0,0,0))
print xyz11.get_coordinates()
xyz21.set_coordinates((80,0,0))
xyz31.set_coordinates((0,0,0))

pp1=IMP.atom.get_leaves(p12[0])[0]
pp2=IMP.atom.get_leaves(p22[0])[0]
pp3=IMP.atom.get_leaves(p32[0])[0]
xyz12=IMP.core.XYZ(pp1.get_particle())
xyz22=IMP.core.XYZ(pp2.get_particle())
xyz32=IMP.core.XYZ(pp3.get_particle())
xyz12.set_coordinates((0,0,0))
xyz22.set_coordinates((80,0,0))
xyz32.set_coordinates((80,0,0))

# <codecell>

eb=IMP.pmi.restraints.basic.ExternalBarrier(representations[0],100)
eb.add_to_model()

eb=IMP.pmi.restraints.basic.ExternalBarrier(representations[1],100)
eb.add_to_model()

# <codecell>

import IMP.pmi.restraints.crosslinking

restraints='''#
particle2 particle3 1 5 1 1
particle1 particle3 1 2 1 2 '''

xl = IMP.pmi.restraints.crosslinking.ISDCrossLinkMS(representations, 
                                   restraints,
                                   length=10.0,
                                   slope=0.0,
 				                   inner_slope=0.01,
                                   resolution=1.0)

psi=xl.get_psi(0.05)

psi[0].set_scale(0.1)

sigma=xl.get_sigma(1.0)

sigma[0].set_scale(5.0)

#psi.set_scale(0.1)
#xl.get_sigma(1.0)

out_dict=xl.get_output()
sorted_keys = sorted(out_dict.keys())
for entry in sorted_keys:
   print entry,out_dict[entry]

print xyz11.get_coordinates()    

xl.add_to_model()
print m.evaluate(False)

# <codecell>

import IMP.pmi.samplers

print xyz11.get_coordinates()

mc = IMP.pmi.samplers.MonteCarlo(m,representations, 1.0)
mc.set_simulated_annealing(min_temp=1.0,
                           max_temp=2.0,
                           min_temp_time=200,
                           max_temp_time=50)

# <codecell>

import IMP.pmi.output

print xyz32.get_coordinates()

o = IMP.pmi.output.Output()
o.init_stat2("modeling.stat",[mc,xl]+representations)

for i in range(0,10000):
   mc.optimize(10)
   o.write_stats2()
   if i%100==0: print i

# <codecell>

import IMP.pmi.output

po=IMP.pmi.output.ProcessOutput("modeling.stat")

# <codecell>

po.get_keys()

# <codecell>

fs=po.get_fields(['ISDCrossLinkMS_Distance_interrb-State:0-1:particle1_2:particle3-1-1-0.05_None',
             'ISDCrossLinkMS_Distance_interrb-State:0-1:particle2_5:particle3-1-1-0.05_None',
             'ISDCrossLinkMS_Distance_interrb-State:1-1:particle1_2:particle3-1-1-0.05_None',
             'ISDCrossLinkMS_Distance_interrb-State:1-1:particle2_5:particle3-1-1-0.05_None',
              'ISDCrossLinkMS_Data_Score_None',
              'ISDCrossLinkMS_Linear_Score_None',
              'ISDCrossLinkMS_Psi_0.05_None'])

# <codecell>

#% matplotlib inline

# <codecell>

IMP.pmi.output.plot_scatter_xy_data(fs['ISDCrossLinkMS_Distance_interrb-State:0-1:particle1_2:particle3-1-1-0.05_None'],
                                    fs['ISDCrossLinkMS_Distance_interrb-State:1-1:particle1_2:particle3-1-1-0.05_None'])

IMP.pmi.output.plot_scatter_xy_data(fs['ISDCrossLinkMS_Distance_interrb-State:0-1:particle2_5:particle3-1-1-0.05_None'],
                                    fs['ISDCrossLinkMS_Distance_interrb-State:0-1:particle1_2:particle3-1-1-0.05_None'])

IMP.pmi.output.plot_scatter_xy_data(fs['ISDCrossLinkMS_Distance_interrb-State:1-1:particle2_5:particle3-1-1-0.05_None'],
                                    fs['ISDCrossLinkMS_Distance_interrb-State:1-1:particle1_2:particle3-1-1-0.05_None'])


IMP.pmi.output.plot_scatter_xy_data(fs['ISDCrossLinkMS_Distance_interrb-State:0-1:particle1_2:particle3-1-1-0.05_None'],
                                    fs['ISDCrossLinkMS_Distance_interrb-State:1-1:particle2_5:particle3-1-1-0.05_None'])

# <codecell>

IMP.pmi.output.plot_fields(fs)

# <codecell>

