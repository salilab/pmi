## \example shared functions to run the example

def get_bead_bits(hierarchy):
    pdbbits=[]
    beadbits=[]
    for h in hierarchy:
       if "_pdb" in h.get_name() and "Res:0" in h.get_parent().get_name() :pdbbits.append(h)
       if "_bead" in h.get_name()                                         :beadbits.append(h)
    return (pdbbits,beadbits)

def create_density(simo,compname,comphier,txtfilename,mrcfilename,read=True):
    #density generation for the EM restraint
   simo.create_component(compname)


   
   if read:
      outhier=simo.add_component_density(compname,
                               inputfile=txtfilename) # read what it was calculated before

   else:
      print comphier
      (pdbbits,beadbits)=get_bead_bits(comphier)
      for h in pdbbits:
          print h.get_name(), len(IMP.pmi.tools.get_residue_indexes(h))    
      number_of_residues=sum(len(IMP.pmi.tools.get_residue_indexes(h)) for h in pdbbits )      
      num_components=int(float(number_of_residues)/em_model_resolution)     
      print compname,num_components,number_of_residues,[h.get_name() for h in comphier]              
      outhier=simo.add_component_density(compname,
                               pdbbits,
                               num_components=num_components, # number of gaussian into which the simulated density is approximated
                               resolution=0,      # resolution that you want to calculate the simulated density
                               outputfile=txtfilename, # do the calculation
                               outputmap=mrcfilename,
                               multiply_by_total_mass=True) # do the calculation and output the mrc                                    
   return outhier

def autobuild(simo,comname,pdbname,chain,resrange,read=True,beadsize=5,color=0.0):
    if pdbname is not None:
      if resrange[-1]==-1: resrange=(resrange[0],len(simo.sequence_dict[comname]))
      if read==False:    
         print "res0",comname, resrange
         outhier=simo.autobuild_model(comname, 
                          pdbname=pdbname,
                          chain=chain,
                          resrange=resrange,
                          resolutions=[0,1,10],
                          color=color, 
                          missingbeadsize=beadsize)
      else:
         print "res1",comname, resrange
         outhier=simo.autobuild_model(comname, 
                          pdbname=pdbname,
                          chain=chain,
                          resrange=resrange,                          
                          resolutions=[1,10], 
                          color=color,
                          missingbeadsize=beadsize)                                 
      return outhier
    else:
      seq_len=len(simo.sequence_dict[comname])
      outhier=simo.add_component_necklace(comname,
                            begin=1,
                            end=seq_len,
                            length=beadsize)  
                            
    return outhier    

