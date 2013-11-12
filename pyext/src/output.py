class Output():
    
    def __init__(self,ascii=True):
        global os,RMF,imprmf,cPickle,impatom,impcore,imp
        import cPickle as cPickle
        import os
        import IMP as imp
        import IMP.atom as impatom
        import IMP.core as impcore
        try:
            import RMF
            import IMP.rmf as imprmf
            self.rmf_library=True
        except ImportError:
            self.rmf_library=False


        self.dictionary_pdbs={}
        self.dictionary_rmfs={}
        self.dictionary_stats={}
        self.dictionary_stats2={}
        self.best_score_list=None
        self.nbestscoring=None
        self.suffix=None
        self.ascii=ascii
        self.initoutput={}
        self.residuetypekey=imp.StringKey("ResidueName")
        self.chainids="ABCDEFGHIJKLMNOPQRSTUVXYWZabcdefghijklmnopqrstuvxywz"
        self.dictchain={}

    def get_pdb_names(self):
        return self.dictionary_pdbs.keys()

    def get_rmf_names(self):
        return self.dictionary_rmfs.keys()

    def get_stat_names(self):
        return self.dictionary_stats.keys()

    def init_pdb(self,name,prot):
        flpdb=open(name,'w')
        flpdb.close()
        self.dictionary_pdbs[name]=prot
        self.dictchain[name]={}
        
        for n,i in enumerate(self.dictionary_pdbs[name].get_children()):
            self.dictchain[name][i.get_name()]=self.chainids[n]

    def write_pdb(self,name,appendmode=True):
        if appendmode:
            flpdb=open(name,'a')
        else:
            flpdb=open(name,'w')
            
        #impatom.write_pdb(self.dictionary_pdbs[name],flpdb)
        
        for n,p in enumerate(impatom.get_leaves(self.dictionary_pdbs[name])):
        
           if p.get_parent().get_name() in self.dictchain[name]:
              protname=p.get_parent().get_name()
           else:
              p0=p.get_parent()
              protname=p0.get_parent().get_name()
        
           resind=impatom.Fragment(p).get_residue_indexes()
        
           if len(resind)==1:
           
           
              flpdb.write(impatom.get_pdb_string(impcore.XYZ(p).get_coordinates(),
                             n,impatom.AT_CA,impatom.ResidueType((p.get_value(self.residuetypekey))),
                             self.dictchain[name][protname],resind[0]))
        flpdb.write("ENDMOL\n")
        
        
        flpdb.close()

    def write_pdbs(self,appendmode=True):
        for pdb in self.dictionary_pdbs.keys():
            self.write_pdb(pdb,appendmode)

    def init_pdb_best_scoring(self,suffix,prot,nbestscoring):
        # save only the nbestscoring conformations
        # create as many pdbs as needed
        self.best_score_list=[]
        self.suffix=suffix
        self.nbestscoring=nbestscoring
        for i in range(self.nbestscoring):
            name=suffix+"."+str(i)+".pdb"
            flpdb=open(name,'w')
            flpdb.close()
            self.dictionary_pdbs[name]=prot
            self.dictchain[name]={}
            for n,i in enumerate(self.dictionary_pdbs[name].get_children()):
                self.dictchain[name][i.get_name()]=self.chainids[n]

    def write_pdb_best_scoring(self,score):
        if self.nbestscoring==None:
            print "Output.write_pdb_best_scoring: init_pdb_best_scoring not run"

        #update the score list
        if len(self.best_score_list)<self.nbestscoring:
            self.best_score_list.append(score)
            self.best_score_list.sort()
            index=self.best_score_list.index(score)
            for i in range(len(self.best_score_list)-2,index-1,-1):
                oldname=self.suffix+"."+str(i)+".pdb"
                newname=self.suffix+"."+str(i+1)+".pdb"
                os.rename(oldname, newname)
            filetoadd=self.suffix+"."+str(index)+".pdb"
            self.write_pdb(filetoadd,appendmode=False)

        else:
            if score<self.best_score_list[-1]:
                self.best_score_list.append(score)
                self.best_score_list.sort()
                self.best_score_list.pop(-1)
                index=self.best_score_list.index(score)
                for i in range(len(self.best_score_list)-1,index-1,-1):
                    oldname=self.suffix+"."+str(i)+".pdb"
                    newname=self.suffix+"."+str(i+1)+".pdb"
                    os.rename(oldname, newname)
                filenametoremove=self.suffix+"."+str(self.nbestscoring)+".pdb"
                os.remove(filenametoremove)
                filetoadd=self.suffix+"."+str(index)+".pdb"
                self.write_pdb(filetoadd,appendmode=False)

    def init_rmf(self,name,prot):
        if not self.rmf_library:
            print "Output error: neet rmf library to init rmf"
            exit()

        rh = RMF.create_rmf_file(name)
        imprmf.add_hierarchy(rh, prot)
        self.dictionary_rmfs[name]=rh

    def add_restraints_to_rmf(self,name,objectlist):
        for o in objectlist:
            try:
               rs=o.get_restraint_for_rmf()
            except:
               rs=o.get_restraint()
            imprmf.add_restraints(self.dictionary_rmfs[name],rs.get_restraints())

    def add_geometries_to_rmf(self,name,objectlist):
        for o in objectlist:
            geos=o.get_geometries()
            imprmf.add_geometries(self.dictionary_rmfs[name],geos)

    
    def add_particle_pair_from_restraints_to_rmf(self,name,objectlist):
        for o in objectlist:
            print "here"
            
            pps=o.get_particle_pairs()
            for pp in pps:
              print type(IMP.core.EdgePairGeometry(pp))
              imprmf.add_geometry(self.dictionary_rmfs[name],IMP.core.EdgePairGeometry(pp))  

    def write_rmf(self,name,nframe):
        imprmf.save_frame(self.dictionary_rmfs[name],nframe)
        self.dictionary_rmfs[name].flush()

    def close_rmf(self,name):
        del self.dictionary_rmfs[name]

    def write_rmfs(self,nframe):
        for rmf in self.dictionary_rmfs.keys():
            self.write_rmf(rmf,nframe)

    def init_stat(self,name,listofobjects):
        if self.ascii:
            flstat=open(name,'w')
            flstat.close()
        else:
            flstat=open(name,'wb')
            flstat.close()

        #check that all objects in listofobjects have a  get_output method

        for l in listofobjects:
            if not "get_output" in dir(l):
                print "Output: object ", l, " doesn't have get_output() method"
                exit()
        self.dictionary_stats[name]=listofobjects

    def set_output_entry(self,key,value):
        self.initoutput.update({key:value})

    def write_stat(self,name,appendmode=True):
        output=self.initoutput
        for obj in self.dictionary_stats[name]:
            d=obj.get_output()
            #remove all entries that begin with _ (private entries)
            dfiltered=dict((k, v) for k, v in d.iteritems() if k[0]!="_")
            output.update(dfiltered)

        if appendmode:
            writeflag='a'
        else:
            writeflag='w'

        if self.ascii:
            flstat=open(name,writeflag)
            flstat.write("%s \n" % output)
            flstat.close()
        else:
            flstat=open(name,writeflag+'b')
            cPickle.dump(output,flstat,2)
            flstat.close()

    def write_stats(self):
        for stat in self.dictionary_stats.keys():
            self.write_stat(stat)

    def get_stat(self,name):
        output={}
        for obj in self.dictionary_stats[name]:
            output.update(obj.get_output())
        return output

    def write_test(self,name,listofobjects):
        '''
        write the test:
        output=output.Output()
        output.write_test("test_modeling11_models.rmf_45492_11Sep13_veena_imp-020713.dat",outputobjects)
        run the test:
        output=output.Output()        
        output.test("test_modeling11_models.rmf_45492_11Sep13_veena_imp-020713.dat",outputobjects)
        '''
        flstat=open(name,'w')    
        output=self.initoutput
        for l in listofobjects:
            if not "get_output" in dir(l):
                print "Output: object ", l, " doesn't have get_output() method"
                exit()
        self.dictionary_stats[name]=listofobjects
        for obj in self.dictionary_stats[name]:
            d=obj.get_output()
            #remove all entries that begin with _ (private entries)
            dfiltered=dict((k, v) for k, v in d.iteritems() if k[0]!="_")
            output.update(dfiltered)
        output.update({"ENVIRONMENT":str(self.get_environment_variables())})
        output.update({"IMP_VERSIONS":str(self.get_versions_of_relevant_modules())})             
        flstat.write("%s \n" % output)
        flstat.close()        
        

    def test(self,name,listofobjects):
        from numpy.testing import assert_approx_equal as aae
        output=self.initoutput
        for l in listofobjects:
            if not "get_output" in dir(l):
                print "Output: object ", l, " doesn't have get_output() method"
                exit()
        for obj in listofobjects:
            output.update(obj.get_output())
        output.update({"ENVIRONMENT":str(self.get_environment_variables())})
        output.update({"IMP_VERSIONS":str(self.get_versions_of_relevant_modules())})   

        flstat=open(name,'r')  
        for l in flstat:
            test_dict=eval(l)
        for k in test_dict:
            if k in output:
               if test_dict[k]!=output[k]: print str(k)+": test failed, old value: "+str(test_dict[k])+" new value "+str(output[k])
               #aae(float(test_dict[k]),
               #    float(output[k]),7,str(k)+": test failed, old value: "+str(test_dict[k])+" new value "+str(output[k]))
            else:
               print str(k)+" from old objects (file "+str(name)+") not in new objects"
               
    def get_environment_variables(self):
        import os 
        return str(os.environ)
    
    def get_versions_of_relevant_modules(self):
        import IMP
        versions={}        
        versions["IMP_VERSION"]=IMP.kernel.get_module_version()     
        try:
           import IMP.pmi
           versions["PMI_VERSION"]=IMP.pmi.get_module_version() 
        except (ImportError):
           pass                      
        try:
           import IMP.isd2
           versions["ISD2_VERSION"]=IMP.isd2.get_module_version()             
        except (ImportError):
           pass         
        return versions

       
    

#-------------------
      
    def init_stat2(self,name,listofobjects,extralabels=None,listofsummedobjects=None):
        #this is a new stat file that should be less 
        #space greedy!
        #listofsummedobjects must be in the form [([obj1,obj2,obj3,obj4...],label)]
        #extralabels
        
        if listofsummedobjects==None: listofsummedobjects=[]
        if extralabels==None: extralabels=[]
        flstat=open(name,'w')
        output={}
        stat2_keywords={"STAT2HEADER":"STAT2HEADER"}
        stat2_keywords.update({"STAT2HEADER_ENVIRON":str(self.get_environment_variables())})
        stat2_keywords.update({"STAT2HEADER_IMP_VERSIONS":str(self.get_versions_of_relevant_modules())})        
        stat2_inverse={}
        
        for l in listofobjects:
            if not "get_output" in dir(l):
                print "Output: object ", l, " doesn't have get_output() method"
                exit()
            else:
                d=l.get_output()
                #remove all entries that begin with _ (private entries)
                dfiltered=dict((k, v) for k, v in d.iteritems() if k[0]!="_")
                output.update(dfiltered)
        
        #check for customizable entries     
        for l in listofsummedobjects:
          for t in l[0]:
            if not "get_output" in dir(t):
                print "Output: object ", t, " doesn't have get_output() method"
                exit()
            else:
                if "_TotalScore" not in t.get_output():
                  print "Output: object ", t, " doesn't have _TotalScore entry to be summed"
                  exit()
                else:
                  output.update({l[1]:0.0})            
        
        for k in extralabels:
            output.update({k:0.0})  
        
        for n,k in enumerate(output):
            stat2_keywords.update({n:k})
            stat2_inverse.update({k:n})
        
        flstat.write("%s \n" % stat2_keywords)
        flstat.close()
        self.dictionary_stats2[name]=(listofobjects,stat2_inverse,listofsummedobjects,extralabels)

    def write_stat2(self,name,appendmode=True):
        output={}
        (listofobjects,stat2_inverse,listofsummedobjects,extralabels)=self.dictionary_stats2[name]

        #writing objects
        for obj in listofobjects:
            od=obj.get_output()
            dfiltered=dict((k, v) for k, v in od.iteritems() if k[0]!="_")
            for k in dfiltered: 
               output.update({stat2_inverse[k]:od[k]})
        
        #writing summedobjects
        for l in listofsummedobjects:
           partial_score=0.0
           for t in l[0]:
             d=t.get_output()
             partial_score+=float(d["_TotalScore"])
           output.update({stat2_inverse[l[1]]:str(partial_score)})
        
        #writing extralabels
        for k in extralabels:
           if k in self.initoutput:
              output.update({stat2_inverse[k]:self.initoutput[k]})
           else:
              output.update({stat2_inverse[k]:"None"})              
        
        if appendmode:
            writeflag='a'
        else:
            writeflag='w'

        flstat=open(name,writeflag)
        flstat.write("%s \n" % output)
        flstat.close()

    def write_stats2(self):
        for stat in self.dictionary_stats2.keys():
            self.write_stat2(stat)

class ProcessOutput():
    
    def __init__(self,filename):
        self.filename=filename
        self.isstat1=False
        self.isstat2=False
        
        #open the file
        if self.filename!=None:
           f=open(self.filename,"r")
        else:
           print "Error: No file name provided. Use -h for help"
           exit()
        
        #get the keys from the first line
        for line in f.readlines():
            d=eval(line)
            self.klist=d.keys()
            #check if it is a stat2 file
            if "STAT2HEADER" in self.klist: 
                import operator
                self.isstat2=True
                for k in self.klist:
                    if "STAT2HEADER" in str(k):
                       #if print_header: print k, d[k]
                       del d[k]
                stat2_dict=d
                #get the list of keys sorted by value
                kkeys=[k[0] for k in sorted(stat2_dict.iteritems(), key=operator.itemgetter(1))]
                self.klist=[k[1] for k in sorted(stat2_dict.iteritems(), key=operator.itemgetter(1))]
                self.invstat2_dict={}
                for k in kkeys:
                    self.invstat2_dict.update({stat2_dict[k]:k})
            else:
                self.isstat1=True
                self.klist.sort()
                
            break
        f.close()    
            
    def get_keys(self):      
        return self.klist
        
    def get_fields(self,fields):
           
           
           outdict={}
           for field in fields:
               outdict[field]=[]
           
           #print fields values
           f=open(self.filename,"r")
           line_number=0
           for line in f.readlines():
              line_number+=1
              try:
                 d=eval(line)
              except:
                 print "# Warning: skipped line number " + str(line_number) + " not a valid line"
                 continue
              if   self.isstat1: [outdict[field].append(d[field]) for field in fields]
              elif self.isstat2: 
                   if line_number==1: continue
                   [outdict[field].append(d[self.invstat2_dict[field]]) for field in fields]
           f.close()
           return outdict        

    def plot_fields(self,fields):
        import matplotlib.pyplot as plt
        
        plt.rc('lines', linewidth=4)
        fig, axs  = plt.subplots(nrows=len(fields))
        fig.set_size_inches(10.5,5.5*len(fields))
        plt.rc('axes', color_cycle=['r'])
        
        n=0
        for key in fields:
           x = range(len(fields[key]))
           y=[float(y) for y in fields[key]]
           if len(fields)>1:
              axs[n].plot(x,y)
              axs[n].set_title(key)
           else:
              axs.plot(x,y)
              axs.set_title(key)
           
           n+=1

        # Tweak spacing between subplots to prevent labels from overlapping
        plt.subplots_adjust(hspace=0.3)
        plt.show()     


    def plot_field_histogram(self,name,values,valuename=None):
        import matplotlib.pyplot as plt
        plt.hist([float(y) for y in values],bins=40,color='#66CCCC',normed=True)
        plt.title(name)
        if valuename==None:
           plt.xlabel(name)
        else:
           plt.xlabel(valuename)
        plt.ylabel("Frequency")
        plt.savefig(name+".png",dpi=150,transparent="True")
        plt.show()
    
    
def plot_fields_box_plots(name,values,positions,
                          valuename="None",positionname="None"):
    '''
    This function plots time series as boxplots
    fields is a list of time series, positions are the x-values
    valuename is the y-label, positionname is the x-label
    '''
    import matplotlib.pyplot as plt
    from matplotlib.patches import Polygon
    #import numpy as np
    
    bps=[]
    fig = plt.figure(figsize=(10,10))
    fig.canvas.set_window_title(name)
    ax1 = fig.add_subplot(111)
    
    plt.subplots_adjust(left=0.161, right=0.990, top=0.95, bottom=0.11)
    

    bps.append(plt.boxplot( values, notch=0, sym='', vert=1, 
                              whis=1.5,positions=positions))  
    
    plt.setp(bps[-1]['boxes'], color='black',lw=1.5)
    plt.setp(bps[-1]['whiskers'], color='black',ls=":",lw=1.5)
    
    plt.xlabel(positionname)
    plt.ylabel(valuename)
    plt.savefig(name+".png",dpi=150,transparent="True")
    plt.show()


def plot_xy_data(x,y):
        import matplotlib.pyplot as plt
        plt.rc('lines', linewidth=4)
        fig, ax  = plt.subplots(nrows=1)
        fig.set_size_inches(10.5,5.5)
        plt.rc('axes', color_cycle=['r'])
        print x
        print y
        ax.plot(x,y)
        plt.show()  


def get_graph_from_hierarchy(hier):
    graph=[]
    depth_dict={}
    depth=0
    (graph,depth,depth_dict)=recursive_graph(hier,graph,depth,depth_dict)
    
    #filters node labels according to depth_dict
    node_labels_dict={}
    node_size_dict={}
    for key in depth_dict:
        node_size_dict=10/depth_dict[key]
        if depth_dict[key]<3: 
           node_labels_dict[key]=key
        else:
           node_labels_dict[key]=""
    draw_graph(graph,labels_dict=node_labels_dict)



def recursive_graph(hier,graph,depth,depth_dict):
    depth=depth+1
    nameh=IMP.atom.Hierarchy(hier).get_name()
    index=str(hier.get_particle().get_index())
    name1=nameh+"|#"+index
    depth_dict[name1]=depth
    
    children=IMP.atom.Hierarchy(hier).get_children()
    
    if len(children)==1 or children==None:
       depth=depth-1
       return (graph,depth,depth_dict)
    
    else:
      for c in children:
        (graph,depth,depth_dict)=recursive_graph(c,graph,depth,depth_dict)
        nameh=IMP.atom.Hierarchy(c).get_name()
        index=str(c.get_particle().get_index())
        namec=nameh+"|#"+index        
        graph.append((name1,namec))
      
      depth=depth-1
      return (graph,depth,depth_dict)


     
def draw_graph(graph, labels_dict=None, graph_layout='spring',
               node_size=5, node_color='blue', node_alpha=0.3,
               node_text_size=11,
               edge_color='blue', edge_alpha=0.3, edge_tickness=1,
               edge_text_pos=0.3,
               text_font='sans-serif'):

    import networkx as nx
    import matplotlib.pyplot as plt


    # create networkx graph
    G=nx.Graph()

    # add edges
    for edge in graph:
        G.add_edge(edge[0], edge[1])

    # these are different layouts for the network you may try
    # shell seems to work best
    if graph_layout == 'spring':
        graph_pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        graph_pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        graph_pos=nx.random_layout(G)
    else:
        graph_pos=nx.shell_layout(G)

    # draw graph
    nx.draw_networkx_nodes(G,graph_pos,node_size=node_size, 
                           alpha=node_alpha, node_color=node_color)
    nx.draw_networkx_edges(G,graph_pos,width=edge_tickness,
                           alpha=edge_alpha,edge_color=edge_color)
    nx.draw_networkx_labels(G, graph_pos,labels=labels_dict,font_size=node_text_size,
                            font_family=text_font)
    plt.show()

    
    
def draw_table():
    
    #still an example!
    
    from ipyD3 import d3object
    from IPython.display import display
    
    d3 = d3object(width=800,
              height=400,
              style='JFTable',
              number=1,
              d3=None,
              title='Example table with d3js',
              desc='An example table created created with d3js with data generated with Python.')
    data=[[1277.0, 654.0, 288.0, 1976.0, 3281.0, 3089.0, 10336.0, 4650.0, 4441.0, 4670.0, 944.0, 110.0],
    [1318.0, 664.0, 418.0, 1952.0, 3581.0, 4574.0, 11457.0, 6139.0, 7078.0, 6561.0, 2354.0, 710.0],
    [1783.0, 774.0, 564.0, 1470.0, 3571.0, 3103.0, 9392.0, 5532.0, 5661.0, 4991.0, 2032.0, 680.0],
    [1301.0, 604.0, 286.0, 2152.0, 3282.0, 3369.0, 10490.0, 5406.0, 4727.0, 3428.0, 1559.0, 620.0],
    [1537.0, 1714.0, 724.0, 4824.0, 5551.0, 8096.0, 16589.0, 13650.0, 9552.0, 13709.0, 2460.0, 720.0],
    [5691.0, 2995.0, 1680.0, 11741.0, 16232.0, 14731.0, 43522.0, 32794.0, 26634.0, 31400.0, 7350.0, 3010.0],
    [1650.0, 2096.0, 60.0, 50.0, 1180.0, 5602.0, 15728.0, 6874.0, 5115.0, 3510.0, 1390.0, 170.0],
    [72.0, 60.0, 60.0, 10.0, 120.0, 172.0, 1092.0, 675.0, 408.0, 360.0, 156.0, 100.0]]
    data=[list(i) for i in zip(*data)]
    sRows=[['January', 'February', 'March', 'April', 'May', 'June', 'July', 'August', 'September', 'October', 'November', 'Deecember']]
    sColumns=[['Prod {0}'.format(i) for i in xrange(1,9)],
          [None, '', None, None, 'Group 1', None, None, 'Group 2']]
    d3.addSimpleTable(   data, 
                     fontSizeCells=[12,],
                     sRows=sRows,
                     sColumns=sColumns,
                     sRowsMargins=[5,50,0],
                     sColsMargins=[5,20,10],
                     spacing=0,
                     addBorders=1,
                     addOutsideBorders=-1,
                     rectWidth=45,
                     rectHeight=0                   
                 )
    html=d3.render(mode=['html', 'show'])
    display(html)
