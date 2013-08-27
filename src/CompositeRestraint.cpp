/**
 *  \file pmi/CompositeRestraint.h
 *  \brief A sigmoid shaped restraint between
 *  residues with discrete classifier
 *  and ambiguous assignment. To be used with
 *  cross-linking mass-spectrometry data.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#include <IMP/pmi/CompositeRestraint.h>
#include <IMP/core/XYZR.h>
#include <math.h>
#include <limits.h>


IMPPMI_BEGIN_NAMESPACE

CompositeRestraint::CompositeRestraint(IMP::kernel::Model *m, 
                          IMP::kernel::ParticleIndexAdaptor handle_particle_index,
                          double coffd, double l, std::string name):
                          Restraint(m, name), 
                          handle_particle_index_(handle_particle_index), 
                          coffd_(coffd), l_(l) {pis_.push_back(handle_particle_index_);}                       


double CompositeRestraint::
                 unprotected_evaluate(DerivativeAccumulator *accum) const
{
    double score=0;
    //std::cout << "here" << std::endl;
    Cache cache;
    CachePot cachepot;    
    Ints excluded_ps;
    excluded_ps.push_back(0);
    double prob=get_probability_per_particle_excluding(0,excluded_ps,cache,cachepot);
    if (prob==0.0){score += std::numeric_limits<double>::max( );}
    else{score+=-log(prob);}
    
    if (accum){};
    
    return score;
}
              
double CompositeRestraint::get_probability_per_particle_excluding(unsigned int ipart, 
                                    Ints excluded_ps, Cache& cache, CachePot& cachepot ) const
{   
    std::sort(excluded_ps.begin(),excluded_ps.end());
    CacheKey key(ipart,excluded_ps);
    if (cache.find(key)!=cache.end()){ 
      return cache.find(key)->second;
    }
    
    IMP::kernel::ParticleIndex ppi=pis_[ipart];  
    core::XYZR di(get_model(), ppi);
    double onemprob=1.0;
    
      
    for(unsigned int k=0;k<get_number_of_particles();++k){
      //std::cout << k << std::endl;
      if (std::find(excluded_ps.begin(), excluded_ps.end(), k) == excluded_ps.end()){
         //std::cout << "ex " << k << std::endl;         
         
         //check that the probability term was already calculated 
         //before
         excluded_ps.push_back(k);
         CacheKeyPot keypot1(ipart,k);

         
         double p;
         
         if (cachepot.find(keypot1)!=cachepot.end()){ 
            p=cachepot.find(keypot1)->second;
         }         
         //otherwise look in the cache of precalculated terms
         else {
            IMP::kernel::ParticleIndex ppk=pis_[k];               
            core::XYZR dk(get_model(), ppk);
            double dist = (di.get_coordinates() -
                dk.get_coordinates()).get_magnitude();
            //p =1.0/(1.0+std::exp(-(coffd_-dist)/l_));
            p=1.0;
            CacheKeyPot keypot2(k,ipart);
            cachepot[keypot1]=p;
            cachepot[keypot2]=p;            
         }
         
         
         if (excluded_ps.size()==get_number_of_particles())
         {
         onemprob *= 1.0-p;
         }
         else
         {
         onemprob *= 1.0-p*get_probability_per_particle_excluding(k,excluded_ps,cache,cachepot);
         }
         excluded_ps.pop_back(); 
      }
    }

  double prob=1.0-onemprob;
  cache[key]=prob;
  return prob;
}





/* Return all particles whose attributes are read by the restraints. To
   do this, ask the pair score what particles it uses.*/
ModelObjectsTemp  CompositeRestraint::do_get_inputs() const
{
  ParticlesTemp ret;
  for(unsigned int k=0;k<get_number_of_particles();++k){
     ret.push_back(get_model()->get_particle(pis_[k]));
  }
  return ret;
}

IMPPMI_END_NAMESPACE
