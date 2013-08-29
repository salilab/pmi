/**
 *  \file pmi/CompositeRestraint.h
 *  \brief A pmf based likelihood function 
 *  with prior knowledge on the flase positive rate.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPPMI_COMPOSITE_RESTRAINT_H
#define IMPPMI_COMPOSITE_RESTRAINT_H
#include "pmi_config.h"
#include <IMP/Restraint.h>
#include <IMP/restraint_macros.h>
#include <IMP/base/tuple_macros.h>
#include <IMP/base/map.h>
//#include <IMP/container/CloseBipartitePairContainer.h>

IMPPMI_BEGIN_NAMESPACE
/** A restraint for ambiguous cross-linking MS data and multiple state approach.
    It marginalizes the false positive rate and depends on the expected fpr and
    an uncertainty parameter beta.
 */

class IMPPMIEXPORT  CompositeRestraint : public Restraint
{
    //particle indexes in the composite
    base::Vector<IMP::kernel::ParticleIndexes> pis_;
    IMP::kernel::ParticleIndexes handle_particle_indexes_;
    double coffd_;
    double l_;
    IMP_NAMED_TUPLE_2(CacheKey, CacheKeys,
                  Int, ipart, Ints, excluded, );

    IMP_NAMED_TUPLE_2(CacheKeyPot, CacheKeyPots,
                  Int, ipart, Int, kpart, );

    typedef base::map<CacheKey, double> Cache;
    typedef base::map<CacheKeyPot, double> CachePot;    
    
    //variables needed to tabulate the exponential
    Floats prob_grid_;
    double invdx_;
    double argmax_;
    double argmin_;
    bool tabprob_;

    inline double calc_prob (double dist) const{
      double argvalue=(dist-coffd_)/l_;  
      double prob;
      if (tabprob_){
         //this prevents something being below the lower value of the array
         double maxarg=std::max(argvalue,argmin_);      
         //this prevents something being above the upper value of the array
         double minarg=std::min(maxarg,argmax_);       
         unsigned k = static_cast<unsigned>( std::floor(minarg*invdx_) );
         prob=prob_grid_[k];
      }
      else{
         prob=1.0/(1.0+std::exp(-argvalue));
      }
      return prob;
    }

    
    //base::map<std::tuple<unsigned int,unsigned int>,
    //          base::Pointer<container::CloseBipartitePairContainer>> map_cont_;

  /* call for probability */
  double get_probability_per_particle_excluding(unsigned int ipart, 
                                    Ints excluded_ps, Cache& cache, CachePot& cachepot) const;

public:


  //! Create the restraint.
  /** Restraints should store the particles they are to act on,
      preferably in a Singleton or PairContainer as appropriate.
   */

  CompositeRestraint(IMP::kernel::Model *m, 
                     IMP::kernel::ParticleIndexesAdaptor handle_particle_indexes, 
                     double coffd, double l, bool tabprob, 
                     std::string name="CompositeRestraint%1%");

  void add_composite_particle(IMP::kernel::ParticleIndexesAdaptor pi){pis_.push_back(pi);}
  
        
  unsigned int get_number_of_elements() const {return pis_.size();}  


  
  //double get_probability() const {return 0.0;}

  /** This macro declares the basic needed methods: evaluate and show
   */
  virtual double
  unprotected_evaluate(IMP::kernel::DerivativeAccumulator *accum)
     const IMP_OVERRIDE;
  virtual IMP::kernel::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
  IMP_OBJECT_METHODS(CompositeRestraint);

  virtual double get_probability() const
  {
    return exp(-unprotected_evaluate(nullptr));
  }


};

IMPPMI_END_NAMESPACE

#endif  /* IMPPMI_COMPOSITE_RESTRAINT_H */
