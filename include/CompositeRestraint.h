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
//#include <IMP/container/CloseBipartitePairContainer.h>

IMPPMI_BEGIN_NAMESPACE
/** A restraint for ambiguous cross-linking MS data and multiple state approach.
    It marginalizes the false positive rate and depends on the expected fpr and
    an uncertainty parameter beta.
 */

class IMPPMIEXPORT  CompositeRestraint : public Restraint
{
    //particle indexes in the composite
    IMP::kernel::ParticleIndexes pis_;
    IMP::kernel::ParticleIndex handle_particle_index_;
    double coffd_;
    double l_;
    //base::map<std::tuple<unsigned int,unsigned int>,
    //          base::Pointer<container::CloseBipartitePairContainer>> map_cont_;

public:


  //! Create the restraint.
  /** Restraints should store the particles they are to act on,
      preferably in a Singleton or PairContainer as appropriate.
   */

  CompositeRestraint(IMP::kernel::Model *m, 
                     IMP::kernel::ParticleIndexAdaptor handle_particle_index, 
                     double coffd, double l, 
                     std::string name="CompositeRestraint%1%");

  void add_composite_particle(IMP::kernel::ParticleIndexAdaptor pi){pis_.push_back(pi);}
  
        
  unsigned int get_number_of_particles() const {return pis_.size();}  

  /* call for probability */
  double get_probability_per_particle_excluding(unsigned int ipart, 
                                    std::vector<unsigned int> excluded_ps ) const;
  
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
