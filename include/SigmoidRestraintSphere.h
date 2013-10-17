/**
 *  \file pmi/SigmoidRestraintSphere.h
 *  \brief Simple Sigmoidal score calculated between sphere surfaces.
 *
 *  Copyright 2007-2013 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPPMI_SIGMOID_RESTRAINT_SPHERE_H
#define IMPPMI_SIGMOID_RESTRAINT_SPHERE_H
#include "pmi_config.h"
#include <IMP/Restraint.h>
#include <IMP/restraint_macros.h>


IMPPMI_BEGIN_NAMESPACE
/** Blah blah.
 */

class IMPPMIEXPORT  SigmoidRestraintSphere : public Restraint
{
    IMP::kernel::ParticleIndex p1_;
    IMP::kernel::ParticleIndex p2_;
    double inflection_;
    double slope_;
    double amplitude_;



public:


  //! Create the restraint.
  /** Restraints should store the particles they are to act on,
      preferably in a Singleton or PairContainer as appropriate.
   */

  SigmoidRestraintSphere(IMP::kernel::Model *m, 
                          IMP::kernel::ParticleIndexAdaptor p1,
                          IMP::kernel::ParticleIndexAdaptor p2,
                          double inflection, double slope, 
                          double amplitude,
                     std::string name="SigmoidRestraintSphere%1%");

  void set_amplitude(double amplitude){amplitude_=amplitude;}
  void increment_amplitude(double amplitude){amplitude_=amplitude_+amplitude;}  
  double get_amplitude(){return amplitude_;}

  /** This macro declares the basic needed methods: evaluate and show
   */
  virtual double
  unprotected_evaluate(IMP::kernel::DerivativeAccumulator *accum)
     const IMP_OVERRIDE;
  virtual IMP::kernel::ModelObjectsTemp do_get_inputs() const IMP_OVERRIDE;
  IMP_OBJECT_METHODS(SigmoidRestraintSphere);


};

IMPPMI_END_NAMESPACE

#endif  /* IMPPMI_SIGMOID_RESTRAINT_SPHERE_H */
