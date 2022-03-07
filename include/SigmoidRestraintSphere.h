/**
 *  \file IMP/pmi1/SigmoidRestraintSphere.h
 *  \brief Simple sigmoidal score calculated between sphere surfaces.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPPMI1_SIGMOID_RESTRAINT_SPHERE_H
#define IMPPMI1_SIGMOID_RESTRAINT_SPHERE_H
#include "pmi1_config.h"
#include <IMP/Restraint.h>
#include <IMP/particle_index.h>


IMPPMI1_BEGIN_NAMESPACE

//! Simple sigmoidal score calculated between sphere surfaces.
class IMPPMI1EXPORT  SigmoidRestraintSphere : public Restraint
{
    IMP::ParticleIndex p1_;
    IMP::ParticleIndex p2_;
    double inflection_;
    double slope_;
    double amplitude_;
    double line_slope_;



public:


  //! Create the restraint.
  SigmoidRestraintSphere(IMP::Model *m, 
                          IMP::ParticleIndexAdaptor p1,
                          IMP::ParticleIndexAdaptor p2,
                          double inflection, double slope, 
                          double amplitude, double line_slope_=0,
                     std::string name="SigmoidRestraintSphere%1%");

  void set_amplitude(double amplitude){amplitude_=amplitude;}
  void increment_amplitude(double amplitude){amplitude_=amplitude_+amplitude;}  
  double get_amplitude(){return amplitude_;}

  virtual double
  unprotected_evaluate(IMP::DerivativeAccumulator *accum) const override;
  virtual IMP::ModelObjectsTemp do_get_inputs() const override;
  IMP_OBJECT_METHODS(SigmoidRestraintSphere);


};

IMPPMI1_END_NAMESPACE

#endif  /* IMPPMI1_SIGMOID_RESTRAINT_SPHERE_H */
