/**
 *  \file IMP/pmi1/MembraneRestraint.h
 *  \brief Favor configurations where target is in the membrane.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 *
 */

#ifndef IMPPMI1_MEMBRANE_RESTRAINT_H
#define IMPPMI1_MEMBRANE_RESTRAINT_H
#include <IMP/pmi1/pmi1_config.h>
#include <IMP/isd/ISDRestraint.h>
#include <IMP/Particle.h>

IMPPMI1_BEGIN_NAMESPACE
//! Membrane Restraint
/** Favors configurations where target is in the membrane
 */
class IMPPMI1EXPORT MembraneRestraint : public isd::ISDRestraint {
  ParticleIndex z_nuisance;
  double thickness;
  double softness;
  double plateau;
  double linear;
  double max_float;
  double log_max_float;

  ParticleIndexes particles_below;
  ParticleIndexes particles_above;
  ParticleIndexes particles_inside;

public:
  MembraneRestraint(Model *m, ParticleIndex z_nuisance, double thickness,
                    double softness, double plateau, double linear);
  void add_particles_below(ParticleIndexes const &particles);
  void add_particles_above(ParticleIndexes const &particles);
  void add_particles_inside(ParticleIndexes const &particles);
  double get_score(double prob) const;
  double get_probability_above(double z, double z_slope_center_upper) const;
  double get_score_above(double z, double z_slope_center_upper) const;
  double get_probability_below(double z, double z_slope_center_lower) const;
  double get_score_below(double z, double z_slope_center_lower) const;
  double get_score_inside(double z, double z_slope_center_lower,
                      double z_slope_center_upper) const;
  virtual double unprotected_evaluate(DerivativeAccumulator *) const override;
  virtual IMP::ModelObjectsTemp do_get_inputs() const override;
  IMP_OBJECT_METHODS(MembraneRestraint);
};

IMPPMI1_END_NAMESPACE

#endif /* IMPPMI1_MEMBRANE_RESTRAINT_H */
