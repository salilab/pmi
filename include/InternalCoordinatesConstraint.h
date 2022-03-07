/**
 *  \file IMP/pmi1/InternalCoordinatesConstraint.h
 *  \brief A modifier that enforce identical internal coordinates for NonRigidMembers.
 *
 *  Copyright 2007-2018 IMP Inventors. All rights reserved.
 */

#ifndef IMPPMI1_INTERNAL_COORDINATES_CONSTRAINT_H
#define IMPPMI1_INTERNAL_COORDINATES_CONSTRAINT_H

#include <IMP/pmi1/pmi1_config.h>
#include <IMP/SingletonModifier.h>
#include <IMP/singleton_macros.h>

IMPPMI1_BEGIN_NAMESPACE

//! An example singleton modifier
/** A simple singleton modifier which ensures the x,y,z coordinates
    stay in a box by wrapping them.

    Such a class could be coupled with an IMP::core::SingletonConstraint
    or IMP::container::SingletonsConstraint to keep a set of particles
    in a box.

    \include range_restrictionExampleSingletonModifier.py

    The source code is as follows:
    \include InternalCoordinatesConstraint.h
    \include InternalCoordinatesConstraint.cpp
 */
class IMPPMI1EXPORT InternalCoordinatesConstraint : public SingletonModifier {

 public:
  InternalCoordinatesConstraint();

  // note, Doxygen wants a semicolon at the end of macro lines
  virtual void apply_index(Model *m, ParticleIndex p) const override;
  virtual ModelObjectsTemp do_get_inputs(
      Model *m, const ParticleIndexes &pis) const override;
  virtual ModelObjectsTemp do_get_outputs(
      Model *m, const ParticleIndexes &pis) const override;
  IMP_SINGLETON_MODIFIER_METHODS(InternalCoordinatesConstraint);
  IMP_OBJECT_METHODS(InternalCoordinatesConstraint);
};

IMPPMI1_END_NAMESPACE

#endif /* IMPPMI1_INTERNAL_COORDINATES_CONSTRAINT_H */
