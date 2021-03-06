// @HEADER
// ************************************************************************
//
//               Rapid Optimization Library (ROL) Package
//                 Copyright (2014) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact lead developers:
//              Drew Kouri   (dpkouri@sandia.gov) and
//              Denis Ridzal (dridzal@sandia.gov)
//
// ************************************************************************
// @HEADER

#ifndef ROL_CONSTRAINTFROMOBJECTIVE_H
#define ROL_CONSTRAINTFROMOBJECTIVE_H

#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_SingletonVector.hpp"

/** @ingroup func_group
    \class ROL::ConstraintFromObjective
    \brief Creates a constraint from an objective function and a offset value

    Example:  Suppose we have an objective function f(x) and we wish to impose,
	      e.g., a condition f(x)-offset = 0, then this class creates the
              scalar constraint c(x) = f(x)-offset 
*/


namespace ROL {

template<typename Real> 
class ConstraintFromObjective : public Constraint<Real> {
private:
  const Ptr<Objective<Real>> obj_;
  Ptr<Vector<Real>>          dualVector_;
  const Real                 offset_;
  bool                       isDualInitialized_;

public:
  ConstraintFromObjective( const Ptr<Objective<Real>> &obj, const Real offset = 0 );

  const Ptr<Objective<Real>> getObjective(void) const;

  void setParameter( const std::vector<Real> &param ) override;

  void update( const Vector<Real>& x, UpdateType type, int iter = -1 ) override;
  void update( const Vector<Real>& x, bool flag = true, int iter = -1 ) override;
  void value( Vector<Real>& c, const Vector<Real>& x, Real& tol ) override;
  void applyJacobian( Vector<Real>& jv, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) override;
  void applyAdjointJacobian( Vector<Real>& ajv, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) override;
  void applyAdjointHessian( Vector<Real>& ahuv, const Vector<Real>& u, const Vector<Real>& v, const Vector<Real>& x, Real& tol ) override;

private:
  Real getValue( const Vector<Real>& x ); 
  void setValue( Vector<Real>& x, Real val );

}; // ConstraintFromObjective

} // namespace ROL

#include "ROL_ConstraintFromObjective_Def.hpp"

#endif // ROL_CONSTRAINTFROMOBJECTIVE_H
