//
// Copyright (c) 2017 CNRS
//
// This file is part of tsid
// tsid is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
// tsid is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// tsid If not, see
// <http://www.gnu.org/licenses/>.
//

#ifndef __invdyn_solvers_hqp_output_hpp__
#define __invdyn_solvers_hqp_output_hpp__

#include "tsid/solvers/fwd.hpp"
#include "tsid/math/fwd.hpp"

#include <vector>


namespace tsid
{
  namespace solvers
  {
    
    class HQPOutput
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      
      typedef math::Vector Vector;
      typedef math::Matrix Matrix;
      typedef math::VectorXi VectorXi;
      
      HQPStatus status;    /// solver status
      Vector x;            /// solution
      Vector lambda;       /// Lagrange multipliers
      VectorXi activeSet;  /// indexes of active inequalities
      Vector activeSetPy;  /// indexes of active inequalities (dtype double)
      Matrix m_A;          /// active constraint normal
      Vector m_b;          /// active constraint bias
      Matrix m_H;          /// Matrices in the KKT equation
      Vector m_g;          /// Hx* + A'lambda = -g
      Matrix m_K;          /// KKT matrix K = [[H A'];[A 0]]
      Matrix m_Kinv;       /// inverse of the KKT matrix
      Vector m_delta;      /// constraint measures: equality -> norm, inequality -> smallest margin
      int iterations;      /// number of iterations performed by the solver
      
      HQPOutput(){}
      
      HQPOutput(int nVars, int nEqCon, int nInCon)
      {
        resize(nVars, nEqCon, nInCon);
      }
      
      void resize(int nVars, int nEqCon, int nInCon)
      {
        x.resize(nVars);
        lambda.resize(nEqCon+nInCon);
        activeSet.resize(nInCon);
      }
    };
  }
}

#endif // ifndef __invdyn_solvers_hqp_output_hpp__
