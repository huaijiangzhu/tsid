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

#ifndef __invdyn_solvers_hqp_eiquadprog_fast_hpp__
#define __invdyn_solvers_hqp_eiquadprog_fast_hpp__

#include "tsid/solvers/solver-HQP-base.hpp"
#include "tsid/solvers/eiquadprog-fast.hpp"
#include <soth/HCOD.hpp>
#include <lexls/lexlsi.h>

namespace tsid
{
  namespace solvers
  {
    /**
     * @brief
     */
    class TSID_DLLAPI SolverHQuadProgFast : public SolverHQPBase
    {
    public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
      
      typedef math::Matrix Matrix;
      typedef math::Vector Vector;
      typedef math::RefVector RefVector;
      typedef math::ConstRefVector ConstRefVector;
      typedef math::ConstRefMatrix ConstRefMatrix;

      SolverHQuadProgFast(const std::string & name);

      void resize(unsigned int n, unsigned int neq, unsigned int nin);
      void initializeSolver(const HQPData & problemData);
      void reinitializeSolver();

      /** Solve the given Hierarchical Quadratic Program
       */
      const HQPOutput & solve(const HQPData & problemData);
      /** Kai hacks
       */
      /** set the initial active set for warm start
       */
      void setInitialActiveSet(const std::vector<soth::cstref_vector_t> & activeSetIn);
      /** get the active set
       */
      const std::vector<soth::cstref_vector_t> getActiveSet();
      /** Resolve the given Hierarchical Quadratic Program with changed right hand side
       */
      const HQPOutput & resolve(const HQPData & problemData);
      const std::vector<Eigen::MatrixXd>& getLagrangeMultipliers();


      /** Get the objective value of the last solved problem. */
      double getObjectiveValue();

      /** Set the current maximum number of iterations performed by the solver. */
      bool setMaximumIterations(unsigned int maxIter);

      void compute_slack(const HQPData & problemData,
        HQPOutput & problemOutput);

      const HQPOutput & solve_local(const HQPData & problemData,
                              const HQPOutput & previousOutput);

    protected:

      void sendMsg(const std::string & s);

      EiquadprogFast m_solver; // <nVars, nEqCon, 2*nIneqCon>

      Matrix m_H;
      Vector m_g;
      Matrix m_CE;
      Vector m_ce0;
      Matrix m_CI;  /// twice the rows because inequality constraints are bilateral
      Vector m_ci0;
      double m_objValue;

      double m_hessian_regularization;

      Eigen::VectorXi m_activeSet;  /// vector containing the indexes of the active inequalities
      int m_activeSetSize;

      unsigned int m_neq;  /// number of equality constraints
      unsigned int m_nin;  /// number of inequality constraints
      unsigned int m_n;    /// number of variables

      // soth::HCOD hcod;
      int p;
      std::vector<Matrix> J;
      std::vector<soth::VectorBound> b;
      // sizes of the hierarchy levels
      struct level {
          unsigned int m_neq = 0;
          unsigned int m_nin = 0;
      };
      std::vector<level> hLvl;
      std::vector<soth::cstref_vector_t> activeSet;
      std::vector<Eigen::MatrixXd> lambda;

      // lexls
      std::vector<Eigen::MatrixXd> A;

      bool coldStart = true;

    };
  }
}

#endif // ifndef __invdyn_solvers_hqp_eiquadprog_fast_hpp__
