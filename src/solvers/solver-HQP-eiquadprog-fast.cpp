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

#include "tsid/solvers/solver-HQP-eiquadprog-fast.hpp"
#include "tsid/math/utils.hpp"
#include "tsid/solvers/eiquadprog-fast.hpp"
#include "tsid/utils/stop-watch.hpp"

//#define PROFILE_EIQUADPROG_FAST



namespace tsid
{
  namespace solvers
  {
    
    using namespace math;
    SolverHQuadProgFast::SolverHQuadProgFast(const std::string & name):
    SolverHQPBase(name),
    m_hessian_regularization(DEFAULT_HESSIAN_REGULARIZATION)
    {
      m_n = 0;
      m_neq = 0;
      m_nin = 0;
    }
    
    void SolverHQuadProgFast::sendMsg(const std::string & s)
    {
      std::cout<<"[SolverHQuadProgFast."<<m_name<<"] "<<s<<std::endl;
    }

    void SolverHQuadProgFast::resize(unsigned int n, unsigned int neq, unsigned int nin)
    {
      hcod.resize(nVar,p);
      J.resize(p);
      b.resize(p);
    }
    
    const HQPOutput & SolverHQuadProgFast::solve(const HQPData & problemData)
    {
      START_PROFILER_EIQUADPROG_FAST(PROFILE_EIQUADPROG_PREPARATION);

      // assign the constraint matrices
      /* SOTH resize */
      p = problemData.size();

      if (p > 0)
      {
        nVar = problemData[0][0].second->cols();
      }
      else
          return m_output;

      hcod.resize(nVar,p);

      int l=0;
      for (const auto & cl : problemData)
      {
        unsigned int neq = 0, nin = 0;
        for(ConstraintLevel::const_iterator it=cl.begin(); it!=cl.end(); it++)
        {
          const ConstraintBase* constr = it->second;
          assert(nVar==constr->cols());
          if(constr->isEquality())
            neq += constr->rows();
          else
            nin += constr->rows();
        }
        // If necessary, resize the constraint matrices
        // resize(1,1,1);

        J[l].resize(neq+nin,m_n);
        J[l].topRows(neq) = m_CI;
        J[l].topRows(nin) = m_CE;
        b[l].resize(neq + nin);

        int i_eq=0, i_in=0;
        for(ConstraintLevel::const_iterator it=cl.begin(); it!=cl.end(); it++)
        {
          const ConstraintBase* constr = it->second;
          if(constr->isEquality())
          {
            m_CE.middleRows(i_eq, constr->rows()) = constr->matrix();
            m_ce0.segment(i_eq, constr->rows())   = -constr->vector();
            J[l].block(i_eq,0,constr->rows(),m_n) = constr->matrix();
            for (int c=0;c<constr->rows();c++)
            {
                b[l][i_eq+c] = soth::Bound(-constr->vector()[c],-constr->vector()[c]);
            }
            i_eq += constr->rows();
        }
          else if(constr->isInequality())
          {
            m_CI.middleRows(i_in, constr->rows()) = constr->matrix();
            m_ci0.segment(i_in, constr->rows())   = -constr->lowerBound();
            i_in += constr->rows();
            m_CI.middleRows(i_in, constr->rows()) = -constr->matrix();
            m_ci0.segment(i_in, constr->rows())   = constr->upperBound();
            J[l].block(neq+i_in,0,constr->rows(),m_n) = -constr->matrix();
            for (int c=0;c<constr->rows();c++)
            {
                b[l][neq+i_in+c] = soth::Bound(-constr->lowerBound()[c],constr->upperBound()[c]);
            }
            i_in += constr->rows();
          }
          else if(constr->isBound())
          {
            m_CI.middleRows(i_in, constr->rows()).setIdentity();
            m_ci0.segment(i_in, constr->rows())   = -constr->lowerBound();
            i_in += constr->rows();
            m_CI.middleRows(i_in, constr->rows()) = -Matrix::Identity(m_n, m_n);
            m_ci0.segment(i_in, constr->rows())   = constr->upperBound();
            J[l].block(neq+i_in,0,constr->rows(),m_n) = -Matrix::Identity(m_n,m_n);
            for (int c=0;c<constr->rows();c++)
            {
                b[l][neq+i_in+c] = soth::Bound(-constr->lowerBound()[c],constr->upperBound()[c]);
            }
            i_in += constr->rows();
          }
        }

        std::cout<<"J["<<l<<"]:\n"<<J[l]<<std::endl;
        hcod.pushBackStage(J[l],b[l]);
        l += 1;
      }


      /* solve HCOD */
      hcod.setNameByOrder("stage_");
      hcod.useDamp(false);
      // const double dampingFactor = 0.0;
      // hcod.setDamping(dampingFactor);
      // hcod.stage(0).damping(0);
      hcod.setInitialActiveSet();

      // if (memory.iteration > 0) {
      //     hcod.setInitialActiveSet(memory.hcodAS);
      // } else {
      //     // hcod.setInitialActiveSet(memory.hcodAS);
      //     hcod.setInitialActiveSet();
      // }
      //   // }
      // hcod.initialize();

      Eigen::VectorXd hcodSolution = Eigen::VectorXd::Zero(m_n);
      hcod.activeSearch(hcodSolution);
      std::cout<<"hcod solution: "<<hcodSolution.transpose()<<std::endl;

      //   // m_output.m_Kinv = m_output.m_K.inverse();
      //   // compute_slack(problemData, m_output);
      
      return m_output;
    }
    
    double SolverHQuadProgFast::getObjectiveValue()
    {
      return m_solver.getObjValue();
    }

    const HQPOutput & SolverHQuadProgFast::resolve(const HQPData & problemData)
    {
        hcod.resetBounds(b);
        Eigen::VectorXd hcodSolution = Eigen::VectorXd::Zero(m_n);
        hcod.activeSearch(hcodSolution);

    }
    
    bool SolverHQuadProgFast::setMaximumIterations(unsigned int maxIter)
    {
      SolverHQPBase::setMaximumIterations(maxIter);
      return m_solver.setMaxIter(maxIter);
    }
  }
}



