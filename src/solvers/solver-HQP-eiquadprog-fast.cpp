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
      // FIXME: make hcod restartable
      hcod = soth::HCOD(m_n,p);
      hcod.setNameByOrder("stage_");
      hcod.useDamp(false);
      hcod.setInitialActiveSet(activeSet);
      hLvl.resize(p);
      J.resize(p);
      b.resize(p);
    }

    void SolverHQuadProgFast::setInitialActiveSet(const std::vector<soth::cstref_vector_t> & activeSetIn) {
    }

    const std::vector<soth::cstref_vector_t> SolverHQuadProgFast::getActiveSet() {
        return activeSet;
    }
    
    const HQPOutput & SolverHQuadProgFast::solve(const HQPData & problemData)
    {
      START_PROFILER_EIQUADPROG_FAST(PROFILE_EIQUADPROG_PREPARATION);

      // assign the constraint matrices
      /* SOTH resize */
      p = problemData.size();

      // number of variables
      if (p > 0)
      {
        m_n = problemData[0][0].second->cols();
        resize(0,0,0);
      }
      else
          return m_output;


      int l=0;
      for (const auto & cl : problemData)
      {
        for(ConstraintLevel::const_iterator it=cl.begin(); it!=cl.end(); it++)
        {
          const ConstraintBase* constr = it->second;
          assert(m_n==constr->cols());
          if(constr->isEquality())
            hLvl[l].m_neq += constr->rows();
          else
            hLvl[l].m_nin += constr->rows();
        }

        J[l].resize(hLvl[l].m_neq + hLvl[l].m_nin,m_n);
        b[l].resize(hLvl[l].m_neq + hLvl[l].m_nin);

        int i_eq=0, i_in=0;
        for(ConstraintLevel::const_iterator it=cl.begin(); it!=cl.end(); it++)
        {
          const ConstraintBase* constr = it->second;
          if(constr->isEquality())
          {
            J[l].block(i_eq,0,constr->rows(),m_n) = constr->matrix();
            for (int c=0;c<constr->rows();c++)
            {
                b[l][i_eq+c] = soth::Bound(constr->vector()[c],constr->vector()[c]); // FIXME: double bound in order to get lagrange multipliers / slack for equality constraints
            }
            i_eq += constr->rows();
          }
          else if(constr->isInequality())
          {
            J[l].block(hLvl[l].m_neq+i_in,0,constr->rows(),m_n) = constr->matrix();
            for (int c=0;c<constr->rows();c++)
            {
                b[l][hLvl[l].m_neq+i_in+c] = soth::Bound(constr->lowerBound()[c],constr->upperBound()[c]);
            }
            i_in += constr->rows();
          }
          else if(constr->isBound())
          {
            J[l].block(hLvl[l].m_neq+i_in,0,constr->rows(),m_n).setIdentity(); // FIXME: need to set the correct variable to 1
            for (int c=0;c<constr->rows();c++)
            {
                b[l][hLvl[l].m_neq+i_in+c] = soth::Bound(constr->lowerBound()[c],constr->upperBound()[c]);
            }
            i_in += constr->rows();
          }
        }

        // std::cout << "tsid::J["<<l<<"]:\n"<<J[l]<<std::endl;
        // std::cout << "tsid::b["<<l<<"]:\n"<<b[l]<<std::endl;

        hcod.pushBackStage(J[l],b[l]);
        l += 1;
      }

      // resize output
      m_output.resize(m_n, hLvl[0].m_neq, 2*hLvl[0].m_nin);

      /* solve HCOD */
      hcod.activeSearch(m_output.x);
      activeSet = hcod.getOptimalActiveSet();

      /* assign rest of m_output */
      // m_output.lambda = hcod.getLagrangeMultipliers();
      // . TODO
      // store active set in output
      int cr=0, actIneqCtr=0;
      m_output.activeSet.resize(hLvl[0].m_nin);
      m_output.activeSetPy.resize(hLvl[0].m_nin);
      for(ConstraintLevel::const_iterator it=problemData[0].begin(); it!=problemData[0].end(); it++)
      {
        const ConstraintBase* constr = it->second;
        for (int cr=0;cr<constr->rows();cr++)
        {
            if (activeSet[0][cr].type > 0)
            {
                if(constr->isInequality() || constr->isBound())
                {
                  m_output.activeSet(actIneqCtr) = activeSet[0][cr].type;
                  m_output.activeSetPy(actIneqCtr) = activeSet[0][cr].type;
                  actIneqCtr++;
                }
            }
        }
        cr += constr->rows();
      }

      compute_slack(problemData, m_output);

      return m_output;
    }
    
    double SolverHQuadProgFast::getObjectiveValue()
    {
      return m_solver.getObjValue();
    }

    const HQPOutput & SolverHQuadProgFast::resolve(const HQPData & problemData)
    {
      // get the changed rhs
      int l=0;
      for (const auto & cl : problemData)
      {
        int i_eq=0, i_in=0;
        for(ConstraintLevel::const_iterator it=cl.begin(); it!=cl.end(); it++)
        {
          const ConstraintBase* constr = it->second;
          if(constr->isEquality())
          {
            for (int c=0;c<constr->rows();c++)
            {
                b[l][i_eq+c] = soth::Bound(constr->vector()[c],constr->vector()[c]);
            }
            i_eq += constr->rows();
          }
          else if(constr->isInequality())
          {
            for (int c=0;c<constr->rows();c++)
            {
                b[l][hLvl[l].m_neq+i_in+c] = soth::Bound(constr->lowerBound()[c],constr->upperBound()[c]);
            }
            i_in += constr->rows();
          }
          else if(constr->isBound())
          {
            for (int c=0;c<constr->rows();c++)
            {
                b[l][hLvl[l].m_neq+i_in+c] = soth::Bound(constr->lowerBound()[c],constr->upperBound()[c]);
            }
            i_in += constr->rows();
          }
        }
        l += 1;
      }

      // reset bounds and resolve
      hcod.resetBounds(b);
      // continue active search with 1 iteration
      hcod.activeSearch_cont(m_output.x,1);

      compute_slack(problemData, m_output);

      return m_output;

    }

    const std::vector<Eigen::MatrixXd>& SolverHQuadProgFast::getLagrangeMultipliers() {
        return hcod.getLagrangeMultipliers();
    }
    
    bool SolverHQuadProgFast::setMaximumIterations(unsigned int maxIter)
    {
      SolverHQPBase::setMaximumIterations(maxIter);
      return m_solver.setMaxIter(maxIter);
    }

    void SolverHQuadProgFast::compute_slack(const HQPData & problemData, 
                                            HQPOutput & problemOutput) {
      const ConstraintLevel & cl0 = problemData[0];
      const Vector & x = problemOutput.x;
      problemOutput.m_slack.resize(p);
      for (int l=0;l<p;l++) {
        problemOutput.m_slack[l] = hcod.getLagrangeMultipliers()[l].col(l).norm();
      }
    }

    const HQPOutput & SolverHQuadProgFast::solve_local(const HQPData & problemData,
                                                       const HQPOutput & previousOutput) {
    }
  }
}



