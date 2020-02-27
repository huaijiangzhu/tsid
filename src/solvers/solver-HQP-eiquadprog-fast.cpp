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
    }

    void SolverHQuadProgFast::initializeSolver(const HQPData & problemData)
    {
      // assign the constraint matrices
      bool dimChange = false;
      /* SOTH resize */
      // FIXME: case of change of hierarchy levels
      p = problemData.size();
      assert(p > 0 && "ERROR: the task hierarchy is empty");

      if (coldStart) {
      	hLvl.resize(p);
      	J.resize(p);
      	b.resize(p);
      }

      if (m_n != problemData[0][0].second->cols()) {
        dimChange = true;
        m_n = problemData[0][0].second->cols();
      }

      int l=0;
      for (const auto & cl : problemData)
      {
	    unsigned int hLvllm_neq_prev = hLvl[l].m_neq;
	    unsigned int hLvllm_nin_prev = hLvl[l].m_nin;
        hLvl[l].m_neq = 0;
        hLvl[l].m_nin = 0;
        for(ConstraintLevel::const_iterator it=cl.begin(); it!=cl.end(); it++)
        {
          const ConstraintBase* constr = it->second;
          assert(m_n==constr->cols());
          if(constr->isEquality())
            hLvl[l].m_neq += constr->rows();
          else
            hLvl[l].m_nin += constr->rows();
        }

        if (hLvllm_neq_prev != hLvl[l].m_neq ||
	        hLvllm_nin_prev != hLvl[l].m_nin) {
            dimChange = true;
            J[l].resize(hLvl[l].m_neq + hLvl[l].m_nin,m_n);
            b[l].resize(hLvl[l].m_neq + hLvl[l].m_nin);
        }

        l += 1;
      }

      // number of variables
      if (coldStart || dimChange)
      {
        coldStart = true; // in case that coldStart = false and dimChange = true
      	hcod = soth::HCOD(m_n,p);
      	hcod.setNameByOrder("stage_");
      	hcod.useDamp(false);
      }
      else
	    reinitializeSolver();
    }

    void SolverHQuadProgFast::reinitializeSolver()
    {
	    hcod.reset();
    }

    void SolverHQuadProgFast::setInitialActiveSet(const std::vector<soth::cstref_vector_t> & activeSetIn) {
    }

    const std::vector<soth::cstref_vector_t> SolverHQuadProgFast::getActiveSet() {
        return activeSet;
    }
    
    const HQPOutput & SolverHQuadProgFast::solve(const HQPData & problemData)
    {
      START_PROFILER_EIQUADPROG_FAST(PROFILE_EIQUADPROG_PREPARATION);

      coldStart = true;
      initializeSolver(problemData);
 
      int l=0;
      for (const auto & cl : problemData)
      {
        int i_eq=0, i_in=0;
        for(ConstraintLevel::const_iterator it=cl.begin(); it!=cl.end(); it++)
        {
          const ConstraintBase* constr = it->second;
          if(constr->isEquality())
          {
            J[l].block(i_eq,0,constr->rows(),m_n) = constr->matrix();
            for (int c=0;c<constr->rows();c++)
            {
                b[l][i_eq+c] = soth::Bound(constr->vector()[c]); // FIXME: double bound in order to get lagrange multipliers / slack for equality constraints
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

        // std::cerr << "tsid::J["<<l<<"]:\n"<<J[l]<<std::endl;
        std::cerr << "tsid::b["<<l<<"]:\n"<<b[l]<<std::endl;

        if (coldStart) {
            hcod.pushBackStage(J[l],b[l]);
        } else
            hcod.resetStage(J[l],b[l],l);
        l += 1;
      }

      if (coldStart)
      {
      	hcod.setInitialActiveSet();
        // coldStart = false;
      }

      // resize output
      m_output.resize(m_n, hLvl[0].m_neq, 2*hLvl[0].m_nin);

      /* solve HCOD */
      hcod.activeSearch(m_output.x);
      activeSet = hcod.getOptimalActiveSet();
      // std::cerr<<"solution:\n"<<m_output.x.transpose()<<std::endl;
      std::cerr << "nrofasiterations "<<hcod.getNrASIterations()<<std::endl;

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
      // std::cout<<"solution: "<<m_output.x.transpose()<<std::endl;

      compute_slack(problemData, m_output);

      // std::cerr<<"recomputeSlack level 0"<<std::endl;
      // Eigen::VectorXd lhs = J[0]*m_output.x;
      // for (int c=0;c<hLvl[0].m_neq+hLvl[0].m_nin;c++)
      //   std::cerr<<b[0][c].distance(lhs[c])<<std::endl;

      // std::cerr<<"recomputeSlack level 1"<<std::endl;
      // Eigen::VectorXd lhs2 = J[1]*m_output.x;
      // for (int c=0;c<hLvl[1].m_neq+hLvl[1].m_nin;c++)
      //   std::cerr<<b[1][c].distance(lhs[c])<<std::endl;

      return m_output;

    }

    const std::vector<Eigen::MatrixXd>& SolverHQuadProgFast::getLagrangeMultipliers() {
        return lambda;
        // return hcod.getLagrangeMultipliers();
    }
    
    bool SolverHQuadProgFast::setMaximumIterations(unsigned int maxIter)
    {
      SolverHQPBase::setMaximumIterations(maxIter);
      return m_solver.setMaxIter(maxIter);
    }

    void SolverHQuadProgFast::compute_slack(const HQPData & problemData, 
                                            HQPOutput & problemOutput) {
      problemOutput.m_slack.resize(p);
      for (int l=0;l<p;l++) {
        // std::cerr<<"w["<<l<<"]:\n"<<hcod.getLagrangeMultipliers()[l]<<std::endl;
        problemOutput.m_slack[l] = hcod.getLagrangeMultipliers()[l].col(l).norm();
      }
    }

    const HQPOutput & SolverHQuadProgFast::solve_local(const HQPData & problemData,
                                                       const HQPOutput & previousOutput) {
    }
  }
}



