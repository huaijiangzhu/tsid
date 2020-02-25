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

#include <iostream>

#include <tsid/solvers/solver-HQP-factory.hxx>
#include <tsid/solvers/solver-HQP-eiquadprog.hpp>
#include <tsid/solvers/solver-HQP-eiquadprog-fast.hpp>
#include <tsid/solvers/solver-HQP-eiquadprog-rt.hpp>
#include <tsid/math/utils.hpp>
#include <tsid/math/constraint-equality.hpp>
#include <tsid/math/constraint-inequality.hpp>
#include <tsid/math/constraint-bound.hpp>
#include <tsid/utils/stop-watch.hpp>
#include <tsid/utils/statistics.hpp>

#define CHECK_LESS_THAN(A,B) BOOST_CHECK_MESSAGE(A<B, #A<<": "<<A<<">"<<B)
#define REQUIRE_FINITE(A) BOOST_REQUIRE_MESSAGE(isFinite(A), #A<<": "<<A)

#define PROFILE_EIQUADPROG "Eiquadprog"
#define PROFILE_EIQUADPROG_RT "Eiquadprog Real Time"
#define PROFILE_EIQUADPROG_FAST "Eiquadprog Fast"

int main()
{
  std::cout << "test_HCOD\n";
  using namespace tsid;
  using namespace math;
  using namespace solvers;

  int p = 2;
  int n = 5;
  int nin = n;
  int neq = 0;
  SolverHQPBase * solver_fast = SolverHQPFactory::createNewSolver(SOLVER_HQP_EIQUADPROG_FAST,
                                                                  "eiquadprog_fast");
  std::cout<<"unittest create new solver"<<std::endl;
  solver_fast->resize(n, neq, nin);

  // CREATE PROBLEM DATA
  HQPData HQPData(p);

  Matrix A1 = Matrix::Identity(n, n);
  Vector b1 = Vector::Zero(n);
  ConstraintEquality cost("c1", A1, b1);
  HQPData[1].push_back(make_pair<double, ConstraintBase*>(1.0, &cost));
  std::cout<<"objective data A:\n"<<A1<<std::endl;
  std::cout<<"objective data b:\n"<<b1<<std::endl;

  Vector x(n);

  Matrix A_in = Matrix::Identity(nin, n);
  Vector A_lb = Vector::Constant(nin,1.);
  Vector A_ub = Vector::Constant(nin,2.);
  ConstraintInequality in_constraint("in1", A_in, A_lb, A_ub);
  HQPData[0].push_back(make_pair<double, ConstraintBase*>(1.0, &in_constraint));
  std::cout<<"constraint data A_in:\n"<<A_in<<std::endl;
  std::cout<<"constraint data A_lb:\n"<<A_lb<<std::endl;
  std::cout<<"constraint data A_ub:\n"<<A_ub<<std::endl;

  Matrix A_eq = Matrix::Identity(neq, n);
  Vector b_eq = 1.5*Vector::Ones(neq);
  ConstraintEquality eq_constraint("eq1", A_eq, b_eq);
  HQPData[0].push_back(make_pair<double, ConstraintBase*>(1.0, &eq_constraint));
  std::cout<<"constraint data A:\n"<<A_eq<<std::endl;
  std::cout<<"constraint data b:\n"<<b_eq<<std::endl;

  const HQPOutput & output = solver_fast->solve(HQPData);
  std::cout << "solution: " << output.x.transpose() << std::endl;
  std::vector<Eigen::MatrixXd> lambda = static_cast<SolverHQuadProgFast*>(solver_fast)->getLagrangeMultipliers();
  std::vector<soth::cstref_vector_t> activeSet = static_cast<SolverHQuadProgFast*>(solver_fast)->getActiveSet();
  std::vector<Eigen::VectorXd> slack;
  slack.resize(p);
  int l = 0;
  for (const auto& lam : lambda)
  {
      std::cout<<"lagrange multipliers:\n"<<lam<<std::endl;
      slack[l] = lam.col(l);
      l += 1;
  }
  std::cout<<"slack:\n"<<output.m_slack<<std::endl;
  for (const auto& as : activeSet)
      for (const auto& c : as)
      std::cout<<"active set:\n"<<c.type<<std::endl;

  b_eq = 2.5*Vector::Ones(neq);
  eq_constraint = ConstraintEquality("eq1", A_eq, b_eq);
  HQPData[0].clear();
  HQPData[0].push_back(make_pair<double, ConstraintBase*>(1.0, &eq_constraint));
  std::cout<<"constraint data A:\n"<<A_eq<<std::endl;
  std::cout<<"constraint data b:\n"<<b_eq<<std::endl;

  const HQPOutput & output2 = static_cast<SolverHQuadProgFast*>(solver_fast)->resolve(HQPData);
  std::cout << "solution: " << output2.x.transpose() << std::endl;
  lambda = static_cast<SolverHQuadProgFast*>(solver_fast)->getLagrangeMultipliers();
  for (const auto& lam : lambda)
      std::cout<<"lagrange multipliers:\n"<<lam<<std::endl;

  activeSet = static_cast<SolverHQuadProgFast*>(solver_fast)->getActiveSet();

}
