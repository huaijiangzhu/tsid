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

// #include <boost/test/unit_test.hpp>
// #include <boost/utility/binary.hpp>

#include <tsid/solvers/solver-HQP-factory.hxx>
#include <tsid/solvers/solver-HQP-eiquadprog.hpp>
#include <tsid/solvers/solver-HQP-eiquadprog-rt.hpp>
#include <tsid/math/utils.hpp>
#include <tsid/math/constraint-equality.hpp>
#include <tsid/math/constraint-inequality.hpp>
#include <tsid/math/constraint-bound.hpp>
#include <tsid/utils/stop-watch.hpp>
#include <tsid/utils/statistics.hpp>

#define CHECK_LESS_THAN(A,B) BOOST_CHECK_MESSAGE(A<B, #A<<": "<<A<<">"<<B)
#define REQUIRE_FINITE(A) BOOST_REQUIRE_MESSAGE(isFinite(A), #A<<": "<<A)

// BOOST_AUTO_TEST_SUITE ( BOOST_TEST_MODULE )

#define PROFILE_EIQUADPROG "Eiquadprog"
#define PROFILE_EIQUADPROG_RT "Eiquadprog Real Time"
#define PROFILE_EIQUADPROG_FAST "Eiquadprog Fast"

//BOOST_AUTO_TEST_CASE ( test_eiquadprog_classic_vs_rt_vs_fast)
int main()
{
  std::cout << "test_HCOD\n";
  using namespace tsid;
  using namespace math;
  using namespace solvers;

  int n = 5;
  int nin = 0;
  int neq = n;
  SolverHQPBase * solver_fast = SolverHQPFactory::createNewSolver(SOLVER_HQP_EIQUADPROG_FAST,
                                                                  "eiquadprog_fast");
  std::cout<<"unittest create new solver"<<std::endl;
  solver_fast->resize(n, neq, nin);

  // CREATE PROBLEM DATA
  HQPData HQPData(2);

  Matrix A1 = Matrix::Identity(n, n);
  Vector b1 = Vector::Zero(n);
  ConstraintEquality cost("c1", A1, b1);
  HQPData[1].push_back(make_pair<double, ConstraintBase*>(1.0, &cost));
  std::cout<<"objective data A:\n"<<A1<<std::endl;
  std::cout<<"objective data b:\n"<<b1<<std::endl;

  Vector x(n);

  Matrix A_in = Matrix::Identity(nin, n);
  Vector A_lb = Vector::Random(nin);
  Vector A_ub = Vector::Random(nin);
  Vector constrVal = A_in*x;
  for(unsigned int i=0; i<nin; i++)
  {
      if(constrVal[i]>A_ub[i])
      {
//        std::cout<<"Inequality constraint "<<i<<" active at first iteration. UB="<<A_ub[i]<<", value="<<constrVal[i]<<endl;
        A_ub[i] = constrVal[i];
      }
      if(constrVal[i]<A_lb[i])
      {
//        std::cout<<"Inequality constraint "<<i<<" active at first iteration. LB="<<A_lb[i]<<", value="<<constrVal[i]<<endl;
        A_lb[i] = constrVal[i];
      }
  }
  ConstraintInequality in_constraint("in1", A_in, A_lb, A_ub);
  HQPData[0].push_back(make_pair<double, ConstraintBase*>(1.0, &in_constraint));

  Matrix A_eq = Matrix::Identity(neq, n);
  Vector b_eq = 1.5*Vector::Ones(neq);
  ConstraintEquality eq_constraint("eq1", A_eq, b_eq);
  HQPData[0].push_back(make_pair<double, ConstraintBase*>(1.0, &eq_constraint));
  std::cout<<"constraint data A:\n"<<A_eq<<std::endl;
  std::cout<<"constraint data b:\n"<<b_eq<<std::endl;

  const HQPOutput & output = solver_fast->solve(HQPData);
  std::cout << "solution: " << output.x.transpose() << std::endl;
  std::vector<Matrix> lambda = solver_fast->getLagrangeMultipliers();
  for (const auto& lam : lambda)
      std::cout<<"lagrange multipliers:\n"<<lam<<std::endl;

  b_eq = 2.5*Vector::Ones(neq);
  eq_constraint = ConstraintEquality("eq1", A_eq, b_eq);
  HQPData[0].clear();
  HQPData[0].push_back(make_pair<double, ConstraintBase*>(1.0, &eq_constraint));
  std::cout<<"constraint data A:\n"<<A_eq<<std::endl;
  std::cout<<"constraint data b:\n"<<b_eq<<std::endl;

  const HQPOutput & output2 = solver_fast->resolve(HQPData);
  std::cout << "solution: " << output2.x.transpose() << std::endl;
  lambda = solver_fast->getLagrangeMultipliers();
  for (const auto& lam : lambda)
      std::cout<<"lagrange multipliers:\n"<<lam<<std::endl;

  // BOOST_REQUIRE_MESSAGE(true,
  //                       " Status FAST "+SolverHQPBase::HQP_status_string[output.status]);

  // std::cout<<"\n### TEST FINISHED ###\n";
  // getProfiler().report_all(3, std::cout);
  // getStatistics().report_all(1, std::cout);
}



// BOOST_AUTO_TEST_SUITE_END ()
