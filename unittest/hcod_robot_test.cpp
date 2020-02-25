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

#include <tsid/contacts/contact-6d.hpp>
#include <tsid/contacts/contact-point.hpp>
#include <tsid/formulations/inverse-dynamics-formulation-acc-force.hpp>
#include <tsid/tasks/task-com-equality.hpp>
#include <tsid/tasks/task-se3-equality.hpp>
#include <tsid/tasks/task-joint-posture.hpp>
#include <tsid/tasks/task-joint-bounds.hpp>
#include <tsid/trajectories/trajectory-euclidian.hpp>
#include <tsid/solvers/solver-HQP-factory.hxx>
#include <tsid/solvers/utils.hpp>
#include <tsid/utils/stop-watch.hpp>
#include <tsid/utils/statistics.hpp>
#include <tsid/math/utils.hpp>

#include <pinocchio/algorithm/joint-configuration.hpp> // integrate
#include <pinocchio/parsers/srdf.hpp>

using namespace tsid;
using namespace tsid::trajectories;
using namespace tsid::math;
using namespace tsid::contacts;
using namespace tsid::tasks;
using namespace tsid::solvers;
using namespace tsid::robots;
using namespace std;

const string romeo_model_path = TSID_SOURCE_DIR"/models/romeo";
const string quadruped_model_path = TSID_SOURCE_DIR"/models/quadruped";

#ifndef NDEBUG
const int max_it = 10;
#else
const int max_it = 1000;
#endif

int main()
{
  std::cout << "test_HCOD\n";
  vector<string> package_dirs;
  package_dirs.push_back(romeo_model_path);
  const string urdfFileName = package_dirs[0] + "/urdf/romeo.urdf";
  RobotWrapper * robot = new RobotWrapper(urdfFileName, package_dirs, pinocchio::JointModelFreeFlyer());
  const unsigned int nv = robot->nv();
  double t = 0.;
  Vector q = neutral(robot->model());
  q(2) += 0.84;
  std::cout << "q: " << q.transpose() << std::endl;
  Vector v = Vector::Zero(nv);

  std::cout<<"create invdyn data"<<std::endl;
  InverseDynamicsFormulationAccForce * tsid = new InverseDynamicsFormulationAccForce("tsid", *robot);
  std::cout<<"compute invdyn data"<<std::endl;
  tsid->computeProblemData(t, q, v);
  std::cout<<"get invdyn data"<<std::endl;
  const pinocchio::Data & data = tsid->data();
  std::cout<<"invdyn data end"<<std::endl;
}
