//
// Copyright (c) 2018 CNRS
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

#ifndef __tsid_python_solver_quadprog_hpp__
#define __tsid_python_solver_quadprog_hpp__

#include <boost/python.hpp>
#include <string>
#include <eigenpy/eigenpy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "tsid/solvers/solver-HQP-eiquadprog.hpp"
#include "tsid/solvers/solver-HQP-eiquadprog-fast.hpp"
#include "tsid/solvers/solver-HQP-output.hpp"
#include "tsid/solvers/fwd.hpp"
#include "tsid/bindings/python/utils/container.hpp"

namespace tsid
{
  namespace python
  {    
    namespace bp = boost::python;
    
    template<typename Solver>
    struct SolverHQuadProgPythonVisitor
    : public boost::python::def_visitor< SolverHQuadProgPythonVisitor<Solver> >
    {
      template<class PyClass>     

      void visit(PyClass& cl) const
      {
        cl
        .def(bp::init<std::string>((bp::arg("name")), "Default Constructor with name"))
        
        .def("resize", &SolverHQuadProgPythonVisitor::resize, bp::args("n", "neq", "nin"))
        .add_property("ObjVal", &Solver::getObjectiveValue, "return obj value")
        .def("solve", &SolverHQuadProgPythonVisitor::solve, bp::args("HQPData"))
        .def("solve_local", &SolverHQuadProgPythonVisitor::solve_local, bp::args("HQPData","previousOutput"))
        .def("compute_slack", &SolverHQuadProgPythonVisitor::compute_slack, bp::args("HQPData","problemOutput"))
        .def("solve", &SolverHQuadProgPythonVisitor::solver_helper, bp::args("HQPData for Python"))
        .def("solve_local", &SolverHQuadProgPythonVisitor::solve_local_helper, bp::args("HQPData for Python","previousOutput"))
        .def("compute_slack", &SolverHQuadProgPythonVisitor::compute_slack_helper, bp::args("HQPData for Python","problemOutput"))
        .def("initializeSolver", &SolverHQuadProgPythonVisitor::initializeSolver_helper, bp::args("HQPData for Python"))
        .def("reinitializeSolver", &SolverHQuadProgPythonVisitor::reinitializeSolver, bp::args("HQPData for Python"))
        .def("setInitialActiveSet", &SolverHQuadProgPythonVisitor::setInitialActiveSet, bp::args("initialActiveSet"))
        ;
      }
       
      static void resize(Solver & self, unsigned int n, unsigned int neq, unsigned int nin){
          self.resize(n, neq, nin);
      }  
      static solvers::HQPOutput solve(Solver & self, const solvers::HQPData & problemData){
          solvers::HQPOutput output;
          output = self.solve(problemData);
          return output;
      }
      static solvers::HQPOutput solve_local(Solver & self, const solvers::HQPData & problemData,
                                      const solvers::HQPOutput & previousOutput){
          solvers::HQPOutput output;
          output = self.solve_local(problemData, previousOutput);
          return output;
      }
      static solvers::HQPOutput solve_local_helper(Solver & self, HQPDatas & HQPDatas,
                                                   const solvers::HQPOutput & previousOutput){
          solvers::HQPOutput output;
          solvers::HQPData data = HQPDatas.get();
          output = self.solve_local(data, previousOutput);
          return output;
      }
      static void compute_slack(Solver & self, solvers::HQPData & problemData,
                                solvers::HQPOutput & problemOutput){
          solvers::HQPOutput output;
          self.compute_slack(problemData, problemOutput);
      }
      static void compute_slack_helper(Solver & self, HQPDatas & HQPDatas,
                                       solvers::HQPOutput & problemOutput){
          solvers::HQPOutput output;
          solvers::HQPData data = HQPDatas.get();
          self.compute_slack(data, problemOutput);
      }
      static void initializeSolver(Solver & self, solvers::HQPData & problemData){
          self.initializeSolver(problemData);
      }
      static void initializeSolver_helper(Solver & self, HQPDatas & HQPDatas){
          solvers::HQPData data = HQPDatas.get();
          self.initializeSolver(data);
      }
      static void reinitializeSolver(Solver & self){
          self.reinitializeSolver();
      }
      static void setInitialActiveSet(Solver & self, const std::vector<soth::cstref_vector_t> & activeSetIn){
          self.setInitialActiveSet(activeSetIn);
      }
      static solvers::HQPOutput solver_helper(Solver & self, HQPDatas & HQPDatas){
          solvers::HQPOutput output;
          solvers::HQPData data = HQPDatas.get();
          output = self.solve(data);
          return output;
      }

      static void expose(const std::string & class_name)
      {
        std::string doc = "Solver EiQuadProg info.";
        bp::class_<Solver>(class_name.c_str(),
                          doc.c_str(),
                          bp::no_init)
        .def(SolverHQuadProgPythonVisitor<Solver>());       
      }
    };
  }
}


#endif // ifndef __tsid_python_solver_quadprog_hpp__
