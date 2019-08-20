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

#ifndef __tsid_python_HQPOutput_hpp__
#define __tsid_python_HQPOutput_hpp__

#include <boost/python.hpp>`
#include <string>
#include <eigenpy/eigenpy.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>

#include "tsid/solvers/solver-HQP-output.hpp"

namespace tsid
{
  namespace python
  {    
    namespace bp = boost::python;
    
    template<typename T>
    struct HQPOutputPythonVisitor
    : public boost::python::def_visitor< HQPOutputPythonVisitor<T> >
    {
      template<class PyClass>     

      void visit(PyClass& cl) const
      {
        cl
        .def(bp::init<>("Defulat Constructor"))
        .def(bp::init<int, int, int>((bp::args("nVars", "nEq", "nInCon"))))
        .add_property("x", &HQPOutputPythonVisitor::x)
        .add_property("status", &HQPOutputPythonVisitor::status)
        .add_property("activeSet", &HQPOutputPythonVisitor::activeSet)
        .add_property("A", &HQPOutputPythonVisitor::A)
        .add_property("b", &HQPOutputPythonVisitor::b)
        .add_property("H", &HQPOutputPythonVisitor::H)
        .add_property("g", &HQPOutputPythonVisitor::g)
        .add_property("K", &HQPOutputPythonVisitor::K)
        .add_property("Kinv", &HQPOutputPythonVisitor::Kinv)
        .add_property("delta", &HQPOutputPythonVisitor::delta)
        ;
      }
      static int status (const T & self) {return self.status;}
      static Eigen::VectorXd x (const T & self) {return self.x;}
      static Eigen::VectorXd activeSet (const T & self) {return self.activeSetPy;}
      static Eigen::MatrixXd A (const T & self) {return self.m_A;}
      static Eigen::MatrixXd H (const T & self) {return self.m_H;}
      static Eigen::MatrixXd K (const T & self) {return self.m_K;}
      static Eigen::MatrixXd Kinv (const T & self) {return self.m_Kinv;}
      static Eigen::VectorXd b (const T & self) {return self.m_b;}
      static Eigen::VectorXd g (const T & self) {return self.m_g;}
      static Eigen::VectorXd delta (const T & self) {return self.m_delta;}
      static void expose(const std::string & class_name)
      {
        std::string doc = "HQPOutput info.";
        bp::class_<T>(class_name.c_str(),
                          doc.c_str(),
                          bp::no_init)
        .def(HQPOutputPythonVisitor<T>());       
      }
    };
  }
}


#endif // ifndef __tsid_python_HQPOutput_hpp__
