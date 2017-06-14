//
// Copyright (c) 2017 CNRS
//
// This file is part of PinInvDyn
// PinInvDyn is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
// PinInvDyn is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// PinInvDyn If not, see
// <http://www.gnu.org/licenses/>.
//

#ifndef __invdyn_solvers_hqp_utils_hpp__
#define __invdyn_solvers_hqp_utils_hpp__

#include "pininvdyn/solvers/fwd.hpp"

#include <string>


namespace pininvdyn
{
  namespace solvers
  {
    
    std::string HQPDataToString(const HQPData & data, bool printMatrices=false);
  }
  
}

#endif // ifndef __invdyn_solvers_hqp_utils_hpp__