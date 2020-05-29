//
//  BoundaryConditions.h
//  CrankNicolson
//
//  Created by Tomás Aquino on 21/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef BoundaryConditions_ODE_h
#define BoundaryConditions_ODE_h

#include <vector>
#include "Grid/Grid.h"

namespace ode
{
  struct BoundaryCondition_Robin_1d
  {
    double derivative_coeff;
    double value_coeff;
    double ind_coeff;
  };

  template <typename BoundaryCondition_1d>
  struct BoundaryConditions_Cell
  {
    // Boundary type of cell along each dimension
    std::vector<int> cell_type;
    // BCs for each dimension and for each boundary
    std::vector<std::vector<BoundaryCondition_1d>> bcs;
    
    // Return the cond BC along dimension dd
    BoundaryCondition_1d operator()(std::size_t dd, std::size_t cond) const
    { return bcs[dd][cond]; }
  };
  
  // Each boundary between a void and a solid is reflecting
  template <typename Grid>
  auto make_reflecting_bcs
  (Grid const& grid,
   std::vector<std::size_t> const& void_idx,
   std::vector<std::size_t> const& solid_idx)
  {
    // BC vector of the size of the full grid
    // Reflecting boundaries implemented as a special case of Robin BC
    std::vector<BoundaryConditions_Cell<BoundaryCondition_Robin_1d>> cell_bcs(grid.size());
    
    // Add info for the voids
    for (auto idx : void_idx)
    {
      auto& bcs = cell_bcs[idx];
      auto& cell_type = bcs.cell_type;
      // Boundary type along each dimension
      cell_type = grid::void_type(grid, solid_idx, idx);
      
      for (std::size_t dd = 0; dd < grid.dim(); ++dd)
      {
        // No boundaries
        bcs.bcs.push_back({});
        auto& current_bcs = bcs.bcs.back();
        // Left or right boundary reflecting condition
        if (cell_type[dd] != 0)
          current_bcs.push_back({ 1., 0., 0. });
        // Remaining boundary reflecting condition when both are solid
        if (cell_type[dd] == 3)
          current_bcs.push_back({ 1., 0., 0. });
      }
    }
    
    return cell_bcs;
  }
}

#endif /* BoundaryConditions_ODE_h */
