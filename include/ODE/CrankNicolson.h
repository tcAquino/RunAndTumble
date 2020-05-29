//
//  CrankNicolson.h
//  CrankNicolson
//
//  Created by Tomás Aquino on 20/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef CrankNicolson_h
#define CrankNicolson_h

#include <exception>
#include <vector>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include "BoundaryConditions.h"
#include "general/Operations.h"

namespace ode
{
  // Crank-Nicolson (CN) solver for diffusion problem
  // with constant coefficients in regular grid
  class CrankNicolson_Diffusion
  {
  public:
    // Eigen objects
    using Matrix = Eigen::SparseMatrix<double>;
    using Vector = Eigen::VectorXd;
    using Triplet = Eigen::Triplet<double, std::size_t>;
    
    // Robin boundary conditions can be set in each wall of each cell
    // Periodic BCs should be handled directly by grid
    using BoundaryCondition_1d = BoundaryCondition_Robin_1d;
    using BoundaryConditions = BoundaryConditions_Cell<BoundaryCondition_Robin_1d>;
    using BoundaryCondition_Container = std::vector<BoundaryConditions>;
    
    // Set up a new solver
    // grid : The full underlying grid
    // points : The indexes of the grid points to solve at
    // bcs : vector with the type of BC in each cell, and the Robin coefficients where needed
    // time_step : discretization time step
    // time : initial time
    // diff : diffusion coefficient along each dimension
    // cell_size : discretization length along each dimension
    template <typename Grid>
    CrankNicolson_Diffusion
    (Grid const& grid, std::vector<std::size_t> const& points,
     BoundaryCondition_Container const& bcs,
     double time_step, double time,
     std::vector<double> diff, std::vector<double> cell_size)
    : time_current{ time }
    , timestep{ time_step }
    , cell_size{ cell_size }
    , points{ points }
    , bcs{ bcs }
    , lhs(grid.size(), grid.size())
    , rhs(grid.size(), grid.size())
    , rhs_ind{ Vector::Zero(grid.size()) }
    , rhs_val{ Vector::Zero(grid.size()) }
    {
      for (std::size_t dd = 0; dd < grid.dim(); ++dd)
        coeff.push_back(diff[dd]*time_step/(2.*cell_size[dd]*cell_size[dd]));
      
      compute_matrices(grid);
      solver.compute(lhs);
    }
    
    // Same discretization length along each dimension
    template <typename Grid>
    CrankNicolson_Diffusion
    (Grid const& grid, std::vector<std::size_t> const& points,
     BoundaryCondition_Container const& bcs,
     double time_step, double time,
     double diff, std::vector<double> const& cell_size)
    : CrankNicolson_Diffusion{
      grid, points, bcs, time_step, time,
      std::vector<double>(grid.dim(), diff), cell_size }
    {}
    
    // Same diffusion coefficient along each dimension
    template <typename Grid>
    CrankNicolson_Diffusion
    (Grid const& grid, std::vector<std::size_t> const& points,
     BoundaryCondition_Container const& bcs,
     double time_step, double time,
     std::vector<double> const& diff, double cell_size)
    : CrankNicolson_Diffusion{
      grid, points, bcs, time_step, time,
      diff, std::vector<double>(grid.dim(), cell_size) }
    {}
    
    // Same discretization length and diffusion coefficient along each dimension
    template <typename Grid>
    CrankNicolson_Diffusion
    (Grid const& grid, std::vector<std::size_t> const& points,
     BoundaryCondition_Container const& bcs,
     double time_step, double time,
     double diff, double cell_size)
    : CrankNicolson_Diffusion{
      grid, points, bcs, time_step, time,
      std::vector<double>(grid.dim(), diff), std::vector<double>(grid.dim(), cell_size) }
    {}
    
    double time()
    { return time_current; }
    
    double time_step()
    { return timestep; }
    
    // Take one time step
    // values : initial values on the grid, takes the new values
    template <typename Container>
    void step(Container& values)
    {
      for (std::size_t idx = 0; idx < points.size(); ++idx)
        rhs_val[points[idx]] = values[points[idx]];
      
      Vector rhs_new = rhs*rhs_val + rhs_ind;
      rhs_val = solver.solveWithGuess(rhs_new, rhs_val);
      
      for (std::size_t idx = 0; idx < points.size(); ++idx)
        values[points[idx]] = rhs_val[points[idx]];
    }
    
    // Take time steps until time >= time_max
    // values : initial values on the grid, takes the new values
    template <typename Container>
    void evolve(double time_max, Container& values)
    {
      for (std::size_t idx = 0; idx < points.size(); ++idx)
        rhs_val[points[idx]] = values[points[idx]];
      
      while (time_current < time_max)
      {
        Vector rhs_new = rhs*rhs_val + rhs_ind;
        rhs_val = solver.solveWithGuess(rhs_new, rhs_val);
        time_current += timestep;
      }
      
      for (std::size_t idx = 0; idx < points.size(); ++idx)
        values[points[idx]] = rhs_val[points[idx]];
    }
    
  private:
    // Prepare the Crank-Nicolson system matrices
    template <typename Grid>
    void compute_matrices(Grid const& grid)
    {
      std::vector<Triplet> triplets_lhs;
      triplets_lhs.reserve((1+3*grid.dim())*points.size());
      
      std::vector<Triplet> triplets_rhs;
      triplets_rhs.reserve((1+3*grid.dim())*points.size());
      
      // Every grid point takes a 1 in the diagonals
      // For the points not to be solved, this is the only contribution,
      // so they don't change
      for (std::size_t idx = 0; idx < grid.size(); ++idx)
      {
        triplets_lhs.push_back({ idx, idx, 1. });
        triplets_rhs.push_back({ idx, idx, 1. });
      }
      
      // Set up coefficients for points to solve
      for (auto idx_current : points)
      {
        // Along each dimension, for each point
        for (std::size_t dd = 0; dd < grid.dim(); ++dd)
        {
          switch (bcs[idx_current].cell_type[dd])
          {
            // No walls : use nearest neighbors
            case 0:
            {
              std::size_t idx_plus = grid.neighbor_cell_up(idx_current, dd);
              std::size_t idx_minus = grid.neighbor_cell_down(idx_current, dd);

              triplets_lhs.push_back({ idx_current, idx_current, 2.*coeff[dd] });
              triplets_lhs.push_back({ idx_current, idx_plus, -coeff[dd] });
              triplets_lhs.push_back({ idx_current, idx_minus, -coeff[dd] });
              triplets_rhs.push_back({ idx_current, idx_current, -2.*coeff[dd] });
              triplets_rhs.push_back({ idx_current, idx_plus, coeff[dd] });
              triplets_rhs.push_back({ idx_current, idx_minus, coeff[dd] });
              break;
            }
            // Wall on left : use BC on left and nearest neighbor on right
            case 1:
            {
              std::size_t idx_other = grid.neighbor_cell_up(idx_current, dd);

              double der_coeff = bcs[idx_current](dd, 0).derivative_coeff;
              double val_coeff = bcs[idx_current](dd, 0).value_coeff;
              double ind_coeff = bcs[idx_current](dd, 0).ind_coeff;

              double current_size = cell_size[dd];

              double aux = 4./(3.*(8.*der_coeff-3.*current_size*val_coeff));
              double coeff_current = -coeff[dd]*aux*(6.*der_coeff+3.*current_size*val_coeff);
              double coeff_other = coeff[dd]*aux*(6.*der_coeff-3.*current_size*val_coeff);
              double coeff_ind = coeff[dd]*aux*3.*current_size*ind_coeff;

              triplets_lhs.push_back({ idx_current, idx_current, -coeff_current });
              triplets_lhs.push_back({ idx_current, idx_other, -coeff_other });
              triplets_rhs.push_back({ idx_current, idx_current, coeff_current });
              triplets_rhs.push_back({ idx_current, idx_other, coeff_other });

              rhs_ind[idx_current] += coeff_ind;
              break;
            }
            // Wall on right : use BC on right and nearest neighbor on left
            case 2:
            {
              std::size_t idx_other = grid.neighbor_cell_down(idx_current, dd);

              double der_coeff = bcs[idx_current](dd, 0).derivative_coeff;
              double val_coeff = bcs[idx_current](dd, 0).value_coeff;
              double ind_coeff = bcs[idx_current](dd, 0).ind_coeff;

              double current_size = -cell_size[dd];

              double aux = 4./(3.*(8.*der_coeff-3.*current_size*val_coeff));
              double coeff_current = -coeff[dd]*aux*(6.*der_coeff+3.*current_size*val_coeff);
              double coeff_other = coeff[dd]*aux*(6.*der_coeff-3.*current_size*val_coeff);
              double coeff_ind = coeff[dd]*aux*3.*current_size*ind_coeff;

              triplets_lhs.push_back({ idx_current, idx_current, -coeff_current });
              triplets_lhs.push_back({ idx_current, idx_other, -coeff_other });
              triplets_rhs.push_back({ idx_current, idx_current, coeff_current });
              triplets_rhs.push_back({ idx_current, idx_other, coeff_other });

              rhs_ind[idx_current] += coeff_ind;
              break;
            }
            // Wall on both sides : use BCs only (first order in space instead of second here)
            case 3:
            {
              double der_coeff_left = bcs[idx_current](dd, 0).derivative_coeff;
              double val_coeff_left = bcs[idx_current](dd, 0).value_coeff;
              double ind_coeff_left = bcs[idx_current](dd, 0).ind_coeff;
              double der_coeff_right = bcs[idx_current](dd, 1).derivative_coeff;
              double val_coeff_right = bcs[idx_current](dd, 1).value_coeff;
              double ind_coeff_right = bcs[idx_current](dd, 1).ind_coeff;

              double aux_left = 1./(2.*der_coeff_left-val_coeff_left*cell_size[dd]);
              double aux_right = 1./(2.*der_coeff_right+val_coeff_right*cell_size[dd]);

              double coeff_current = coeff[dd]*cell_size[dd]*
                (der_coeff_left*aux_left - der_coeff_right*aux_right);
              double coeff_ind = coeff[dd]*cell_size[dd]*
                (ind_coeff_left*aux_left - ind_coeff_right*aux_right);

              triplets_lhs.push_back({ idx_current, idx_current, -coeff_current });
              triplets_rhs.push_back({ idx_current, idx_current, coeff_current });

              rhs_ind[idx_current] += coeff_ind;
              break;
            }
            default:
            {
              throw std::runtime_error{ "Undefined cell type." };
            }
          }
        }
      }
      
      // Add coefficients to CN matrices
      lhs.setFromTriplets(triplets_lhs.cbegin(), triplets_lhs.cend());
      rhs.setFromTriplets(triplets_rhs.cbegin(), triplets_rhs.cend());
    }
    
    double time_current;
    double timestep;
    std::vector<double> cell_size;
    std::vector<double> coeff;
    std::vector<std::size_t> const& points;
    BoundaryCondition_Container const& bcs;
    Matrix lhs;
    Matrix rhs;
    Vector rhs_ind;
    Vector rhs_val;
    Eigen::BiCGSTAB<Matrix, Eigen::IncompleteLUT<double>> solver;
  };
}


#endif /* CrankNicolson_h */
