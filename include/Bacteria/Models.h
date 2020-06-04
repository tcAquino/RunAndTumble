//
//  Models.h
//  RunAndTumble
//
//  Created by Tomás Aquino on 22/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef Models_Bacteria_h
#define Models_Bacteria_h

#include <exception>
#include <string>
#include <fstream>
#include <iostream>
#include <vector>
#include "general/Constants.h"
#include "Field/ScalarField.h"
#include "Grid/Grid.h"
#include "ODE/CrankNicolson.h"
#include "Stochastic/CTRW/Boundary.h"
#include "Stochastic/CTRW/CTRW.h"
#include "Stochastic/CTRW/Grid.h"
#include "Stochastic/CTRW/JumpGenerator.h"
#include "Stochastic/CTRW/Particle.h"
#include "Stochastic/CTRW/PTRW.h"
#include "Stochastic/CTRW/Reaction.h"
#include "Stochastic/CTRW/State.h"
#include "Stochastic/CTRW/StateSwitcher.h"
#include "Stochastic/CTRW/Transitions_State.h"

namespace bacteria
{
  // Discretization grid
  using Grid = grid::RegularGrid<useful::Empty>;

  // Eulerian fields
  using Concentration = field::ScalarField<Grid, double>;
  using Solver = ode::CrankNicolson_Diffusion;
  using BC_Container = Solver::BoundaryCondition_Container;
  using Gradient = field::Gradient<
    Concentration, Grid, BC_Container>;

  // Eulerian-Lagrangian coupling
  using Concentration_Particle = ctrw::GridParticles_Quantity<
    ctrw::GridParticles, Concentration>;
  using Gradient_Particle = ctrw::GridParticles_Quantity<
    ctrw::GridParticles, Gradient>;
  using Reaction =
    ctrw::Reaction_GridConcentrationDecay_ParticleNr_GridConcentration;

  // Particle tracking
  using State = ctrw::State_RunTumble_PTRW<double, std::size_t>;
  using Particle = ctrw::Particle<State>;
  using JumpGenerator = ctrw::JumpGenerator_JumpAngle_2d;
  using CTRW = ctrw::CTRW<State>;
  
  namespace model_test
  {
    char model_name[64] = "test";
    
    // Spatial dimension
    const std::size_t dim = 2;
    
    // Particle tracking
    using OrientationGenerator = ctrw::OrientationGenerator_Gradient_1d<
      Concentration_Particle, Gradient_Particle>;
    using OrientationGenerator_Wall =
      ctrw::OrientationGenerator_Flip;
    using StateSwitcher =
      ctrw::StateSwitcher_ConstantRate;
    using Boundary =
      boundary::ClosestVoidCenter_Tumble<Grid, dim>;
    
    // Eulerian fields
    const bool chemoattractant = false;
    
    struct Parameters_Discretization
    {
      double time_step;
      std::vector<std::size_t> nr_grid_cells;
      
      void get(std::string const& filename)
      {
        std::ifstream file(filename);
        if (!file.is_open())
          throw std::runtime_error{ "Could not open " +
            filename + " for reading" };
        file >> time_step;
        for (std::size_t dd = 0; dd < dim; ++dd )
          file >> nr_grid_cells[dd];
        if (file.fail())
          throw std::invalid_argument{ "Inappropriate input file " +
            filename };
      }
      
      static void print_parameter_list()
      {
        std::cout << "Discretization parameters:\n"
                  << "\ttime_step :     Discretization time step\n"
                  << "\tnr_grid_cells : Number of grid cells along each dimension (dim components)\n";
      }
    };
    
    struct Parameters_EulerianFields
    {
      double nutrient_diffusion;
      
      void get(std::string const& filename)
      {
        std::ifstream file(filename);
        if (!file.is_open())
          throw std::runtime_error{ "Could not open " +
            filename + " for reading" };
        file >> nutrient_diffusion;
        if (file.fail())
          throw std::invalid_argument{ "Inappropriate input file " +
            filename };
      }
      
      static void print_parameter_list()
      {
        printf("EulerianFields parameters:\n"
               "\tnutrient_diffusion : Nutrient diffusion coefficient\n");
      }
    };
    
    struct Parameters_Bacteria
    {
      double run_velocity;
      double rate_run_to_tumble;
      double rate_tumble_to_run;
      double rate_wall_tumble_to_run;
      double angle_variance_factor;
      double reaction_rate_volume_per_mass;
      double preferred_concentration;
      
      double angle_variance;
      
      void get(std::string const& filename)
      {
        std::ifstream file(filename);
        if (!file.is_open())
          throw std::runtime_error{ "Could not open " +
            filename + " for reading" };
        file >> run_velocity
             >> rate_run_to_tumble >> rate_tumble_to_run
             >> rate_wall_tumble_to_run >> angle_variance_factor
             >> reaction_rate_volume_per_mass
             >> preferred_concentration;
        if (file.fail())
          throw std::invalid_argument{ "Inappropriate input file " +
            filename };
        angle_variance = 2.*constants::pi*angle_variance_factor;
      }
      
      static void print_parameter_list()
      {
        std::cout << "Bacteria parameters:\n"
                  << "\trun_velocity :                  Bacteria velocity during runs\n"
                  << "\trate_run_to_tumble :            Rate of switching from run to tumble\n"
                  << "\trate_tumble_to_run :            Rate of switching from tumble to run\n"
                  << "\trate_wall_tumble_to_run :       Rate of switching from wall-tumble to run\n"
                  << "\tangle_variance_factor :         Variance of new angle relative to gradient\n"
                  << "\t                                direction during tumble, as fraction of 2pi\n"
                  << "\treaction_rate_volume_per_mass : Rate of nutrient consumption [volume/mass/time]\n"
                  << "\tpreferred_concentration :       Nutrient concentration preferred by bacteria\n";
      }
    };
  }
}

#endif /* Models_Bacteria_h */
