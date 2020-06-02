//
//  main.cpp
//  RunAndTumble
//
//  Created by Tomás Aquino on 02/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#include <cmath>
#include <exception>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "Bacteria/Domains.h"
#include "Bacteria/InitialConditions.h"
#include "Bacteria/Models.h"
#include "Bacteria/Output.h"
#include "general/Constants.h"
#include "general/Ranges.h"
#include "general/useful.h"

int main(int argc, const char * argv[])
{
  using namespace bacteria;
	using namespace model_test;
  using namespace domain_test_ring;
  using namespace initial_condition_particles_line_ring;
  using namespace initial_condition_fields_constant;
  using namespace output_standard;
  
  if (argc == 1)
  {
    Parameters_EulerianFields::print_parameter_list();
    Parameters_Bacteria::print_parameter_list();
    Parameters_Domain::print_parameter_list();
    Parameters_Discretization::print_parameter_list();
    Parameters_ICBacteria::print_parameter_list();
    Parameters_ICEulerianFields::print_parameter_list();
    Parameters_Output::print_parameter_list();
    
    return 0;
  }
  
  std::cout << "Getting parameters...\n";
  if (argc < 3 || argc > 5)
    throw std::invalid_argument{ "Inappropriate parameters" };
  
  std::size_t arg = 1;
  std::string parameter_set_name = argv[arg++];
  
  std::size_t run_nr = argc <= 2
  ? 0 : strtoul(argv[arg++], NULL, 0);
  
  std::string directory_params = argc <= 3
  ? "../parameters" : argv[arg++];
  std::string directory_data = argc <= 4
  ? "../data" : argv[arg++];
  
  Parameters_EulerianFields parameters_EulerianFields;
  parameters_EulerianFields.get(directory_params +
    "/parameters_EulerianFields_" +
    parameter_set_name + ".dat");
  
  Parameters_Bacteria parameters_Bacteria;
  parameters_Bacteria.get(directory_params +
    "/parameters_Bacteria_" +
    parameter_set_name + ".dat");
  
  Parameters_Domain parameters_Domain;
  parameters_Domain.get(directory_params +
    "/parameters_Output_" +
    parameter_set_name + ".dat");
  
  Parameters_Discretization parameters_Discretization;
  parameters_Discretization.get(directory_params +
    "/parameters_Discretization_" + parameter_set_name + ".dat");
  
  Parameters_ICBacteria parameters_ICBacteria;
  parameters_ICBacteria.get(directory_params +
    "/parameters_ICBacteria_" +
    parameter_set_name + ".dat");
  
  Parameters_ICEulerianFields parameters_ICEulerianFields;
  parameters_ICEulerianFields.get(directory_params +
    "/parameters_ICEulerianFields_" + parameter_set_name + ".dat");
  
  Parameters_Output parameters_Output;
  parameters_Output.get(directory_params +
    "/parameters_Output_" + parameter_set_name + ".dat");
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up grid...\n";
  auto domain = make_domain(parameters_Domain);
  std::vector<double> discretization_length = operation::div(domain_dimensions,
    parameters_Discretization.nr_grid_cells);
  std::vector<std::size_t> discretization_nr_padded =
    operation::plus_scalar(parameters_Discretization.nr_grid_cells, 2);
  std::vector<double> domain_corner = operation::minus(domain.box.center,
    operation::div_scalar(domain_dimensions, 2.));
  std::vector<double> grid_corner = operation::minus(domain_corner,
    discretization_length);
  Grid grid{ discretization_nr_padded, discretization_length, grid_corner };
  grid::KDTree_Grid<Grid, dim> kdtree_grid{ grid };
  Grid_void_solid void_solid{ domain, kdtree_grid };
  grid::KDTree_Grid_Mask<Grid, dim> kdtree_void{ grid, void_solid.voids() };
  grid::KDTree_Grid_Mask<Grid, dim> kdtree_solid{ grid, void_solid.solids() };
  auto bcs{
    ode::make_reflecting_bcs(grid,
      void_solid.voids(), void_solid.solids()) };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up Eulerian fields...\n";
  std::vector<Concentration> concentration_fields{
    { grid,
      NutrientValue{ parameters_ICEulerianFields },
      void_solid.voids() } };
  if (chemoattractant)
    concentration_fields.push_back({
      grid,
      ChemoAttractantValue{ parameters_ICEulerianFields },
      void_solid.voids() });
  
  std::vector<Gradient> gradient_fields{
    { concentration_fields[0], grid, bcs } };
  if (chemoattractant)
    gradient_fields.push_back({
      concentration_fields[1], grid, bcs });
  
  std::cout << "\tSetting up Eulerian solvers...\n";
  std::vector<std::unique_ptr<Solver>> solvers;
  solvers.push_back(std::make_unique<Solver>(
    grid, void_solid.voids(), bcs,
    parameters_Discretization.time_step, 0.,
    parameters_EulerianFields.nutrient_diffusion,
    discretization_length));
  if (chemoattractant)
    solvers.emplace_back(std::make_unique<Solver>(
      grid, void_solid.voids(), bcs,
      parameters_Discretization.time_step, 0.,
      parameters_EulerianFields.nutrient_diffusion,
      discretization_length));
  std::cout << "\t\tDone!\n";
  std::cout << "\tDone!\n";

  std::cout << "Setting up particle tracking...\n";
  CTRW ctrw{ ParticleMaker<Particle>{ parameters_ICBacteria }(domain) };
  
  ctrw::GridParticles particle_grid{ grid, ctrw.particles() };
  std::vector<Concentration_Particle> concentrations_particle{
    { particle_grid, concentration_fields[0],
      Concentration_Particle::Grid_update{} } };
  if (chemoattractant)
    concentrations_particle.push_back({
      particle_grid, concentration_fields[1],
      Concentration_Particle::Grid_update{} });
  
  std::vector<Gradient_Particle> gradients_particle{
    { particle_grid, gradient_fields[0],
      Gradient_Particle::Grid_update{} } };
  if (chemoattractant)
    gradients_particle.push_back({
      particle_grid, gradient_fields[1],
      Gradient_Particle::Grid_update{} });

  double particle_mass = parameters_ICBacteria.total_mass
    /parameters_ICBacteria.nr_particles;
  double reaction_rate_volume_per_particle = particle_mass*
  parameters_Bacteria.reaction_rate_volume_per_mass;
  Reaction reaction{
    reaction_rate_volume_per_particle,
    parameters_Discretization.time_step };
  
  ctrw::Transitions_RunAndTumble_PTRW transitions{
    JumpGenerator{ parameters_Bacteria.run_velocity*
      parameters_Discretization.time_step },
    OrientationGenerator{ concentrations_particle[0],
      gradients_particle[0],
      parameters_Bacteria.angle_variance,
      parameters_Bacteria.preferred_concentration },
    OrientationGenerator_Wall{},
    StateSwitcher{ parameters_Discretization.time_step,
      parameters_Bacteria.rate_run_to_tumble,
      parameters_Bacteria.rate_tumble_to_run,
      parameters_Bacteria.rate_wall_tumble_to_run },
    Boundary{ grid, kdtree_void, kdtree_solid } };
  ctrw::PTRW ptrw{
    ctrw, transitions,
    parameters_Discretization.time_step, 0. };
  std::cout << "\tDone!\n";
  
  std::cout << "Setting up output...\n";
  std::vector<double> measure_times;
  switch (parameters_Output.measure_spacing)
  {
    case 0:
      measure_times = range::linspace(parameters_Output.time_min,
        parameters_Output.time_max, parameters_Output.nr_measures);
      break;
    case 1:
      if (parameters_Output.time_min <= 0.)
        throw std::invalid_argument{ "Nonpositive logspaced measure points" };
      measure_times = range::logspace(parameters_Output.time_min,
        parameters_Output.time_max, parameters_Output.nr_measures);
      break;
    default:
      throw std::runtime_error{ "Undefined measure spacing." };
  }
  
  std::vector<std::unique_ptr<Measurer_Particle>> measurer_particle;
  if (parameters_Output.output_particle)
    measurer_particle.push_back(std::make_unique<Measurer_Particle>(
        directory_data + "/Data_particle_state_" + model_name + domain_name +
        initial_condition_particles_name + initial_condition_fields_name +
        parameter_set_name + std::to_string(run_nr) + ".dat"));
  
  std::vector<std::unique_ptr<Measurer_Field>> measurer_field;
  if (parameters_Output.output_chemoattractant
      && chemoattractant)
    measurer_field.push_back(std::make_unique<Measurer_Field>(
        directory_data + "/Data_nutrient_" + model_name + domain_name +
        initial_condition_particles_name + initial_condition_fields_name +
        parameter_set_name + std::to_string(run_nr) + ".dat"));
  if (parameters_Output.output_chemoattractant
      && chemoattractant)
    measurer_field.push_back(std::make_unique<Measurer_Field>(
        directory_data + "/Data_chemoattractant_" + model_name + domain_name +
        initial_condition_particles_name + initial_condition_fields_name +
        parameter_set_name + std::to_string(run_nr) + ".dat"));
  
  std::cout << std::setprecision(2) << std::scientific;
  std::cout << "\tDone!\n";
  
  if (parameters_Output.output_grid)
  {
    std::string filename_grid_void{
      directory_data + "/Data_grid_void_" + directory_data +
      domain_name + parameter_set_name + ".dat" };
    std::string filename_grid_solid{
      directory_data + "/Data_grid_solid_" + directory_data +
      domain_name + parameter_set_name + ".dat" };
    std::cout << "Outputing grid structure...\n";
    grid::print_centers(grid, void_solid.voids(), filename_grid_void);
    grid::print_centers(grid, void_solid.solids(), filename_grid_solid);
    std::cout << "\tDone!\n";
  }
  
  std::cout << "Starting dynamics...\n";
  for (std::size_t tt = 0; tt < measure_times.size(); ++tt)
  {
    double time = measure_times[tt];
    std::cout << "\ttime = " << time
              << "maximum = " << measure_times.back();
    while (ptrw.time() < time)
    {
      ptrw.step();
      for (std::size_t ii = 0; ii < solvers.size(); ++ii)
        solvers[ii]->step(concentration_fields[ii]);
      particle_grid.update(grid, ctrw.particles());
      for (std::size_t ii = 0; ii < solvers.size(); ++ii)
      {
        concentrations_particle[ii].update_by_grid();
        gradients_particle[ii].update_by_grid();
      }
      reaction(concentration_fields[0], particle_grid, grid);
    }
    for (std::size_t ii = 0; ii < measurer_particle.size(); ++ii)
      (*measurer_particle[ii])(ptrw, OutputState{}, time,
        Measurer_Particle::New{});
    for (std::size_t ii = 0; ii < measurer_field.size(); ++ii)
      (*measurer_field[ii])(concentration_fields[ii], void_solid.voids(), time);
  }
  std::cout << "\tDone!\n";
  
  return 0;
}
