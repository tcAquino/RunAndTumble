//
//  InitialConditions.h
//  RunAndTumble
//
//  Created by Tomás Aquino on 22/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef InitialConditions_Bacteria_h
#define InitialConditions_Bacteria_h

#include <fstream>
#include <iostream>
#include <string>
#include <random>
#include <utility>
#include "general/Constants.h"

namespace bacteria
{
  template <typename Particle>
  auto make_particles_line_2d
  (std::size_t nr_particles, double position_xx,
   std::pair<double, double> const& bounds_yy)
  {
    using State = typename Particle::State;
    std::vector<Particle> particles;
    particles.reserve(nr_particles);
    
    std::mt19937 rng{ std::random_device{}() };
    std::uniform_real_distribution<double> dist;
    for (std::size_t part = 0; part < nr_particles; ++part)
    {
      double position_yy = bounds_yy.first +
        dist(rng)*(bounds_yy.second - bounds_yy.first);
      double orientation = constants::pi*(2.*dist(rng) - 1.);
      int state = 0;
      particles.push_back(State{
        { position_xx, position_yy },
        orientation,
        state,
        part
      });
    }
    
    return particles;
  };
  
  namespace initial_condition_particles_line_rectangle
  {
    char initial_condition_particles_name[16] = "line_rectangle";
    
    struct Parameters_ICBacteria
    {
      double total_mass;
      std::size_t nr_particles;
      double position_xx;
      
      void get(std::string const& filename)
      {
        std::ifstream file(filename);
        if (!file.is_open())
          throw useful::open_read_error(filename);
        file >> total_mass >> nr_particles >> position_xx;
        if (file.fail())
           throw useful::bad_file_contents(filename);
      }
      
      static void print_parameter_list()
      {
        std::cout << "ICBacteria parameters:\n"
                  << "\ttotal_mass :   Total initial mass of bacteria\n"
                  << "\tnr_particles : Initial number of bacteria\n"
                  << "\tposition_xx :  Initial longitudinal position of line injection\n";
      }
    };
    
    template <typename Particle>
    struct ParticleMaker
    {
      std::size_t nr_particles;
      double position_xx;
      
      ParticleMaker(Parameters_ICBacteria const& parameters)
      : nr_particles{ parameters.nr_particles }
      , position_xx{ parameters.position_xx }
      {}
      
      template <typename Domain>
      auto operator()(Domain const& domain) const
      {
        return make_particles_line_2d<Particle>(
         nr_particles,
         position_xx,
         { domain.box.corner[1],
           domain.box.corner[1]+domain.box.dimensions[1] });
      }
    };
  }
  
  namespace initial_condition_particles_line_ring
  {
    char initial_condition_particles_name[16] = "line_ring";
    
    struct Parameters_ICBacteria
    {
      double total_mass;
      std::size_t nr_particles;
      
      void get(std::string const& filename)
      {
        std::ifstream file(filename);
        if (!file.is_open())
          throw useful::open_read_error(filename);
        file >> total_mass >> nr_particles;
        if (file.fail())
          throw useful::bad_file_contents(filename);
      }
      
      static void print_parameter_list()
      {
        std::cout << "ICBacteria parameters:\n"
                  << "\ttotal_mass :   Total initial mass of bacteria\n"
                  << "\tnr_particles : Initial number of bacteria\n";
      }
    };
    
    template <typename Particle>
    struct ParticleMaker
    {
      std::size_t nr_particles;
      
      ParticleMaker(Parameters_ICBacteria const& parameters)
      : nr_particles{ parameters.nr_particles }
      {}
      
      template <typename Domain>
      auto operator()(Domain const& domain) const
      {
        return make_particles_line_2d<Particle>(
          nr_particles,
          0.,
          { domain.spheres[0].radius, domain.box.radius });
      }
    };
  }
  
  namespace initial_condition_fields_constant
  {
    char initial_condition_fields_name[16] = "constant";
    
    struct Parameters_ICEulerianFields
    {
      double nutrient_concentration;
      
      void get(std::string const& filename)
      {
        std::ifstream file(filename);
        if (!file.is_open())
          throw useful::open_read_error(filename);
        file >> nutrient_concentration;
        if (file.fail())
          throw useful::bad_file_contents(filename);
      }
      
      static void print_parameter_list()
      {
        std::cout << "ICEulerianFields parameters:\n"
                  << "\tnutrient_concentration : Homogeneous initial nutrient concentration\n";
      }
    };
    
    struct NutrientValue
    {
      NutrientValue
      (Parameters_ICEulerianFields const& parameters)
      : concentration{ parameters.nutrient_concentration }
      {}
      
      template <typename Position>
      auto operator()(Position const&) const
      { return concentration; }
      
      double concentration;
    };
    
    struct ChemoAttractantValue
    {
      ChemoAttractantValue
      (Parameters_ICEulerianFields const&)
      {}
      
      template <typename Position>
      auto operator()(Position const&) const
      { return 0.; }
    };
  }
}

#endif /* InitialConditions_Bacteria_h */
