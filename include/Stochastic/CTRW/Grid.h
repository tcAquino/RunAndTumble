//
//  Grid.h
//  RunAndTumble
//
//  Created by Tomás Aquino on 02/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef Grid_ctrw_h
#define Grid_ctrw_h

#include "general/useful.h"
#include <unordered_map>
#include <utility>
#include <vector>

namespace ctrw
{
  class GridParticles
  {
  public:
    using Container_particle = std::vector<std::size_t>;
    using Container_grid = std::unordered_map<std::size_t, std::vector<std::size_t>>;
    
    GridParticles()
    {}
    
    template <typename Grid, typename Particles>
    GridParticles(Grid const& grid, Particles const& particles)
    { update(grid, particles); }
    
    template <typename Grid, typename Particles>
    void update(Grid const& grid, Particles const& particles)
    {
      particle_grid.resize(particles.size());
      
      for (std::size_t part = 0; part < particles.size(); ++part)
        particle_grid[part] = grid.cell(particles[part].state_new().position);
      
      occupied_grid.clear();
      for (std::size_t part = 0; part < particles.size(); ++part)
      {
        auto it = occupied_grid.insert({ grid.cell(particles[part].state_new().position), {} }).first;
        it->second.push_back(part);
      }
    }
    
    template <typename Quantity>
    auto compute_by_grid(Quantity const& quantity) const
    {
      std::vector<typename Quantity::value_type> container;
      compute_by_grid(quantity, container);
      
      return container;
    }
    
    template <typename Quantity, typename Container>
    void compute_by_grid(Quantity const& quantity, Container& container) const
    {
      for (auto const& val : occupied_grid)
      {
        auto quantity_val = quantity(val.first);
        for (auto const& part : val.second)
          container[part] = quantity_val;
      }
    }
    
    template <typename Quantity>
    auto compute_by_particle(Quantity const& quantity) const
    {
      std::vector<typename Quantity::value_type> container;
      container.reserve(particle_grid.size());
      for (std::size_t part = 0; part < particle_grid.size(); ++part)
        container.push_back(quantity(particle_grid[part]));
      
      return container;
    }
    
    template <typename Quantity, typename Container>
    void compute_by_particle(Quantity const& quantity, Container& container) const
    {
      for (std::size_t part = 0; part < particle_grid.size(); ++part)
        container[part] = quantity(particle_grid[part]);
    }
                            
    template <typename Quantity>
    auto compute(Quantity const& quantity, std::size_t part) const
    { return quantity(particle_grid[part]); }
    
    std::size_t size() const
    { return particle_grid.size(); }
    
    Container_particle const& grid_idx_by_particle() const
    { return particle_grid; }
    
    Container_grid const& particles_by_grid_idx() const
    { return occupied_grid; }
    
  private:
    Container_particle particle_grid;
    Container_grid occupied_grid;
  };
  
  template <typename GridParticles, typename Quantity, typename Value_t = typename Quantity::value_type>
  class GridParticles_Quantity
  {
  public:
    using value_type = Value_t;
    struct Grid_update{};
    struct Particle_update{};
    
    GridParticles_Quantity
    (GridParticles const& particle_grid, Quantity const& quantity, Grid_update)
    : particle_grid{ particle_grid }
    , quantity{ std::move(quantity) }
    { update_by_grid(); }
    
    GridParticles_Quantity
    (GridParticles const& particle_grid, Quantity const& quantity, Particle_update)
    : particle_grid{ particle_grid }
    , quantity{ std::move(quantity) }
    { update_by_particle(); }
    
    void update_by_grid()
    {
      quantity_container.resize(particle_grid.size());
      particle_grid.compute_by_grid(quantity, quantity_container);
    }
    
    void update_by_particle()
    {
      quantity_container.resize(particle_grid.size());
      particle_grid.compute_by_particle(quantity, quantity_container);
    }
    
    template <typename State>
    value_type operator()(State const& state) const
    {
      return quantity_container[state.tag];
    }
    
  private:
    using Container_type = std::vector<value_type>;
    
    GridParticles const& particle_grid;
    Quantity const& quantity;
    Container_type quantity_container;
  };
}

#endif /* Grid_ctrw_h */
