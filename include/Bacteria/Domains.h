//
//  Domains.h
//  RunAndTumble
//
//  Created by Tomás Aquino on 22/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef Domains_Bacteria_h
#define Domains_Bacteria_h

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <utility>
#include <vector>
#include "general/Shape.h"
#include "general/useful.h"

namespace bacteria
{
  template <typename Shape>
  struct Domain
  {
    using DomainShape = Shape;
    
    Shape box;
    std::vector<shape::Paralleliped<>> parallelipeds;
    std::vector<shape::Sphere<>> spheres;
  };
  
  class Grid_void_solid
  {
  public:
    template <typename Domain, typename KDTree>
    Grid_void_solid(Domain const& domain, KDTree const& kdtree_grid)
    {
      for (auto const& shape : domain.parallelipeds)
        kdtree_grid.inside(shape, solid_idx);
      for (auto const& shape : domain.spheres)
        kdtree_grid.inside(shape, solid_idx);
      kdtree_grid.outside(domain.box, solid_idx);
      std::sort(solid_idx.begin(), solid_idx.end());
      
      for (std::size_t idx = 0; idx < kdtree_grid.size(); ++idx)
        if (!useful::contains(solid_idx, idx))
          void_idx.push_back(idx);
    }
    
    std::size_t nr_voids() const
    { return void_idx.size(); }
    
    std::size_t nr_solids() const
    { return solid_idx.size(); }
    
    std::size_t voids(std::size_t idx) const
    { return void_idx[idx]; }
    
    std::size_t solids(std::size_t idx) const
    { return solid_idx[idx]; }
    
    auto const& voids() const
    { return void_idx; }
    
    auto const& solids() const
    { return solid_idx; }
    
  private:
    std::vector<std::size_t> void_idx;
    std::vector<std::size_t> solid_idx;
  };
  
  namespace domain_test_rectangle
  {
    char domain_name[64] = "test_rectangle";
    
    struct Parameters_Domain
    {
      void get(std::string const&){}
      
      static void print_parameter_list()
      {
        std::cout << "Domain parameters:\n"
                  << "\tNo parameters\n";
      }
    };
    
    using DomainShape = shape::Paralleliped<>;
    std::vector<double> domain_dimensions = { 200., 100. };
    
    std::vector<shape::Paralleliped<>> make_rectangles(Parameters_Domain const& = {})
    {
      std::vector<shape::Paralleliped<>> rectangles;
      
      std::vector<double> corners_x{ 50., 60. };
      std::vector<double> corners_y_even{ 10. };
      std::vector<double> corners_y_odd{ 12.5, 22.5, 32.5, 42.5, 52.5, 62.5, 72.5, 82.5 };
      std::pair<double, double> column_width_height_even{ 5., 80. };
      std::pair<double, double> column_width_height_odd{ 120., 5. };
      for (std::size_t column = 0; column < corners_x.size(); ++column)
      {
          if (!(column%2))
            for (std::size_t row = 0; row < corners_y_even.size(); ++row)
              rectangles.push_back({
                { corners_x[column], corners_y_even[row] },
                { column_width_height_even.first, column_width_height_even.second } });
          else
            for (std::size_t row = 0; row < corners_y_odd.size(); ++row)
              rectangles.push_back({
                { corners_x[column], corners_y_odd[row] },
                { column_width_height_odd.first, column_width_height_odd.second } });
      }
      
      return rectangles;
    }
    
    Domain<DomainShape> make_domain(Parameters_Domain const& = {})
    {
      return{
        { { 0., 0. }, domain_dimensions },
        make_rectangles(),
        {} };
    }
  }
  
  namespace domain_test_ring
  {
    char domain_name[64] = "test_ring";
    
    struct Parameters_Domain
    {
      void get(std::string const&){}
      
      static void print_parameter_list()
      {
        std::cout << "Domain parameters:\n"
                  << "\tNo parameters\n";
      }
    };
    
    using DomainShape = shape::Sphere<>;
    double outer_radius = 50.;
    double inner_radius = 40.;
    std::vector<double> domain_dimensions{ 2.*outer_radius, 2.*outer_radius };
    
    Domain<DomainShape> make_domain(Parameters_Domain const& = {})
    {
      return {
        { { 0., 0. }, outer_radius },
        {},
        { { { 0., 0. }, inner_radius } } };
    }
  }
}

#endif /* Domains_Bacteria_h */
