//
//  Shape.h
//  RunAndTumble
//
//  Created by Tomás Aquino on 07/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef Shape_h
#define Shape_h

#include <cmath>
#include <vector>
#include "general/Operations.h"

namespace shape
{
  template <typename Container_t = std::vector<double>>
  struct Paralleliped
  {
    using Position_type = Container_t;
    Position_type corner;
    Position_type dimensions;
    
    Position_type half_dimensions{ operation::div_scalar(dimensions,2.) };
    Position_type center{ operation::plus(corner, half_dimensions) };
    
    Paralleliped(Position_type corner, Position_type dimensions)
    : corner{ corner }
    , dimensions{ dimensions }
    {}
      
    std::size_t dim() const
    { return dimensions.size(); }
    
    template <typename Position>
    bool inside(Position const& position) const
    {
      operation::minus(position, center, position_helper);
      for (std::size_t dd = 0; dd < dim(); ++dd)
        if (std::abs(position_helper[dd]) > half_dimensions[dd])
          return false;
      return true;
    }
      
  private:
    mutable Position_type position_helper{ Position_type(dimensions.size()) };
  };
  
  template <typename Container_t = std::vector<double>>
  struct Sphere
  {
    using Position_type = Container_t;
    Position_type center;
    double radius;

    Sphere(Position_type center, double radius)
    : center{ center }
    , radius{ radius }
    {}

    std::size_t dim() const
    { return center.size(); }

    template <typename Position>
    bool inside(Position const& position) const
    {
      operation::minus(position, center, position_helper);
      if (operation::abs_sq(position_helper) > radius*radius)
        return false;
      return true;
    }

  private:
      mutable Position_type position_helper{ Position_type(center.size()) };
  };
  
  template <typename Position, typename Sphere>
  void ReflectOffSphere_OutsideToInside
  (Position const& position_new, Position const& position_old, Sphere const& sphere,
   Position& reflected, Position& contact_point)
  {
    auto insideminusoutside = operation::minus(position_new, position_old);
    auto outsideminuscenter = operation::minus(position_old, sphere.center);
    
    double coeff_a = operation::abs_sq(insideminusoutside);
    double coeff_b = 2.*operation::dot(insideminusoutside, outsideminuscenter);
    double coeff_c = operation::dot(outsideminuscenter, outsideminuscenter) - sphere.radius*sphere.radius;
    double fraction_to_intersection = -(coeff_b +
                                        std::sqrt(coeff_b*coeff_b - 4.*coeff_a*coeff_c))/(2.*coeff_a);
    operation::linearOp(fraction_to_intersection, insideminusoutside, position_old, contact_point);
    
    auto leftover = operation::times_scalar(1.-fraction_to_intersection, insideminusoutside);
    auto normal = operation::minus(contact_point, sphere.center);
    operation::div_scalar_InPlace(normal, sphere.radius);
    
    operation::linearOp(-2.*operation::dot(leftover, normal), normal, position_new, reflected);
  }
  
  template <typename Position, typename Sphere>
  void ReflectOffSphere_OutsideToInside
  (Position const& position_new, Position const& position_old, Sphere const& sphere,
   Position& reflected)
  {
    auto contact_point(position_new.size());
    ReflectOffSphere_OutsideToInside(position_new, position_old, sphere, reflected, contact_point);
  }
}

#endif /* Shape_h */
