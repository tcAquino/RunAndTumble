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
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include "general/useful.h"
#include "Geometry/Shape.h"

namespace bacteria
{
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
    
    using DomainShape = geometry::Parallelepiped<>;
    
    std::vector<geometry::Parallelepiped<>> make_rectangles
    (Parameters_Domain const& = {})
    {
      std::vector<geometry::Parallelepiped<>> rectangles;
      
      std::vector<double> corners_x{ 50., 60. };
      std::vector<double> corners_y_even{ 10. };
      std::vector<double> corners_y_odd{
        12.5, 22.5, 32.5, 42.5, 52.5, 62.5, 72.5, 82.5 };
      std::pair<double, double> column_width_height_even{ 5., 80. };
      std::pair<double, double> column_width_height_odd{ 120., 5. };
      for (std::size_t column = 0; column < corners_x.size(); ++column)
      {
          if (!(column%2))
            for (std::size_t row = 0; row < corners_y_even.size(); ++row)
              rectangles.push_back({
                { corners_x[column], corners_y_even[row] },
                { column_width_height_even.first,
                  column_width_height_even.second } });
          else
            for (std::size_t row = 0; row < corners_y_odd.size(); ++row)
              rectangles.push_back({
                { corners_x[column], corners_y_odd[row] },
                { column_width_height_odd.first,
                  column_width_height_odd.second } });
      }
      
      return rectangles;
    }
    
    geometry::Domain<DomainShape> make_domain(Parameters_Domain const& = {})
    {
      return{
        { { 0., 0. }, { 200., 100. } },
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
    
    using DomainShape = geometry::Sphere<>;
    double outer_radius = 50.;
    double inner_radius = 40.;
    
    geometry::Domain<DomainShape> make_domain(Parameters_Domain const& = {})
    {
      return {
        { { 0., 0. }, outer_radius },
        {},
        { { { 0., 0. }, inner_radius } } };
    }
  }
  
  namespace cif
  {
    std::string circle_tag = "R";
    std::string rectangle_tag = "B";
    std::string layer_tag = "L";
    std::string begin_tag = "DS";
    std::string end_tag = "DF";
    
    auto no_box
    (std::string const& filename, std::string const& expected)
    {
      return std::runtime_error{ "Bad domain box: Expected " +
        expected + ", found nothing" };
    }
    
    auto bad_box
    (std::string const& filename,
     std::string const& expected, std::string const& found)
    {
      return std::runtime_error{ "Bad domain box: Expected " +
        expected + ", found " + found };
    }
    
    struct Parameters_Domain
    {
      std::string filename_domain;
      double unit_conversion;
      std::size_t layer_inclusions;
      std::size_t layer_box;
      
      void get(std::string const& filename)
      {
        std::ifstream file(filename);
        if (!file.is_open())
          throw useful::open_read_error(filename);
        file >> filename_domain >> unit_conversion
             >> layer_inclusions >> layer_box;
        if (file.fail())
          throw useful::bad_file_contents(filename);
      }
      
      static void print_parameter_list()
      {
        std::cout << "Domain parameters:\n"
                  << "\tfilename_domain  : Name (with path) of the .cif file with domain info\n"
                  << "\tunit_conversion  : To be multiplied by values to convert units\n"
                  << "\tlayer_inclusions : Number of the layer declaring inclusions\n"
                  << "\tlayer_box        : Number of the layer declaring domain box\n";
      }
    };
    
    std::size_t get_layer
    (std::ifstream& file, std::string line)
    {
      std::size_t layer;
      std::stringstream stream(line.substr(3, line.length()));
      stream >> layer;
      if (stream.fail())
        throw useful::parse_error_line(line);
      return layer;
    }
    
    std::pair<std::size_t, bool> find_layer(std::ifstream& file)
    {
      std::string line;
      while(getline(file, line))
      {
        if (line.compare(0, layer_tag.length(), layer_tag) == 0)
          return { get_layer(file, line), 1 };
      }
      return { 0, 0 };
    }
    
    bool find_start(std::ifstream& file)
    {
      std::string line;
      while (getline(file, line))
        if (line.compare(0, begin_tag.length(), begin_tag) == 0)
          return 1;
      return 0;
    }
    
    template <typename Domain>
    std::pair<std::size_t, bool> read_inclusions
    (std::ifstream& file, Domain& domain,
     Parameters_Domain const& parameters)
    {
      std::string line;
      while (getline(file, line))
      {
        if (line.compare(0, rectangle_tag.length(), rectangle_tag) == 0)
        {
          std::cout << line << "\n";
          std::stringstream stream(line.substr(1, line.length()));
          double length_xx;
          double length_yy;
          double center_xx;
          double center_yy;
          stream >> length_xx >> length_yy >> center_xx >> center_yy;
          if (stream.fail())
            throw useful::parse_error(parameters.filename_domain, line);
          domain.parallelepipeds.push_back({
            { parameters.unit_conversion*(center_xx-length_xx/2.),
              parameters.unit_conversion*(center_yy-length_yy/2.) },
            { parameters.unit_conversion*length_xx,
              parameters.unit_conversion*length_yy } });
        }
        else if (line.compare(0, circle_tag.length(), circle_tag) == 0)
        {
          std::stringstream stream(line.substr(1, line.length()));
          double diameter;
          double center_xx;
          double center_yy;
          stream >> diameter >> center_xx >> center_yy;
          if (stream.fail())
            throw useful::parse_error(parameters.filename_domain, line);
          domain.spheres.push_back({
            { parameters.unit_conversion*center_xx,
              parameters.unit_conversion*center_yy },
            parameters.unit_conversion*diameter/2 });
        }
        else if (line.compare(0, layer_tag.length(), layer_tag) == 0)
          return { get_layer(file, line), 1 };
        else if (line.compare(0, end_tag.length(), end_tag) == 0)
          return { 0, 0 };
      }
      
      throw useful::bad_eof(parameters.filename_domain, end_tag);
    }
    
    template <typename Domain>
    std::pair<std::size_t, bool> read_box
    (std::ifstream&, Domain&, Parameters_Domain const&);
    
    template <>
    std::pair<std::size_t, bool> read_box
    (std::ifstream& file, geometry::Domain<geometry::Sphere<>>& domain,
     Parameters_Domain const& parameters)
    {
      std::string line;
      bool found_rectangle = 0;
      bool found_circle = 0;
      while(getline(file, line))
      {
        if (line.compare(0, rectangle_tag.length(), rectangle_tag) == 0)
          found_rectangle = 1;
        else if (line.compare(0, circle_tag.length(), circle_tag) == 0)
        {
          std::stringstream stream(line.substr(1, line.length()));
          double diameter;
          double center_xx;
          double center_yy;
          stream >> diameter >> center_xx >> center_yy;
          if (stream.fail())
            throw useful::parse_error(parameters.filename_domain, line);
          domain.box = {
            { parameters.unit_conversion*center_xx,
              parameters.unit_conversion*center_yy },
            parameters.unit_conversion*diameter/2 };
          found_circle = 1;
        }
        else if (line.compare(0, layer_tag.length(), layer_tag) == 0)
        {
          if (!found_rectangle)
            throw found_circle
            ? no_box(parameters.filename_domain, "rectangle")
            : bad_box(parameters.filename_domain, "rectangle", "circle");
          return { get_layer(file, line), 1 };
        }
        else if (line.compare(0, end_tag.length(), end_tag) == 0)
        {
          if (!found_rectangle)
            throw found_circle
            ? no_box(parameters.filename_domain, "rectangle")
            : bad_box(parameters.filename_domain, "rectangle", "circle");
          return { 0, 0 };
        }
      }
        
      throw useful::bad_eof(parameters.filename_domain, end_tag);
    }
    
    template <>
    std::pair<std::size_t, bool> read_box
    (std::ifstream& file, geometry::Domain<geometry::Parallelepiped<>>& domain,
     Parameters_Domain const& parameters)
    {
      std::string line;
      bool found_rectangle = 0;
      bool found_circle = 0;
      while(getline(file, line))
      {
        if (line.compare(0, rectangle_tag.length(), rectangle_tag) == 0)
        {
          std::stringstream stream(line.substr(1, line.length()));
          double length_xx;
          double length_yy;
          double center_xx;
          double center_yy;
          stream >> length_xx >> length_yy >> center_xx >> center_yy;
          domain.box = {
            { parameters.unit_conversion*(center_xx-length_xx/2),
              parameters.unit_conversion*(center_yy-length_yy/2.) },
            { parameters.unit_conversion*length_xx,
              parameters.unit_conversion*length_yy } };
          found_rectangle = 1;
        }
        else if (line.compare(0, circle_tag.length(), circle_tag) == 0)
          found_circle = 1;
        else if (line.compare(0, layer_tag.length(), layer_tag) == 0)
        {
          if (!found_rectangle)
            throw found_circle
            ? no_box(parameters.filename_domain, "rectangle")
            : bad_box(parameters.filename_domain, "rectangle", "circle");
          return { get_layer(file, line), 1 };
        }
        else if (line.compare(0, end_tag.length(), end_tag) == 0)
        {
          if (!found_rectangle)
          throw found_circle
          ? no_box(parameters.filename_domain, "rectangle")
          : bad_box(parameters.filename_domain, "rectangle", "circle");
          return { 0, 0 };
        }
      }
        
      throw useful::bad_eof(parameters.filename_domain, end_tag);
    }
    
    template <typename DomainShape>
    geometry::Domain<DomainShape> make_domain(Parameters_Domain const& parameters)
    {
      geometry::Domain<DomainShape> domain;

      std::ifstream file{ parameters.filename_domain };
      if (!file.is_open())
        throw useful::open_read_error(parameters.filename_domain);
      if (!find_start(file))
        throw useful::bad_eof(parameters.filename_domain, begin_tag);

      bool found_box = 0;
      auto layer = find_layer(file);
      while (layer.second == 1)
      {
       if (layer.first == parameters.layer_inclusions)
         layer = read_inclusions(file, domain, parameters);
       else if (layer.first == parameters.layer_box)
       {
         layer = read_box(file, domain, parameters);
         found_box = 1;
       }
      }

      if (!found_box)
        throw useful::bad_eof(parameters.filename_domain, "box layer");

      return domain;
    }
  }
  
  namespace domain_cif_rectangle
  {
    using namespace cif;
    
    char domain_name[64] = "cif_rectangle";
    
    using DomainShape = geometry::Parallelepiped<>;
    
    geometry::Domain<DomainShape> make_domain
    (cif::Parameters_Domain const& parameters)
    {
      return cif::make_domain<DomainShape>(parameters);
    }
  }
  
  namespace domain_cif_circle
  {
    using namespace cif;
    
    char domain_name[64] = "cif_circle";
    
    using DomainShape = geometry::Sphere<>;
    
    geometry::Domain<DomainShape> make_domain
    (cif::Parameters_Domain const& parameters)
    {
      return cif::make_domain<DomainShape>(parameters);
    }
  }
}

#endif /* Domains_Bacteria_h */
