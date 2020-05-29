//
//  Measurer.h
//
//  Created by Tomas Aquino on 29/5/20.
//  Copyright © 2020 Tomas Aquino. All rights reserved.
//

#ifndef Measurer_Fields_h
#define Measurer_Fields_h

#include <valarray>
#include <fstream>
#include <string>
#include "general/useful.h"

namespace field
{
  class Measurer_scalar
  {
  public:
    Measurer_scalar
    (std::string filename,
     int precision = 8, std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      output << std::setprecision(precision)
             << std::scientific;
      if (!output.is_open())
        throw std::runtime_error{
          "Could not open file " + filename + "for writing" };
    }
      
    ~Measurer_scalar()
    { output.close(); }
      
    template <typename Field>
    void operator()
    (double time, Field const& field,
     std::vector<std::size_t> const& points)
    {
      output << time;
      for (auto idx : points)
        output << delimiter << field[idx];
      output << "\n";
    }
    
  private:
    std::fstream output;
    std::string delimiter;
  };
}

#endif /* Measurer_Fields_h */
