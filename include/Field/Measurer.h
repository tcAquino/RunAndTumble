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
  class Measurer
  {
  public:
    Measurer
    (std::string const& filename,
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
      
    ~Measurer()
    { output.close(); }
      
    template <typename Field, typename Container>
    void operator()
    (Field const& field, Container const& points)
    {
      bool delim = 0;
      for (auto idx : points)
      {
        useful::print(output, field[idx], delim, delimiter);
        delim = 1;
      }
      output << "\n";
    }
      
    template <typename Field, typename Value, typename Container>
    void operator()
    (Field const& field, Container const& points, Value tag)
    {
      output << tag << delimiter;
      (*this)(field, points);
    }
    
  private:
    std::ofstream output;
    std::string delimiter;
  };
}

#endif /* Measurer_Fields_h */
