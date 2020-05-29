//
//  Output.h
//  RunAndTumble
//
//  Created by Tomás Aquino on 29/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef Output_Bacteria_h
#define Output_Bacteria_h

#include <fstream>
#include <string>
#include "general/useful.h"
#include "Stochastic/CTRW/Measurer.h"
#include "Field/Measurer.h"

namespace bacteria
{
  namespace output_standard
  {
    using Measurer_Particle = ctrw::Measurer_state;
    using Measurer_Field = field::Measurer_scalar;
    
    struct OutputState
    {
      template <typename State, typename Stream>
      void operator()
      (State const& state, Stream& output,
       std::string const& delimiter = "\t")
      {
        useful::print(output, state.position, 0, delimiter);
        output << delimiter << state.orientation;
        output << delimiter << state.state;
      }
    };
    
    struct Parameters_Output
    {
      double time_min;
      double time_max;
      std::size_t nr_measures;
      int measure_spacing;
      bool output_grid;
      bool output_particle;
      bool output_nutrient;
      bool output_chemoattractant;
      
      void get(std::string const& filename)
      {
        std::ifstream file(filename);
        if (!file.is_open())
          throw std::runtime_error{ "Could not open " +
            filename + " for reading" };
        file >> time_min >> time_max
             >> nr_measures >> measure_spacing
             >> output_grid >> output_particle
             >> output_nutrient >> output_chemoattractant;
        if (file.fail())
          throw std::invalid_argument{ "Inappropriate input file " +
            filename };
      }
      
      static void print_parameter_list()
      {
        std::cout << "Output parameters:\n"
                  << "\ttime_min :               Time of first output\n"
                  << "\ttime_max :               Time of last output\n"
                  << "\tnr_measures :            Number of outputs\n"
                  << "\tmeasure_spacing :        0 - Linearly-spaced outputs;\n"
                  << "\t                         1 - Logarithmically spaced-outputs\n"
                  << "\toutput_grid :            0 - Do not output grid voids and solid positions\n"
                  << "\t                         1 - Output grid void and solid positions\n"
                  << "\toutput_particle_state :  0 - Do not output particle info\n"
                  << "\t                         1 - Output particle information\n"
                  << "\toutput_nutrient :        0 - Do not output nutrient info\n"
                  << "\t                         1 - Output nutrient info\n"
                  << "\toutput_chemoattractant : 0 - Do not output chemoattractant info\n"
                  << "\t                         1 - Output chemoattractant info (if chemoattractant on)\n";
      }
    };
  }
}

#endif /* Output_Bacteria_h */
