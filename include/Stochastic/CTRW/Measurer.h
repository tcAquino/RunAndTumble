//
//  Measurer.h
//  CTRW_velocity
//
//  Created by Tomas Aquino on 14/1/20.
//  Copyright © 2020 Tomas Aquino. All rights reserved.
//

#ifndef Measurer_CTRW_h
#define Measurer_CTRW_h

#include <valarray>
#include <fstream>
#include <string>
#include "general/useful.h"

namespace ctrw
{
  class Measurer_state
  {
  public:
    Measurer_state
    (std::string filename,
     int precision = 8, std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      output << std::setprecision(precision)
             << std::scientific;
      if (!output.is_open())
        throw std::runtime_error{
          "Could not open file " + filename + " for writing" };
    }
    
    ~Measurer_state()
    { output.close(); }
      
    template <typename CTRW, typename OutputState>
    void operator()
    (CTRW const& ctrw, OutputState output_state)
    {
      bool delim;
      for (auto const& part : ctrw.particles())
      {
        if (delim)
          output << delimiter;
        output_state(part.state_new(), output, delimiter);
        delim = 1;
      }
      output << "\n";
    }
      
    template <typename CTRW, typename OutputState, typename Value>
    void operator()
    (CTRW const& ctrw, OutputState output_state, Value tag)
    {
      output << tag << delimiter;
      (*this)(ctrw, output_state);
    }
    
  private:
    std::fstream output;
    std::string delimiter;
  };
      
  class Measurer_collection
  {
  public:
    Measurer_collection
    (std::string const& filename, int precision = 8, std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      output << std::setprecision(precision)
             << std::scientific;
      if (!output.is_open())
      throw std::runtime_error{
        "Could not open file " + filename + "for writing" };
    }
      
    ~Measurer_collection()
    { output.close(); }
    
    template <typename Subject>
    void operator()(Subject const& subject)
    {
      std::string delim = "";
      for (std::size_t ii = 0; ii < subject.size(); ++ii)
      {
        output << delim << subject.size(ii);
        for (auto it = subject.cbegin(ii); it != subject.cend(ii); ++it)
          output << delimiter << *it;
        delim = delimiter;
      }
      output << "\n";
    }
      
    template <typename Subject, typename Value>
    void operator()(Subject const& subject, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject);
    }
    
  private:
    std::ofstream output;
    const std::string delimiter;
  };
  
  class Measurer_container
  {
  public:
    Measurer_container
    (std::string const& filename, int precision = 8, std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      output << std::setprecision(precision)
             << std::scientific;
      if (!output.is_open())
      throw std::runtime_error{
        "Could not open file " + filename + "for writing" };
    }
    
    ~Measurer_container()
    { output.close(); }
    
    template <typename Subject, typename Getter>
    void operator()(Subject const& subject, Getter get)
    {
      bool delim = 0;
      for (auto const& part : subject.particles())
      {
        auto values = get(part);
        useful::print(output, values, delim, delimiter);
        delim = 1;
      }
      output << "\n";
    }
    
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
  private:
    std::ofstream output;
    const std::string delimiter;
  };
  
  class Measurer_value
  {
  public:
    Measurer_value
    (std::string const& filename, int precision = 8, std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      output << std::setprecision(precision);
      output << std::scientific;
      if (!output.is_open())
      throw std::runtime_error{
        "Could not open file " + filename + "for writing" };
    }
    
    ~Measurer_value()
    { output.close(); }
    
    template <typename Subject, typename Getter>
    void operator()(Subject const& subject, Getter get)
    {
      std::string delim = "";
      for (auto const& part : subject.particles())
      {
        output << delim << get(part);
        delim = delimiter;
      }
      output << "\n";
    }
    
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
  private:
    std::ofstream output;
    const std::string delimiter;
  };
      
  class Measurer_value_total
  {
  public:
    Measurer_value_total
    (std::string const& filename, int precision = 8, std::string delimiter = "\t")
    : output{ filename }
    , delimiter{ delimiter }
    {
      output << std::setprecision(precision)
             << std::scientific;
      if (!output.is_open())
      throw std::runtime_error{
        "Could not open file " + filename + "for writing" };
    }
    
    ~Measurer_value_total()
    { output.close(); }
    
    template <typename Subject, typename Getter>
    void operator()(Subject const& subject, Getter get)
    {
      decltype(get(subject.particles().back())) val{ 0. };
      for (auto const& part : subject.particles())
        val += get(part);
      output << val << "\n";
    }
    
    template <typename Subject, typename Getter, typename Value>
    void operator()(Subject const& subject, Getter get, Value tag)
    {
      output << tag << delimiter;
      (*this)(subject, get);
    }
    
    private:
      std::ofstream output;
      const std::string delimiter;
  };
  
//  class Measurer_crossing_dist_cprint
//  {
//  public:
//    Measurer_crossing_dist_cprint
//    (std::valarray<double> measure_at, std::size_t nr_particles = 0)
//    : measure_at{ measure_at }
//    , measure_value(measure_at.size())
//    {
//      for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//        measure_value[mm].reserve(nr_particles);
//    }
//
//    template <typename Subject, typename Getter>
//    void update(Subject const& subject, Getter const& get)
//    {
//      for (auto const& part : subject.particles())
//        for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//          if ((part.state_new().position[0]-measure_at[mm]) *
//              (part.state_old().position[0]-measure_at[mm]) <= 0.)
//            measure_value[mm].push_back(get(part));
//    }
//
//    std::size_t size(std::size_t mm) const
//    { return measure_value[mm].size(); }
//
//    void print(FILE* fout) const
//    {
//      for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//      {
//        fprintf(fout, "%15.6e", measure_at[mm]);
//        for (auto const& val : measure_value[mm])
//          fprintf(fout, "%15.6e", val);
//        fprintf(fout, "\n");
//      }
//    }
//
//    const std::valarray<double> measure_at;
//
//  private:
//    std::vector<std::vector<double>> measure_value;
//  };
//
//  class Measurer_dist_cprint
//  {
//  public:
//    Measurer_dist_cprint
//    (std::valarray<double> measure_at, std::size_t nr_particles = 0)
//    : measure_at(measure_at)
//    , measure_value(measure_at.size())
//    {
//      for (std::size_t tt = 0; tt < measure_at.size(); ++tt)
//        measure_value[tt].reserve(nr_particles);
//    }
//
//    template <typename Subject, typename Getter>
//    void update(Subject const& subject, std::size_t mm, Getter get)
//    {
//      for (auto const& part : subject.particles())
//        measure_value[mm].push_back(get(part));
//    }
//
//    std::size_t size(std::size_t mm) const
//    { return measure_value[mm].size(); }
//
//    void print(FILE* fout) const
//    {
//      for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//      {
//        fprintf(fout, "%15.6e", measure_at[mm]);
//        for (auto const& val : measure_value[mm])
//          fprintf(fout, "%15.6e", val);
//        fprintf(fout, "\n");
//      }
//    }
//
//    const std::valarray<double> measure_at;
//
//  private:
//    std::vector<std::vector<double>> measure_value;
//  };
//
//  class Measurer_mean_cprint
//  {
//  public:
//    Measurer_mean_cprint(std::valarray<double> measure_at)
//    : measure_at(measure_at)
//    , mean(measure_at.size())
//    {}
//
//    template <typename Subject, typename Getter>
//    void update(Subject const& subject, std::size_t mm, Getter get)
//    {
//      double val = 0.;
//      for (auto const& part : subject.particles())
//        val += get(part);
//      mean[mm] = val/subject.nr_particles();
//    }
//
//    std::size_t size(std::size_t mm) const
//    { return mean.size(); }
//
//    void print(FILE* fout) const
//    {
//      for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//        fprintf(fout, "%15.6e%15.6e\n", measure_at[mm], mean[mm]);
//    }
//
//    const std::valarray<double> measure_at;
//
//  private:
//    std::valarray<double> mean;
//  };
//
//  class Measurer_crossing_total_cprint
//  {
//  public:
//    Measurer_crossing_total_cprint
//    (std::valarray<double> measure_at)
//    : measure_at(measure_at)
//    , measure_value(measure_at.size())
//    {}
//
//    template <typename Subject, typename Getter>
//    void update(Subject const& subject, Getter const& get)
//    {
//      for (auto const& part : subject.particles())
//        for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//          if ((part.state_new().position[0]-measure_at[mm]) *
//              (part.state_old().position[0]-measure_at[mm]) <= 0.)
//          {
//            measure_value[mm] += get(part);
//            ++counts[mm];
//          }
//    }
//
//    std::size_t size(std::size_t mm) const
//    { return counts[mm]; }
//
//    void print(FILE* fout) const
//    {
//      for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//        fprintf(fout, "%15.6e%15.6e\n", measure_at[mm],measure_value[mm]);
//    }
//
//    const std::valarray<double> measure_at;
//
//  private:
//    std::vector<double> measure_value;
//    std::vector<std::size_t> counts{ std::vector<std::size_t>(measure_value.size()) };
//  };
//
//  class Measurer_total_cprint
//  {
//  public:
//    Measurer_total(std::valarray<double> measure_at)
//    : measure_at(measure_at)
//    , total(measure_at.size())
//    {}
//
//    template <typename Subject, typename Getter>
//    void update(Subject const& subject, std::size_t mm, Getter get)
//    {
//      for (auto const& part : subject.particles())
//        total[mm] += get(part);
//    }
//
//    std::size_t size(std::size_t mm) const
//    { return total.size(); }
//
//    void print(FILE* fout) const
//    {
//      for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//        fprintf(fout, "%15.6e%15.6e\n", measure_at[mm], total[mm]);
//    }
//
//    const std::valarray<double> measure_at;
//
//  private:
//    std::valarray<double> total;
//  };
//
//  class Measurer_mean_variance_cprint
//  {
//  public:
//    Measurer_mean_variance_cprint
//    (std::valarray<double> measure_at)
//    : measure_at(measure_at)
//    , mean(measure_at.size())
//    , variance(measure_at.size())
//    {}
//
//    template <typename Subject, typename Getter>
//    void update(Subject const& subject, std::size_t mm, Getter get)
//    {
//      double val;
//      double val_mean = 0.;
//      double val_mean_sq = 0.;
//      for (auto const& part : subject.particles())
//      {
//        val = get(part);
//        val_mean += val;
//        val_mean_sq += val*val;
//      }
//      std::size_t nr_particles = subject.nr_particles();
//      mean[mm] = val_mean/nr_particles;
//      variance[mm] = val_mean_sq/nr_particles-mean[mm]*mean[mm];
//    }
//
//    std::size_t size(std::size_t mm) const
//    { return mean.size(); }
//
//    void print(FILE* fout) const
//    {
//      for (std::size_t mm = 0; mm < measure_at.size(); ++mm)
//        fprintf(fout, "%20.12e%20.12e%20.12e\n", measure_at[mm], mean[mm], variance[mm]);
//    }
//
//    const std::valarray<double> measure_at;
//
//  private:
//    std::valarray<double> mean;
//    std::valarray<double> variance;
//  };
//
//  template<typename Target_type, typename Return_type>
//  class Measurer_first_return_cprint
//  {
//  public:
//    Measurer_first_return_cprint
//    (Target_type target, std::size_t nr_particles,
//     std::size_t nr_measures = 0, Return_type initial = 0.)
//    : target{ target }
//    , return_last(nr_particles, initial)
//    {
//      return_data.reserve(nr_measures);
//    }
//
//    Measurer_first_return
//    (Target_type target, std::size_t nr_particles,
//     std::size_t nr_measures, std::vector<Return_type> initial)
//    : target{ target }
//    , return_last{ initial }
//    {
//      return_data.reserve(nr_measures);
//    }
//
//    template
//    <typename Subject, typename Getter_new, typename Getter_old,
//    typename Getter_return_new, typename Getter_return_old>
//    void update
//    (Subject const& subject, Getter_new getter_new, Getter_old getter_old,
//      Getter_return_new getter_return_new, Getter_return_old getter_return_old)
//    {
//      for (std::size_t pp = 0; pp < subject.nr_particles(); ++pp)
//      {
//        auto const& part = subject.particles(pp);
//        auto val_new = getter_new(part);
//        auto val_old = getter_old(part);
//        if (val_new == target && val_old != target)
//        {
//          return_data.push_back(getter_return_new(part) - return_last[pp]);
//          ++nr_counts;
//        }
//        else if (val_new != target && val_old == target)
//        {
//          return_last[pp] = getter_return_new(part);
//        }
//      }
//    }
//
//    std::size_t counts()
//    { return nr_counts; }
//
//    void print(FILE* fout) const
//    {
//      for (auto const& val : return_data)
//        fprintf(fout, "%15.6e\n", val);
//    }
//
//  private:
//    const Target_type target;
//    std::vector<Return_type> return_last;
//    std::vector<Return_type> return_data;
//    std::size_t nr_counts{ 0 };
//  };
}

#endif /* Measurer_CTRW_h */
