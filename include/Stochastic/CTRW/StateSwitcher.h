//
//  StateSwitcher.h
//  RunAndTumble
//
//  Created by Tomás Aquino on 07/04/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef StateSwitcher_h
#define StateSwitcher_h

#include <cmath>

namespace ctrw
{
  class StateSwitcher_ConstantRate
  {
  public:
    double prob_run_to_tumble;
    double prob_tumble_to_run;
    double prob_wall_tumble_to_run;
    
    StateSwitcher_ConstantRate
    (double time_step,
     double rate_run_to_tumble, double rate_tumble_to_run,
     double rate_wall_tumble_to_run)
    : prob_run_to_tumble{ 1.-std::exp(-rate_run_to_tumble*time_step) }
    , prob_tumble_to_run{ 1.-std::exp(-rate_tumble_to_run*time_step) }
    , prob_wall_tumble_to_run{ 1.-std::exp(-rate_wall_tumble_to_run*time_step) }
    {}
    
    StateSwitcher_ConstantRate
    (double time_step, double rate_run_to_tumble, double rate_tumble_to_run)
    : StateSwitcher_ConstantRate{
      time_step, rate_run_to_tumble, rate_tumble_to_run, rate_tumble_to_run }
    {}
    
    template <typename State = useful::Empty>
    int run(State const& state = {})
    {
      if (std::bernoulli_distribution{prob_run_to_tumble}(rng))
        return 1;
      return 0;
    }
    
    template <typename State = useful::Empty>
    int tumble(State const& state = {})
    {
      if (std::bernoulli_distribution{prob_tumble_to_run}(rng))
        return 0;
      return 1;
    }
    
    template <typename State = useful::Empty>
    int wall_tumble(State const& state = {})
    {
      if (std::bernoulli_distribution{prob_wall_tumble_to_run}(rng))
        return 0;
      return 2;
    }
    
  private:
    std::mt19937 rng{ std::random_device{}() };
  };

  template <typename Concentration_particle>
  class StateSwitcher_RateConcentrationBracket
  {
  public:
    // Lower and upper bound of bracker
    std::pair<double, double> concentration_bracket;
    // Probabilities outside and inside bracket
    std::pair<double, double> prob_run_to_tumble;
    std::pair<double, double> prob_tumble_to_run;
    std::pair<double, double> prob_wall_tumble_to_run;
    
    StateSwitcher_RateConcentrationBracket
    (double time_step, std::pair<double, double> concentration_bracket,
     std::pair<double, double> rate_run_to_tumble,
     std::pair<double, double> rate_tumble_to_run,
     std::pair<double, double> rate_wall_tumble_to_run,
     Concentration_particle const& concentration_particle)
    : concentration_bracket{ concentration_bracket }
    , prob_run_to_tumble{
      1.-std::exp(-rate_run_to_tumble.first*time_step),
      1.-std::exp(-rate_run_to_tumble.second*time_step) }
    , prob_tumble_to_run{
      1.-std::exp(-rate_tumble_to_run.first*time_step),
      1.-std::exp(-rate_tumble_to_run.second*time_step) }
    , prob_wall_tumble_to_run{
      1.-std::exp(-rate_wall_tumble_to_run.first*time_step),
      1.-std::exp(-rate_wall_tumble_to_run.second*time_step) }
    , concentration_particle{ concentration_particle }
    {}
    
    StateSwitcher_RateConcentrationBracket
    (double time_step, std::pair<double, double> concentration_bracket,
     std::pair<double, double> rate_run_to_tumble,
     std::pair<double, double> rate_tumble_to_run,
     Concentration_particle const& concentration_particle)
    : StateSwitcher_RateConcentrationBracket{
      time_step, concentration_bracket,
      rate_run_to_tumble, rate_tumble_to_run, rate_tumble_to_run,
      concentration_particle }
    {}
    
    template <typename State>
    int run(State const& state)
    {
      double prob = inside_bracket(state)
      ? prob_run_to_tumble.second
      : prob_run_to_tumble.first;
      if (std::bernoulli_distribution{prob}(rng))
        return 1;
      return 0;
    }
    
    template <typename State>
    int tumble(State const& state)
    {
      double prob = inside_bracket(state)
      ? prob_tumble_to_run.second
      : prob_tumble_to_run.first;
      if (std::bernoulli_distribution{prob}(rng))
        return 0;
      return 1;
    }
    
    template <typename State>
    int wall_tumble(State const& state)
    {
      double prob = inside_bracket(state)
      ? prob_wall_tumble_to_run.second
      : prob_wall_tumble_to_run.first;
      if (std::bernoulli_distribution{prob}(rng))
        return 0;
      return 2;
    }
    
    template <typename State>
    bool inside_bracket(State const& state)
    {
      return
      concentration_particle(state) > concentration_bracket.first &&
      concentration_particle(state) < concentration_bracket.second;
    }
    
  private:
    std::mt19937 rng{ std::random_device{}() };
    Concentration_particle const& concentration_particle;
  };
}

#endif /* StateSwitcher_h */
