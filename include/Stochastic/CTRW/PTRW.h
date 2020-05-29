//
//  PTRW.h
//  CTRW
//
//  Created by Tomas Aquino on 1/9/19.
//  Copyright © 2019 Tomas Aquino. All rights reserved.
//

// Interface for PTRW using CTRW object framework

#ifndef PTRW_h
#define PTRW_h

#include <list>

namespace ctrw
{
	template <typename CTRW, typename Transitions, typename Time_t = double>
	class PTRW
	{
	public:
		using Particle = typename CTRW::Particle;
		using Container = typename CTRW::Container;
    using Time = Time_t;

    PTRW(CTRW& ctrw, Transitions& transitions, Time time = 0)
    : ctrw{ ctrw }
    , transitions{ transitions }
    , current_time{ time }
		{
      set_step(transitions.timestep());
    }
    
    PTRW(CTRW& ctrw, Transitions& transitions, Time time_step, Time time)
    : ctrw{ ctrw }
    , transitions{ transitions }
    , current_time{ time }
    , time_step{ time_step }
    {}

		void inject(Particle particle)
		{ ctrw.inject(particle); }

		template <typename Criterium>
		void remove(Criterium criterium)
		{ ctrw.remove(criterium); }

		template <typename IntegerType>
		void remove(std::list<IntegerType>& list)
		{ ctrw.remove(list); }

		void step(Time time_step)
		{
			transitions.timestep(time_step);
			ctrw.step(transitions);
      current_time += time_step;
		}

    void step()
    {
      ctrw.step(transitions);
      current_time += time_step;
    }

    void evolve_by_particle(Time time_max)
    {
      ctrw.evolve_time(time_max, transitions);
      current_time = time_max;
    }
    
    void evolve(Time time_max)
    {
      while (current_time < time_max)
        step();
      current_time = time_max;
    }

    void set_step(Time time_step)
    {
      this->time_step = time_step;
      transitions.timestep(time_step);
    }

    double time() const
    { return current_time; }

		Container const& particles() const
		{ return ctrw.particles(); }

		Container const& particles(std::size_t part) const
		{ return ctrw.particles(part); }

		std::size_t size() const
		{ return ctrw.size(); }

	private:
		CTRW& ctrw;
		Transitions& transitions;

    Time current_time;
    Time time_step;
  };
}

#endif /* PTRW_h */