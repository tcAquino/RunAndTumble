//
//  StateGetter.h
//  BeadPack
//
//  Created by Tomás Aquino on 17/05/2020.
//  Copyright © 2020 Tomás Aquino. All rights reserved.
//

#ifndef StateGetter_h
#define StateGetter_h

#include <utility>
#include "general/Operations.h"

namespace ctrw
{
  template <typename Getter>
  struct Get_new_from_particle
  {
    Getter get;
    
    Get_new_from_particle(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    
    template <typename Particle>
    auto operator()(Particle const& particle) const
    {
      return get(particle.state_new());
    }
  };
  template <typename Getter> Get_new_from_particle(Getter&&) -> Get_new_from_particle<Getter>;
  
  template <typename Getter>
  struct Get_old_from_particle
  {
    Getter get;
    
    Get_old_from_particle(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    
    template <typename Particle>
    auto operator()(Particle const& particle) const
    {
      return get(particle.state_old());
    }
  };
  template <typename Getter> Get_old_from_particle(Getter&&) -> Get_old_from_particle<Getter>;
  
  template <typename Getter>
  struct Get_from_particle
  {
    Getter get;
    
    Get_from_particle(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    template <typename Particle>
    auto operator()(Particle const& particle) const
    {
      return get(particle.state_new(), particle.state_old());
    }
  };
  template <typename Getter> Get_from_particle(Getter&&) -> Get_from_particle<Getter>;
  
  struct Get_position
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.position; }
  };
  
  template <std::size_t dd>
  struct Get_position_component
  {
    template <typename State>
    auto operator()(State const& state) const
    { return operation::project<dd>(state.position); }
  };
  
  struct Get_time
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.time; }
  };
  
  struct Get_mass
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.mass; }
  };
  
  struct Get_velocity
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.velocity; }
  };
  
  template <typename std::size_t dd>
  struct Get_velocity_component
  {
    template <typename State>
    auto operator()(State const& state) const
    { return operation::project<dd>(state.velocity); }
  };
  
  template <typename Property>
  struct Get_position_property
  {
    Property property;
    
    Get_position_property(Property&& property = {})
    : property{ std::forward<Property>(property) }
    {}

    template <typename State>
    auto operator()(State const& state) const
    { return property(state.Position()); }
  };
  template <typename Property>
  Get_position_property(Property&&) ->
  Get_position_property<Property>;
  
  template <typename Property>
  struct Get_velocity_property
  {
    Property property;
    
    Get_velocity_property(Property&& property = {})
    : property{ std::forward<Property>(property) }
    {}

    template <typename State>
    auto operator()(State const& state) const
    { return property(state.velocity); }
  };
  template <typename Property>
  Get_velocity_property(Property&&) ->
  Get_velocity_property<Property>;
  
  struct Get_position_periodic
  {
    std::vector<double> domain_dimensions;
    
    Get_position_periodic(std::vector<double> domain_dimensions)
    : domain_dimensions{ domain_dimensions }
    {}
    
    template <typename State>
    auto operator()(State const& state) const
    {
      return operation::plus(state.position,
                             operation::times(domain_dimensions,
                                              state.periodicity));
    }
  };
  
  struct Get_tag
  {
    template <typename State>
    auto operator()(State const& state) const
    { return state.tag; }
  };
  
  template <typename Getter>
  struct Get_new
  {
    Getter get;
    
    Get_new(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    
    template <typename State>
    auto operator()(State const& state_new, State const&) const
    { return get(state_new); }
  };
  template <typename Getter> Get_new(Getter&&) -> Get_new<Getter>;
  
  
  template <typename Getter>
  struct Get_old
  {
    Getter get;
    
    Get_old(Getter&& get = {})
    : get{ std::forward<Getter>(get) }
    {}
    
    template <typename State>
    auto operator()(State const&, State const& state_old) const
    { return get(state_old); }
  };
  template <typename Getter> Get_old(Getter&&) -> Get_old<Getter>;
  
  template <typename Getter = Get_position,
  typename Getter_velocity = Get_velocity>
  struct Get_position_interp_velocity
  {
    double time;
    Getter get;
    Getter_velocity get_velocity;

    Get_position_interp_velocity
    (double time, Getter&& get = {}, Getter_velocity&& get_velocity = {})
    : time{ time }
    , get{ std::forward<Getter>(get) }
    , get_velocity{ get_velocity }
    {}

    template <typename State>
    auto operator()
    (State const& state_new, State const& state_old) const
    {
      double delta_t = time-state_old.time;
      auto vel = get_velocity(state_old.velocity);
      return operation::plus(get(state_old.position), vel*delta_t);
    }
  };
  template <typename Getter, typename Getter_velocity>
  Get_position_interp_velocity
  (double, Getter&&, Getter_velocity&&) ->
  Get_position_interp_velocity<Getter, Getter_velocity>;

  template <typename Getter = Get_time,
  typename VelocityMapper = useful::Forward<double>>
  struct Get_time_interp_velocity
  {
    double position;
    Getter get;
    VelocityMapper velocity_mapper;

    Get_time_interp_velocity
    (double position, Getter&& get = {},
     VelocityMapper&& velocity_mapper = {})
    : position{ position }
    , Getter{ std::forward<Getter>(get) }
    , velocity_mapper{ std::forward<VelocityMapper>(velocity_mapper) }
    {}

    template <typename State>
    auto operator()
    (State const& state_new, State const& state_old) const
    {
      double delta_x = position-state_old.position;
      double vel = velocity_mapper(state_old.velocity);
      return state_old.time+delta_x/vel;
    }
  };
  template <typename Getter, typename VelocityMapper>
  Get_time_interp_velocity
  (double, Getter&&, VelocityMapper&&) ->
  Get_time_interp_velocity<Getter, VelocityMapper>;
}

#endif /* StateGetter_h */
