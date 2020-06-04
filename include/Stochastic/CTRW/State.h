//
//  State.h
//  CTRW_2
//
//  Created by Tomas Aquino on 9/26/17.
//  Copyright © 2017 Tomas Aquino. All rights reserved.
//

#ifndef State_h
#define State_h

#include <valarray>
#include <vector>
#include "general/Operations.h"
#include "general/useful.h"

// Particle states for tdrw/ptrw/ctrw models

namespace ctrw
{
  template
  <typename Position_t, typename Periodicity_t = std::vector<int>, typename Mass_t = double,
  typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_periodic
  {
    using Position_type = Position_t;
    using Periodicity_type = Periodicity_t;
    using Mass_type = Mass_t;
    using Time_type = Time_t;
    using Tag_type = Tag_t;
    
    State_periodic()
    {}
    
    State_periodic
    (Position_t position, Periodicity_type periodicity, Mass_type mass = 1,
     Time_t time = 0, Tag_t tag = {})
    : position{ position }
    , mass{ mass }
    , periodicity{ periodicity }
    , time{ time }
    , tag{ tag }
    {}

    Position_type position;
    Periodicity_type periodicity;
    Mass_type mass;
    Time_type time;
    Tag_type tag;
  };
  
  template
  <typename Position_t, typename Velocity_t,
  typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_velocity
  {
    using Position_type = Position_t;
    using Velocity_type = Velocity_t;
    using Time_type = Time_t;
    using Tag_type = Tag_t;
    
    State_velocity()
    {}
    
    State_velocity
    (Position_t position, Velocity_t velocity,
     Time_t time = {}, Tag_t tag = {})
    : position(position)
    , velocity(velocity)
    , time(time)
    , tag(tag)
    {}

    Position_type position;
    Velocity_type velocity;
    Time_type time;
    Tag_type tag;
  };

  template <typename Position_t, typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_position
  {
    using Position_type = Position_t;
    using Time_type = Time_t;
    using Tag_type = Tag_t;

    State_position()
    {}

    State_position
    (Position_t position, Time_t time = {}, Tag_t tag = {})
    : position(position)
    , time(time)
    {}

    Position_type position{};
    Time_type time{};
    Tag_type tag{};
  };
  
  template <typename Position_t, typename Mass_t = double, typename Time_t = double, typename Tag_t = useful::Empty>
  struct State_position_mass
  {
    using Position_type = Position_t;
    using Mass_type = Mass_t;
    using Time_type = Time_t;
    using Tag_type = Tag_t;

    State_position_mass()
    {}

    State_position_mass
    (Position_t position, Mass_type mass = 1., Time_t time = {}, Tag_t tag = {})
    : position{ position }
    , mass{ mass }
    , time{ time }
    {}

    Position_type position;
    Mass_type mass;
    Time_type time;
    Tag_type tag;
  };

  template <typename Position_t>
  struct State_PTRW_position
  {
    using Position_type = Position_t;

    State_PTRW_position()
    {}

    State_PTRW_position(Position_t position)
    : position(position)
    {}

    Position_type position{};
  };

  template <typename Position_t>
  struct State_PTRW_position_type
  {
    using Position_type = Position_t;

    State_PTRW_position_type()
    {}

    State_PTRW_position_type(std::size_t type)
    : type(type)
    {}

    State_PTRW_position_type
    (Position_type position, std::size_t type = 0)
    : position(position)
    , type(type)
    {}

    Position_type position{};
    std::size_t type{0};
  };

  template <typename Node_t>
  struct State_Node
  {
    Node_t const* node;
    double time{ 0. };

    auto Position() const
    { return node->position; }
  };

  template <typename Node_t>
  struct State_Node_Mass
  {
    Node_t const* node;
    double time{ 0. };
    double mass{ 1. };

    auto Position() const
    { return node->position; }
  };

  template <typename Node_t>
  struct State_Node_ReactionTime
  {
    Node_t const* node;
    double time{ 0. };
    double reaction_time{ -1. };

    auto Position() const
    { return node->position; }
  };

  template <typename Node_t, typename Position_t>
  struct State_Node_Position
  {
    using Position_type = Position_t;
    using Time_type = double;
    Node_t const* node;
    Position_t position;
    double time;

    State_Node_Position(Node_t const* node, double time = 0.)
    : node(node)
    , position(node->position)
    , time(time)
    {}

    Position_type Position() const
    { return position; }
  };

  template <typename Position_t>
  struct State_Position_ReactionTime
  {
    using Position_type = Position_t;
    using Time_type = double;

    State_Position_ReactionTime()
    {}

    State_Position_ReactionTime(Position_t position, Time_type time = {})
    : position(position)
    , time(time)
    {}
    
    Position_type position{};
    Time_type time{};
    Time_type reaction_time{ -1. };
  };
  
  template <typename Orientation_t, typename Tag_t = useful::Empty>
  struct State_RunTumble
  {
    using Position_type = std::vector<double>;
    using Orientation_type = Orientation_t;
    using Time_type = double;
    using Tag_type = Tag_t;
    
    State_RunTumble()
    {}
    
    State_RunTumble
    (Position_type position, Orientation_type orientation, bool run = 0,
     Time_type time = {}, Tag_type tag = {})
    : position{ position }
    , orientation{ orientation }
    , run{ run }
    , time{ time }
    , tag{ tag }
    {}
    
    Position_type position{};
    Orientation_type orientation{};
    bool run;
    Time_type time{};
    Tag_type tag;
  };
  
  template <typename Orientation_t = double, typename Tag_t = std::size_t,
  typename Position_t = std::vector<double>>
  struct State_RunTumble_PTRW
  {
    using Position_type = Position_t;
    using Orientation_type = Orientation_t;
    using Tag_type = Tag_t;
    
    State_RunTumble_PTRW()
    {}
    
    State_RunTumble_PTRW
    (Position_type position, Orientation_type orientation,
     int state, Tag_type tag = {})
    : position{ position }
    , orientation{ orientation }
    , state{ state }
    , tag{ tag }
    {}
    
    Position_type position{};
    Orientation_type orientation{};
    int state;
    Tag_type tag;
  };
}

#endif /* State_h */
