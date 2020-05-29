//
// useful.h
//
// Created on: Mar 15, 2011
// Author: tomas
//

// Miscelaneous collection of useful objects and algorithms

#ifndef USEFUL_H_
#define USEFUL_H_

#include <algorithm>
#include <cstdlib>
#include <cstring>
#include <string>
#include <functional>
#include <utility>

namespace useful
{
  std::size_t factorial(std::size_t nn)
  {
    if (nn == 0)	 return 1;
    return nn * factorial(nn - 1);
  }

  std::size_t factorial_incomplete(std::size_t nn, std::size_t mm)
  {
    std::size_t result = 1.;
    ++nn;
    for (std::size_t ii = 0; ii < mm && nn --> 0; ++ii) result *= nn;

    return result;
  }

  // Check if container contains val
  // Warning: container must be sorted
  template <typename T, typename U>
  bool contains(const std::vector<T>& container, U const& val)
  {
    auto it = std::lower_bound(
          container.begin(),
          container.end(),
          val,
          [](T const& elem, U const& val){ return elem < val; });
    
    return it != container.end() && *it == val;
  }
  
  template <typename T, typename U, typename Comp_less, typename Comp_eq>
  bool contains
  (const std::vector<T>& container, U const& val,
   Comp_less comp_less, Comp_eq comp_eq)
  {
    auto it = std::lower_bound(
          container.begin(),
          container.end(),
          val,
          comp_less);
    
    return it != container.end() && comp_eq(*it, val);
  }

  // Sign of val
  template <typename T> int sgn(T val)
  { return (T(0) < val) - (val < T(0)); }

  template <typename Object_Type, typename Return_Type = Object_Type>
  struct StoreConst
  {
    using value_type = Return_Type;
    
    Return_Type operator()() const
    { return obj; }
    const Object_Type obj;
  };

  template <typename Object_Type, typename Return_Type = Object_Type>
  struct Store
  {
    using value_type = Return_Type;
    
    Return_Type operator()() const
    { return obj; }
    Object_Type obj;
  };

  struct Empty
  {
    template <typename ...Args>
    Empty(Args...){}
    Empty(){}
  };
  
  template <typename Stream, typename Container>
  void print(Stream& stream, Container const& container, bool delimit_first = 0, std::string delimiter = "\t")
  {
    std::string delim = delimit_first ? delimiter : "";
    for (auto const& val : container)
    {
      stream << delim << val;
      delim = delimiter;
    }
  }
  
  struct DoNothing
  {
    template <typename ...Args>
    void operator()(Args...) const
    {}
    
  };

  template <typename ...Args>
  struct DoFalse
  { bool operator()(Args...) const { return false; } };

  template <typename Object, typename Params>
  Object Create(Params const& parameters = useful::Empty())
  {
    Object tracker;
    tracker.Initialize(parameters);

    return tracker;
  }
  
  template <typename Object>
  class Maker
  {
    template <typename ...Args>
    Object operator()(Args... args)
    { return Object{ args... }; }
  };

  template <typename T>
  struct Creator
  {
    typedef T value_type;
    typedef T* pointer;

    pointer operator() (void) const
    { return (new T ()); }

    pointer operator() (double param) const
    { return (new T(param)); }
  };

  template <typename Object_Type>
  struct Forward
  {
    Object_Type operator()(Object_Type const& object)
    { return object; }
  };

  template <typename Object_Type>
  struct Forward_ref
  {
    Object_Type const& operator()(Object_Type const& object)
    { return object; }
  };

  struct Bin
  {
    double left_edge;
    double right_edge;

    double mean() const
    { return (left_edge + right_edge)/2.; }

    double width() const
    { return right_edge - left_edge; }
  };
  
  template
  <typename Value, typename Container_boundaries, typename Container>
  void bin(Value const& value, Container_boundaries const& boundaries, Container& hist)
  {
    if (value < boundaries[0])
    {
      ++hist[0];
      return;
    }
    if (value >= boundaries[boundaries.size()-1])
    {
      ++hist[boundaries.size()-2];
      return;
    }
    
    auto upper_edge = std::upper_bound(std::begin(boundaries),
                                       std::end(boundaries),
                                       value);
    std::size_t idx = upper_edge - std::begin(boundaries);
    ++hist[idx-1];
  }
  
  template
  <typename Value, typename Container_boundaries, typename Container, typename Compare>
  void bin(Value const& value, Container_boundaries const& boundaries, Container const& hist, Compare comp)
  {
    if (comp(value, boundaries[0]))
    {
      ++hist[0];
      return;
    }
    if (!comp(value, boundaries[boundaries.size()-1]))
    {
      ++hist[boundaries.size()-2];
      return;
    }
    
    auto upper_edge = std::upper_bound(std::begin(boundaries),
                                       std::end(boundaries),
                                       value);
    std::size_t idx = upper_edge - std::begin(boundaries);
    ++hist[idx-1];
  }

  template <typename TT> struct Selector_t{};
  template <typename TT, TT(val)> struct Selector{};

  template <typename T>
  T& ensure_ref(T const* obj)
  { return *obj; }

  template <typename T>
  T& ensure_ref(T const& obj)
  { return obj; }

  // Simple hash for a container by combining element hashes
  template<typename Container>
  struct hash_container
  {
    std::size_t operator()(Container const& container) const
    {
      using value_type = typename Container::value_type;
      std::size_t hash_val = std::hash<value_type>()(container[ 0 ]);
      for (size_t ii = 1; ii < container.size() - 1; ++ii)
      {
        hash_val = (hash_val ^ (std::hash<value_type>()(container[ii]) << 1)) >> 1;
      }
      hash_val = hash_val ^ (std::hash<value_type>()(container[container.size()-1]) << 1);

      return hash_val;
    }
  };

  template <typename T>
  void hash_combine(std::size_t & seed, T const& v)
  {
    std::hash<T> hasher;
    seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
  }

  template <typename S, typename T>
  struct hash_pair
  {
    std::size_t operator()(const std::pair< S, T > & v) const
    {
      std::size_t seed = 0;
      hash_combine(seed, v.first);
      hash_combine(seed, v.second);
      return seed;
    }
  };

  // Replace NaNs
  template < typename T >
  void deNaN(T& to_replace, T replace_with = T())
  {
    if (to_replace != to_replace) to_replace = replace_with;
  }

  // Replace values below tolerance
  template <typename T>
  void chop(T& to_chop, T tolerance, T replace_with = T())
  {
    if (std::abs(to_chop - replace_with) < tolerance) to_chop = replace_with;
  }

  // Copy end element of container to position, then erase end element
  template <typename Container>
  void swap_erase(Container& container, std::size_t position)
  {
    if (position != container.size() - 1)
      container[position] = container.back();
    container.pop_back();
  }

  // Copy pointer to end element of vector to position, then delete end element
  template < typename Container >
  void swap_delete(Container& container, std::size_t position)
  {
    if (position != container.size() - 1)
      container[position] = container.back();
    delete container.back();
    container.pop_back();
  }
  
  std::size_t countlines(FILE *fin)
  {
    std::size_t lines = 0;
    int maxlength = 255;
    char buffer[ maxlength + 1 ];

    while (fgets(buffer , maxlength + 1 , fin) != NULL)
    {
      if (buffer[std::strlen(buffer) - 1] != '\n')
        throw "Line too long.";
      ++lines;
    }
    if (!feof(fin))
      throw "Could not reach end of file.";

    rewind(fin);

    return lines;
  }

  //	To check if a class has a method
  //	From here: https://stackoverflow.com/questions/29772601/why-is-sfinae-causing-failure-when-there-are-two-functions-with-different-signat
  //	bundle of types
  template<class...>struct types{using type=types;};
  template<class...>struct voider{using type=void;};
  //	Some dists still do not to have void_t
  template<class...Ts>using void_t=typename voider<Ts...>::type;
  //	hide the SFINAE stuff in a details namespace:
  namespace details
  {
    template<template<class...>class Z, class types, class=void>
    struct has_method : std::false_type {};
    template<template<class...>class Z, class...Ts>
    struct has_method<Z,types<Ts...>,void_t<Z<Ts...>>>:std::true_type{};
  }
  // has_method<template, types...> is true iff template<types...> is valid
  template<template<class...>class Z, class...Ts>
  using has_method=details::has_method<Z,types<Ts...>>;

  // Implemented here because some versions of gcc have trouble with the standard one
  template <typename T>
  bool isnan(T const& val)
  { return val != val; }

  template <typename Tuple, typename F, std::size_t ...Indices>
  void for_each_impl(Tuple const& tuple, F f, std::index_sequence< Indices... >)
  {
    using swallow = int[];
    (void)swallow{ 1, (f(std::get<Indices>(tuple)), void(), int{})... };
  }

  template <typename Tuple, typename F>
  void for_each(Tuple const& tuple, F f)
  {
    constexpr std::size_t N = std::tuple_size< std::remove_reference_t< Tuple > >::value;
    for_each_impl(tuple, f, std::make_index_sequence< N >{});
  }

  template <typename Tuple, typename F, std::size_t ...Indices>
  void for_each_impl(Tuple&& tuple, F&& f, std::index_sequence< Indices... >)
  {
    using swallow = int[];
    (void)swallow{ 1, (f(std::get< Indices >(std::forward< Tuple >(tuple))), void(), int{})... };
  }

  template <typename Tuple, typename F>
  void for_each(Tuple&& tuple, F&& f)
  {
    constexpr std::size_t N = std::tuple_size< std::remove_reference_t< Tuple > >::value;
    for_each_impl(std::forward<Tuple>(tuple), std::forward<F>(f),  std::make_index_sequence<N>{});
  }

  // Indices for template metamagic
  template <std::size_t... Indices>
  struct indices
  { using next = indices< Indices..., sizeof...(Indices) >; };

  template <std::size_t size>
  struct build_indices
  { using type = typename build_indices< size - 1 >::type::next; };

  template <>
  struct build_indices< 0 >
  { using type = indices<>; };
}

#endif /* USEFUL_H_ */
