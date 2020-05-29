/*
 * ScalarField.h
 *
 *  Created on: Sep 10, 2014
 *      Author: tomas
 */

#ifndef SCALARFIELD_H_
#define SCALARFIELD_H_

#include <cstdio>

namespace field
{
  template <typename Grid, typename Value_t>
  class ScalarField
  {
  public:
    Grid const& grid;
    using value_type = Value_t;
    using container_type = std::vector<value_type>;
    using position_type = typename Grid::position_type;
    using index_type = typename Grid::index_type;
    container_type field{ container_type(grid.size()) };
    
    ScalarField(Grid const& grid, value_type value = {})
    : grid{ grid }
    , field(grid.size(), value)
    {}
    
    template <typename Container_idx>
    ScalarField(Grid const& grid, value_type value, Container_idx const& index)
    : grid{ grid }
    { set(value, index); }
    
    ScalarField(Grid const& grid, container_type const& field)
    : grid{ grid }
    , field{ field }
    {}
    
    template <typename Container_idx>
    ScalarField(Grid const& grid, container_type const& field, Container_idx const& index)
    : ScalarField{ grid }
    { set(field, index); }
    
    template <typename Value_center>
    ScalarField(Grid const& grid, Value_center const& value_center)
    : ScalarField{ grid }
    { set(value_center); }
    
    template <typename Value_center, typename Container_idx>
    ScalarField(Grid const& grid, Value_center const& value_center, Container_idx const& index)
    : ScalarField{ grid }
    { set(value_center, index); }
    
    void set(container_type const& field)
    {
      for (std::size_t idx = 0; idx < field.size(); ++idx)
        this->field[idx] = field[idx];
    }
    
    template <typename Container_idx>
    void set(container_type const& field, Container_idx const& index)
    {
      for (std::size_t idx = 0; idx < index.size(); ++idx)
        this->field[index[idx]] = field[idx];
    }
    
    template <typename Value_center>
    void set(Value_center value_center)
    {
      for (std::size_t idx = 0; idx < field.size(); ++idx)
        field[idx] = value_center(grid.cell_centers(idx));
    }
    
    template <typename Value_center, typename Container_idx>
    void set(Value_center value_center, Container_idx const& index)
    {
      for (std::size_t idx = 0; idx < index.size(); ++idx)
        field[index[idx]] = value_center(grid.cell_center(index[idx]));
    }
    
    void set(value_type value)
    {
      for (std::size_t idx = 0; idx < grid.size(); ++idx)
        field[idx] = value;
    }
    
    template <typename Container_idx>
    void set(value_type value, Container_idx const& index)
    {
      for (std::size_t idx = 0; idx < index.size(); ++idx)
        field[index[idx]] = value;
    }
    
    value_type operator()(std::size_t idx) const
    { return field[idx]; }
    
    template <typename Position>
    value_type operator()(Position const& position) const
    { return field[grid.cell(position)]; }
    
    value_type operator()(index_type const& indexes) const
    { return field[grid.cell(indexes)]; }
    
    value_type const& operator[] (std::vector<std::size_t> const& indexes) const
    { return field[grid.cell(indexes)]; }
    
    value_type& operator[] (std::vector<std::size_t> const& indexes)
    { return field[grid.cell(indexes)]; }
    
    value_type const& operator[] (std::size_t index) const
    { return field[index]; }
    
    value_type& operator[] (std::size_t index)
    { return field[index]; }
    
    value_type value_cell_up
    (index_type const& index, std::size_t dim, std::size_t offset = 1) const
    {
      return field[grid.neighbor_cell_up(index, dim, offset)];
    }

    value_type value_cell_down
    (index_type const& index, std::size_t dim, std::size_t offset = 1) const
    {
      return field[grid.neighbor_cell_down(index, dim, offset)];
    }
    
    value_type value_cell_up
    (position_type const& position, std::size_t dim, std::size_t offset = 1) const
    {
      return field[grid.neighbor_cell_up(position, dim, offset)];
    }

    value_type value_cell_down
    (position_type const& position, std::size_t dim, std::size_t offset = 1) const
    {
      return field[grid.neighbor_cell_down(position, dim, offset)];
    }
    
    value_type value_cell_up
    (std::size_t idx, std::size_t dim, std::size_t offset = 1) const
    {
      return field[grid.neighbor_cell_up(idx, dim, offset)];
    }

    value_type value_cell_down
    (std::size_t idx, std::size_t dim, std::size_t offset = 1) const
    {
      return field[grid.neighbor_cell_down(idx, dim, offset)];
    }
  };
  
  template <typename Concentration, typename Grid, typename BC_Container>
  class Gradient
  {
  public:
    using value_type = std::vector<double>;
    
    Concentration const& concentration;
    Grid const& grid;
    BC_Container const& bc_container;
    
    Gradient
    (Concentration const& concentration, Grid const& grid,
     BC_Container const& bc_container)
    : concentration{ concentration }
    , grid{ grid }
    , bc_container{ bc_container }
    {}
    
    value_type operator()(std::size_t idx) const
    {
      value_type gradient;
      gradient.reserve(grid.dim());
      
      for (std::size_t dd = 0; dd < grid.dim(); ++dd)
      {
        switch (bc_container[idx].cell_type[dd])
        {
          case 0:
            gradient.push_back((concentration.value_cell_up(idx, dd) -
                                concentration.value_cell_down(idx, dd))/
                               (2.*grid.lengths[dd]));
            break;
          case 1:
            gradient.push_back((concentration.value_cell_up(idx, dd) -
                                concentration(idx))/
                               grid.lengths[dd]);
            break;
          case 2:
            gradient.push_back(
                               (concentration(idx) -
                                concentration.value_cell_down(idx, dd))/
                               grid.lengths[dd]);
            break;
          case 3:
            gradient.push_back(0.);
            break;
          default:
            throw std::runtime_error("Undefined void type.");
        }
      }
      
      return gradient;
    }
  };
}



#endif /* SCALARFIELD_H_ */
