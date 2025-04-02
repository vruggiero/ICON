// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#include <unistd.h>

#include "tests.h"
#include "test_common.h"
#include "grid2vtk.h"
#include "io_utils.h"

int main(void) {

  {
    double lon_vertices[] =
      {-10,-9,-8,-7,-6,-5,-4,-3,-2,-1,0,1,2,3,4,5,6,7,8,9,10};
    double lat_vertices[] =
      {-5,-4,-3,-2,-1,0,1,2,3,4,5};
    size_t nbr_vertices[2] = {sizeof(lon_vertices) / sizeof(lon_vertices[0]),
                              sizeof(lat_vertices) / sizeof(lat_vertices[0])};
    int cyclic[2] = {0, 0};
    struct yac_basic_grid_data reg2d_grid =
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);

    yac_write_basic_grid_data_to_file(&reg2d_grid, "test_grid2vtk_1");

    yac_basic_grid_data_free(reg2d_grid);

    if (!yac_file_exists("test_grid2vtk_1.vtk")) {
      PUT_ERR("ERROR in yac_write_basic_grid_data_to_file");
    } else {
      unlink("test_grid2vtk_1.vtk");
    }
  }

  {
    double lon_vertices[] =
      {0,20,40,60,80,100,120,140,160,180,200,220,240,260,280,300,320,340,360};
    double lat_vertices[] =
      {-90,-70,-50,-30,-10,10,30,50,70,90};
    size_t nbr_vertices[2] = {sizeof(lon_vertices) / sizeof(lon_vertices[0]),
                              sizeof(lat_vertices) / sizeof(lat_vertices[0])};
    int cyclic[2] = {0, 0};
    struct yac_basic_grid_data reg2d_grid =
      yac_generate_basic_grid_data_reg_2d_deg(
        nbr_vertices, cyclic, lon_vertices, lat_vertices);

    yac_write_basic_grid_data_to_file(&reg2d_grid, "test_grid2vtk_2");

    if (!yac_file_exists("test_grid2vtk_2.vtk")) {
      PUT_ERR("ERROR in yac_write_basic_grid_data_to_file");
    } else {
      unlink("test_grid2vtk_2.vtk");
    }

    // add dummy cell core masks
    reg2d_grid.core_cell_mask =
      xmalloc(reg2d_grid.num_cells * sizeof(*reg2d_grid.core_cell_mask));
    for (size_t i = 0; i < reg2d_grid.num_cells; ++i)
      reg2d_grid.core_cell_mask[i] = (i & 1) == 1;

    yac_write_basic_grid_data_to_file(&reg2d_grid, "test_grid2vtk_3");

    if (!yac_file_exists("test_grid2vtk_3.vtk")) {
      PUT_ERR("ERROR in yac_write_basic_grid_data_to_file");
    } else {
      unlink("test_grid2vtk_3.vtk");
    }

    yac_basic_grid_data_free(reg2d_grid);
  }

  {
    double icon_cell_coords[2][3] = {{0,120,240}, {85,85,85}};
    enum yac_edge_type icon_cell_edge_types[3] =
      {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};
    double reg2d_cell_coords[2][4] = {{-30,30,30,-30},{87,87,84,84}};
    enum yac_edge_type reg2d_cell_edge_types[4] =
      {YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE, YAC_LAT_CIRCLE_EDGE, YAC_LON_CIRCLE_EDGE};
    double fesom_cell_coords[2][6] =
      {{65,125,185,245,305,365}, {86,86,86,86,86,86}};
    enum yac_edge_type fesom_cell_edge_types[6] =
      {YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE,
       YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE, YAC_GREAT_CIRCLE_EDGE};
    struct yac_grid_cell cells[] = {
      generate_cell_deg(
        icon_cell_coords[0], icon_cell_coords[1], icon_cell_edge_types, 3),
      generate_cell_deg(
        reg2d_cell_coords[0], reg2d_cell_coords[1], reg2d_cell_edge_types, 4),
      generate_cell_deg(
        fesom_cell_coords[0], fesom_cell_coords[1], fesom_cell_edge_types, 6)};

    yac_write_grid_cells_to_file(
      cells, sizeof(cells) / sizeof(cells[0]), "test_grid2vtk_4", 5);

    if (!yac_file_exists("test_grid2vtk_4.vtk")) {
      PUT_ERR("ERROR in yac_write_grid_cells_to_file");
    } else {
      unlink("test_grid2vtk_4.vtk");
    }

    for (size_t i = 0; i < sizeof(cells)/sizeof(cells[0]); ++i)
      yac_free_grid_cell(&cells[i]);
  }

  return TEST_EXIT_CODE;
}

