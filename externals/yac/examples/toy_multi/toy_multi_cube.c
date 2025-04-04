// Copyright  (C)  2013 DKRZ, MPI-M
// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

// #define VERBOSE

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "yac.h"
#include "yac_utils.h"

#include "toy_multi_common.h"

// redefine YAC assert macros
#undef YAC_ASSERT
#define STR_USAGE "Usage: %s -n cube edge length\n"
#define YAC_ASSERT(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments( int argc, char ** argv, size_t * cube_n);

int main (int argc, char *argv[]) {

  double tic;

  // Initialisation of MPI

  MPI_Init (0, NULL);

  MPI_Barrier(MPI_COMM_WORLD);
  tic=MPI_Wtime();

  size_t cube_n = 50; // default cube edge length
  parse_arguments(argc, argv, &cube_n);
  yac_cinit();
  yac_cdef_calendar(YAC_PROLEPTIC_GREGORIAN);
  yac_cdef_datetime("2008-03-09T16:05:07", "2008-03-10T16:05:07");

  int comp_id;
  char * comp_name = "CUBE";

  yac_cdef_comp(comp_name, &comp_id);

  MPI_Comm comp_comm;
  yac_cget_comp_comm(comp_id, &comp_comm);

  int rank, comp_rank, comp_size;
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_rank(comp_comm, &comp_rank);
  MPI_Comm_size(comp_comm, &comp_size);
  MPI_Comm_free(&comp_comm);

  unsigned nbr_vertices;
  unsigned nbr_cells;
  unsigned * num_vertices_per_cell;
  unsigned * cell_to_vertex;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * cell_core_mask;
  int * corner_core_mask;
  int * global_cell_id;
  int * global_corner_id;

  yac_generate_part_cube_grid_information(cube_n, &nbr_vertices, &nbr_cells,
                                          &num_vertices_per_cell, &cell_to_vertex,
                                          &x_vertices, &y_vertices, &x_cells,
                                          &y_cells, &global_cell_id,
                                          &cell_core_mask, &global_corner_id,
                                          &corner_core_mask, comp_rank, comp_size);

  int grid_id;
  char grid_name[256];
  sprintf(grid_name, "%s_grid", comp_name);

  yac_cdef_grid_unstruct(
    grid_name, nbr_vertices, nbr_cells, (int*)num_vertices_per_cell,
    x_vertices, y_vertices, (int*)cell_to_vertex, &grid_id );

  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);
  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);

  int cell_point_id;
  int corner_point_id;

  yac_cdef_points_unstruct(
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_cells, y_cells, &cell_point_id);
  yac_cdef_points_unstruct(
    grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_vertices, y_vertices,
    &corner_point_id);

  double * cell_point_data =
    malloc(nbr_cells * sizeof(*cell_point_data));
  double * corner_point_data =
    malloc(nbr_vertices * sizeof(*corner_point_data));

  for (unsigned i = 0; i < nbr_cells; ++i)
    cell_point_data[i] = yac_test_func(x_cells[i], y_cells[i]);

  for (unsigned i = 0; i < nbr_vertices; ++i)
    corner_point_data[i] = yac_test_func(x_vertices[i], y_vertices[i]);

  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0)
    printf(
      "toy_multi_common(%s): Time for initialisation %f\n",
      comp_name, MPI_Wtime()-tic);

  // initialize VTK file
  char vtk_filename[32];
  sprintf(vtk_filename, "%s_out_%d.vtk", comp_name, comp_rank);

  yac_coordinate_pointer point_data =
    malloc(nbr_vertices * sizeof(*point_data));
  for (unsigned i = 0; i < nbr_vertices; ++i)
    LLtoXYZ(x_vertices[i], y_vertices[i], point_data[i]);

  YAC_VTK_FILE *vtk_file = yac_vtk_open(vtk_filename, comp_name);
  yac_vtk_write_point_data(vtk_file, (double *)point_data, nbr_vertices);
  yac_vtk_write_cell_data(
    vtk_file, (unsigned *)cell_to_vertex,
    (unsigned*)num_vertices_per_cell, nbr_cells);
  yac_vtk_write_point_scalars_int(
    vtk_file, corner_core_mask, nbr_vertices, "corner_core_mask");
  yac_vtk_write_point_scalars_int(
    vtk_file, global_corner_id, nbr_vertices, "global_corner_id");
  yac_vtk_write_cell_scalars_int(
    vtk_file, cell_core_mask, nbr_cells, "cell_core_mask");
  yac_vtk_write_cell_scalars_int(
    vtk_file, global_cell_id, nbr_cells, "global_cell_id");

  int ret =
    run_toy_multi_common(
      comp_name, comp_id, grid_id, cell_point_id, corner_point_id,
      cell_point_data, corner_point_data, vtk_file);

  free(point_data);
  free(corner_point_data);
  free(cell_point_data);
  free(cell_core_mask);
  free(corner_core_mask);
  free(global_cell_id);
  free(global_corner_id);
  free(x_vertices);
  free(y_vertices);
  free(x_cells);
  free(y_cells);
  free(num_vertices_per_cell);
  free(cell_to_vertex);

  return ret;
}

static void parse_arguments(int argc, char ** argv, size_t * cube_n) {

  int opt;
  while ((opt = getopt(argc, argv, "n:")) != -1) {
    YAC_ASSERT((opt == 'c') || (opt == 'n'), "invalid command argument")
    switch (opt) {
      default:
      case 'n':
        *cube_n = atoi(optarg);
        YAC_ASSERT(*cube_n > 0, "invalid cube edge length");
        break;
    }
  }
}
