// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

// #define VERBOSE

#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <math.h>
#include "yac.h"
#include "yac_utils.h"

/* ------------------------------------------------- */

const char * fieldName[] = {"icon_out", "cube_out"};

// redefine YAC assert macros
#undef YAC_ASSERT
#define STR_USAGE "Usage: %s -c configFilename -n cube edge length\n"
#define YAC_ASSERT(exp, msg) \
  { \
    if(!((exp))) { \
      fprintf(stderr, "ERROR: %s\n" STR_USAGE, msg, argv[0]); \
      exit(EXIT_FAILURE); \
    } \
  }

static void parse_arguments(
  int argc, char ** argv, char const ** configFilename, size_t * cube_n);

int main (int argc, char *argv[]) {

  // Initialisation of MPI

  MPI_Init (NULL, NULL);

  char const * configFilename = "perf_toy.yaml"; // default configuration file
  size_t cube_n = 405; //default cube edge length
  parse_arguments(argc, argv, &configFilename, &cube_n);

  xt_initialize(MPI_COMM_WORLD);

  /* The initialisation phase includes the reading of the
   * coupling configuration */
#ifdef VERBOSE
  printf (". main: calling yac_cinit\n");
#endif

  double tic, toc, time;

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  // move yaml file to node cache
  {
    FILE *f = fopen(configFilename, "rb");
    fseek(f, 0, SEEK_END);
    long fsize = ftell(f);
    fseek(f, 0, SEEK_SET);  /* same as rewind(f); */

    char *string = malloc(fsize + 1);
    size_t dummy = fread(string, 1, fsize, f);
    (void)(dummy); //UNUSED
    fclose(f);
    free(string);
  }

  yac_cinit ();
  yac_cread_config_yaml(configFilename);

  /* The usual component definition, here for two sequential components on the
     same process */

#ifdef VERBOSE
  printf (". main: calling yac_cdef_comp\n");
#endif

  int comp_id;
  char * comp_name = "CUBE";

  yac_cdef_comp ( comp_name, &comp_id );
#ifdef VERBOSE
  printf(". main: defined %s with local comp ID %i \n", comp_name, comp_id);
#endif

  MPI_Comm local_comm;

  yac_cget_comp_comm(comp_id, &local_comm);

  int rank, size;

  MPI_Comm_rank(local_comm,&rank);
  MPI_Comm_size(local_comm,&size);

  MPI_Comm_free(&local_comm);


  int cell_point_id;
  int corner_point_id;

  int field_ids[2];
  int grid_id;

  unsigned nbr_vertices;
  unsigned nbr_cells;
  unsigned * num_vertices_per_cell;
  unsigned * cell_to_vertex;
  double * x_vertices;
  double * y_vertices;
  double * x_cells;
  double * y_cells;

  int * global_cell_id;
  int * global_cell_id_rank;
  int * global_corner_id;
  int * global_corner_id_rank;

  yac_generate_part_cube_grid_information(
    (unsigned)cube_n, &nbr_vertices, &nbr_cells,
    &num_vertices_per_cell, &cell_to_vertex,
    &x_vertices, &y_vertices, &x_cells,
    &y_cells, &global_cell_id,
    &global_cell_id_rank, &global_corner_id,
    &global_corner_id_rank, rank, size);

  double * x_points, * y_points;

  x_points = x_vertices;
  y_points = y_vertices;

  double * x_center, * y_center;

  x_center = x_cells;
  y_center = y_cells;

  yac_cdef_grid_unstruct(
    "cube_grid", nbr_vertices, nbr_cells, (int*)num_vertices_per_cell,
    x_vertices, y_vertices, (int*)cell_to_vertex, &grid_id);

  int * cell_core_mask = malloc(nbr_cells * sizeof(*cell_core_mask));
  int * corner_core_mask = malloc(nbr_vertices * sizeof(*corner_core_mask));
  for (unsigned i = 0; i < nbr_cells; ++i)
    cell_core_mask[i] = global_cell_id_rank[i] == -1;
  for (unsigned i = 0; i < nbr_vertices; ++i)
    corner_core_mask[i] = global_corner_id_rank[i] == -1;

  yac_cset_global_index(global_cell_id, YAC_LOCATION_CELL, grid_id);
  yac_cset_core_mask(cell_core_mask, YAC_LOCATION_CELL, grid_id);
  yac_cset_global_index(global_corner_id, YAC_LOCATION_CORNER, grid_id);
  yac_cset_core_mask(corner_core_mask, YAC_LOCATION_CORNER, grid_id);

  yac_cdef_points_unstruct(
    grid_id, nbr_cells, YAC_LOCATION_CELL, x_center, y_center, &cell_point_id);
  yac_cdef_points_unstruct(
    grid_id, nbr_vertices, YAC_LOCATION_CORNER, x_points, y_points, &corner_point_id);

  /* Field definition for the component */

  yac_cdef_field(
    fieldName[0], comp_id, &cell_point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
    &field_ids[0]);
  yac_cdef_field(
    fieldName[1], comp_id, &cell_point_id, 1, 1, "2", YAC_TIME_UNIT_SECOND,
    &field_ids[1]);

  toc=MPI_Wtime();
  time = toc-tic;
  printf ("CUBE: Time for initialisation %f\n", time ); 

  /* Search. */

#ifdef VERBOSE
  printf (". main: calling yac_cenddef\n");
#endif

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  yac_cenddef( );

  toc=MPI_Wtime();
  time = toc-tic;
  printf("CUBE: Time for search %f\n", time ); 

  double * conserv_in = malloc(nbr_cells * sizeof(*conserv_in));
  double * avg_in = malloc(nbr_vertices * sizeof(*avg_in));

  for (unsigned i = 0; i < nbr_cells; ++i)
    conserv_in[i] = -10;
  for (unsigned i = 0; i < nbr_vertices; ++i) avg_in[i] = -10;

  double * cell_out = malloc(nbr_cells * sizeof(*cell_out));
  double * vertex_out = malloc(nbr_vertices * sizeof(*vertex_out));

  int err;
  int info;

  for (unsigned i = 0; i < nbr_cells; ++i)
    cell_out[i] = yac_test_func(x_center[i], y_center[i]);
  for (unsigned i = 0; i < nbr_vertices; ++i)
    vertex_out[i] = yac_test_func(x_vertices[i], y_vertices[i]);

  MPI_Barrier(MPI_COMM_WORLD);

  tic=MPI_Wtime();

  {
    double *point_set_data[1];
    double **collection_data[1] = {point_set_data};

    point_set_data[0] = cell_out;
    yac_cput(field_ids[1], 1, collection_data, &info, &err);
  }

  {
    double *collection_data[1] = {conserv_in};

    yac_cget(field_ids[0], 1, collection_data, &info, &err);
  }

  toc=MPI_Wtime();
  time = toc-tic;
  printf ("CUBE: Time for ping-pong %f\n", time );

#ifdef VTK_OUTPUT
  //----------------------------------------------------------
  // write field to vtk output file
  //----------------------------------------------------------

  char vtk_filename[32];

  sprintf(vtk_filename, "perf_toy_cube_%d.vtk", rank);

  double point_data[nbr_vertices][3];
  for (int i = 0; i < nbr_vertices; ++i) {
   LLtoXYZ(x_vertices[i], y_vertices[i], point_data[i]);
  }

  YAC_VTK_FILE *vtk_file = yac_vtk_open(vtk_filename, "cube_out");
  yac_vtk_write_point_data(vtk_file, (double *)point_data, nbr_vertices);
  yac_vtk_write_cell_data(vtk_file, (unsigned *)cell_to_vertex,
                      (unsigned*)num_vertices_per_cell, nbr_cells);
  yac_vtk_write_point_scalars_int(
    vtk_file, corner_core_mask, nbr_vertices, "corner_core_mask");
  yac_vtk_write_point_scalars_int(
    vtk_file, global_corner_id, nbr_vertices, "global_corner_id");
  yac_vtk_write_cell_scalars_int(
    vtk_file, cell_core_mask, nbr_cells, "cell_core_mask");
  yac_vtk_write_cell_scalars_int(
    vtk_file, global_cell_id, nbr_cells, "global_cell_id");

  yac_vtk_write_cell_scalars_double(
    vtk_file, conserv_in, nbr_cells, "conserv_in");
  yac_vtk_write_cell_scalars_double(
    vtk_file, cell_out, nbr_cells, "cell_out");
  yac_vtk_write_point_scalars_double(
    vtk_file, avg_in, nbr_vertices, "avg_in");
  yac_vtk_write_point_scalars_double(
    vtk_file, vertex_out, nbr_vertices, "vertex_out");

  yac_vtk_close(vtk_file);
#endif // VTK_OUTPUT

  free(corner_core_mask);
  free(global_corner_id_rank);
  free(global_corner_id);
  free(cell_core_mask);
  free(global_cell_id_rank);
  free(global_cell_id);
  free(x_cells);
  free(y_cells);
  free(x_vertices);
  free(y_vertices);
  free(num_vertices_per_cell);
  free(cell_to_vertex);

  yac_cfinalize();

  xt_finalize();

  MPI_Finalize();

  free(conserv_in);
  free(cell_out);
  free(avg_in);
  free(vertex_out);

  return EXIT_SUCCESS;
}

static void parse_arguments(
  int argc, char ** argv,
  char const ** configFilename, size_t * cube_n) {

  int opt;
  while ((opt = getopt(argc, argv, "c:n:")) != -1) {
    YAC_ASSERT((opt == 'c') || (opt == 'n'), "invalid command argument")
    switch (opt) {
      default:
      case 'c':
        *configFilename = optarg;
        break;
      case 'n':
        *cube_n = atoi(optarg);
        YAC_ASSERT(*cube_n > 0, "invalid cube edge length");
        break;
    }
  }
}
