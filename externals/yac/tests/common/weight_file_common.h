// Copyright (c) 2024 The YAC Authors
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef WEIGHT_FILE_COMMON_H
#define WEIGHT_FILE_COMMON_H

#ifndef YAC_CORE_H
#include "location.h"
#endif // YAC_CORE_H

/**
 * creates a weight file
 * @param[in] file_name name of the weight file
 * @param[in] src_id source indices for all links
 * @param[in] tgt_id target indices for all links
 * @param[in] weights weights for all links
 * @param[in] num_links number of links
 * @param[in] src_locations location for each source field
 * @param[in] num_src_fields number of source fields
 * @param[in] num_links_per_src_field number of links per source field
 * @param[in] tgt_id_fixed target indices for all points that receive a fixed
 *                         value
 * @param[in] num_fixed_tgt number of target points that receive a fixed value
 * @param[in] fixed_values fixed values to be used
 * @param[in] num_tgt_per_fixed_value number of target points per fixed value
 * @param[in] num_fixed_values number of fixed values
 * @param[in] tgt_location location of the target field
 * @param[in] src_grid_name name of the source grid
 * @param[in] tgt_grid_name name of the target grid
 * @remark the sum of all entries in num_links_per_src_field must be equal to
 *         num_links
 * @remark the sum of all entries in num_tgt_per_fixed_value must be equal to
 *         num_fixed_tgt
 * @remark the links have to be ordered by the source field index (ascending
 *         order)
 */
void write_weight_file(char const * file_name, int const * src_id,
                       int const * tgt_id, double const * weights,
                       unsigned num_links,
                       enum yac_location const * src_locations,
                       unsigned num_src_fields,
                       int const * num_links_per_src_field,
                       int * tgt_id_fixed, unsigned num_fixed_tgt,
                       double * fixed_values,
                       int * num_tgt_per_fixed_value,
                       unsigned num_fixed_values,
                       enum yac_location tgt_location,
                       char const * src_grid_name,
                       char const * tgt_grid_name);

/**
 * checks whether data in a weight file is identical to the reference data
 * provided
 * @param[in] file_name name of the weight file that is to be checked
 * @param[in] ref_src_address reference source indices for all links
 * @param[in] ref_tgt_address reference target indices for all links
 * @param[in] ref_weights reference weights for all links
 * @param[in] ref_num_links reference number of links
 * @param[in] ref_src_locations reference location for each source field
 * @param[in] ref_num_src_fields reference number of source fields
 * @param[in] ref_num_links_per_src_field reference number of links per source
 *                                        field
 * @param[in] ref_tgt_id_fixed reference target indices for all points that
 *                             receive a fixed value
 * @param[in] ref_fixed_values reference fixed values to be used
 * @param[in] ref_num_tgt_per_fixed_value number of target points per fixed
 *                                        value
 * @param[in] ref_num_fixed_values reference number of fixed values
 * @param[in] ref_tgt_location reference location of the target field
 * @param[in] ref_src_grid_name name of the source grid
 * @param[in] ref_tgt_grid_name name of the target grid
 * @remark routine adds the number of differences found to the global variable
 *         \ref err_count__ (see \ref tests.h)
 */
void check_weight_file(char const * file_name, int const * ref_src_address,
                       int const * ref_tgt_address, double const * ref_weights,
                       unsigned ref_num_links,
                       enum yac_location const * ref_src_locations,
                       unsigned ref_num_src_fields,
                       int const * ref_num_links_per_src_field,
                       int const * ref_tgt_id_fixed,
                       double const * ref_fixed_values,
                       int const * ref_num_tgt_per_fixed_value,
                       unsigned ref_num_fixed_values,
                       enum yac_location ref_tgt_location,
                       char const * ref_src_grid_name,
                       char const * ref_tgt_grid_name);

#endif // WEIGHT_FILE_COMMON_H

