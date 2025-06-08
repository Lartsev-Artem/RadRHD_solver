/**
 * @file gen_illum_by_dirs.cpp
 * @author Artem
 * @brief Генерация излучения на сфере направлений по закону (1/r)
 * @version 0.1
 * @date 2025-06-07
 *
 * @copyright Copyright (c) 2025
 *
 */

#include "reader_json.h"
#include "reader_txt.h"
#include "writer_bin.h"

int main(int argc, char *argv[])
{
    MPI_START(argc, argv);

    std::string file_config = "config/directories_cfg.json";

    if (argc != 4)
    {
        printf("Error input data!\n");
        printf("Input:\n");
        printf("path\\cfg.json\n");
        printf("path\\output_file.bin\n");
        printf("path\\ref_direction\n");
        return e_completion_fail;
    }

    file_config = argv[1];
    std::string out_file = argv[2];
    int ref_id = std::stoi(argv[3]);

    if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
        return e_completion_fail;

    grid_directions_t grid;
    files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid);

    Vector3 dir0 = grid.directions[ref_id].dir;
    std::vector<Type> vals(grid.size, 0);
    for (size_t i = 0; i < grid.size; i++)
    {
        if (ref_id != i)
        {
            vals[i] = 1.0 / (grid.directions[i].dir - dir0).norm();
        }
    }
    vals[ref_id] = 1.1 * (*std::max_element(vals.begin(), vals.end()));

    MPI_END;
    return files_sys::bin::WriteSimple(out_file, vals);
}