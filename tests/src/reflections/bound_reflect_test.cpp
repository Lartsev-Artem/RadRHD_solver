/**
 * @file bound_reflect_test.cpp
 * @author Artem
 * @brief Вычисление отражения по всей сфере направлений для одной ячейки сетки
 * @version 0.1
 * @date 2025-06-08
 *
 * @copyright Copyright (c) 2025
 *
 */
#include "reader_json.h"
#include "reader_bin.h"
#include "reader_txt.h"
#include "writer_bin.h"

#include "intersections.h"
#include "linear_alg.h"

const Vector3 k_bound_face_norm = Vector3(0, 0, 1); ///< внешняя нормаль пространственной грани

/// \brief Прямое отражение
static void direct_reflection(
    const grid_directions_t &grid_direction,
    const std::vector<Face> &surface,
    const std::vector<Type> &illum_src,
    std::vector<Type> &reflection_val)
{
    reflection_val.assign(grid_direction.size, 0);

    for (size_t i = 0; i < grid_direction.size; i++)
    {
        const Vector3 &dir = grid_direction.directions[i].dir;

        if (k_bound_face_norm.dot(dir) < 0)
        {
            continue; // пропускаем только отражаемые лучи
        }

        Vector3 refl_dir = reflection_dir(dir, k_bound_face_norm);

        for (size_t face = 0; face < grid_direction.size; face++)
        {
            Vector3 intersection;
            Type dist = intersection::RayIntersectsTriangle(Vector3::Zero(), refl_dir, surface[face], intersection);
            if (dist > 0)
            {
                reflection_val[face] = illum_src[i];
                break;
            }
        }
    }

    return;
}

/// \brief Обратное отражение
static void inverse_reflection(
    const grid_directions_t &grid_direction,
    const std::vector<Face> &surface,
    const std::vector<Type> &illum_src,
    std::vector<Type> &reflection_val)
{
    reflection_val.assign(grid_direction.size, 0);

    for (size_t i = 0; i < grid_direction.size; i++)
    {
        const Vector3 &dir = grid_direction.directions[i].dir;

        /// \note т.к. трассируем в обратном направлении, знак меняем
        if (k_bound_face_norm.dot(dir) > 0)
        {
            continue; // пропускаем только отражаемые лучи (изнутри)
        }

        // это на самом деле исходное направление
        Vector3 refl_dir = reflection_dir(dir, k_bound_face_norm);

        for (size_t face = 0; face < grid_direction.size; face++)
        {
            Vector3 intersection;
            Type dist = intersection::RayIntersectsTriangle(Vector3::Zero(), refl_dir, surface[face], intersection);
            if (dist > 0)
            {
                reflection_val[i] = illum_src[face];
                // reflection_id[i] = face; /// \note для основного расчета нужна ячейка
                break;
            }
        }
    }

    return;
}

static int read_mesh_direction(const std::string &base_address, interpolation_mesh_direction_t &grid_st)
{
    if (files_sys::bin::ReadSimple(base_address + F_NUM_NEIGHS_SURFACE, grid_st.num_of_neighs))
        return e_completion_fail;

    if (files_sys::bin::ReadSimple(base_address + F_NODE_CELLS_ID_SURFACE, grid_st.node_cells_id))
        return e_completion_fail;

    if (files_sys::bin::ReadSimple(base_address + F_NODES_COORD_SURFACE, grid_st.nodes_coord))
        return e_completion_fail;

    return e_completion_success;
}

static int init_reflection(
    const grid_directions_t &grid_direction,
    const std::vector<Face> &surface,
    std::vector<int> &reflection_id,
    grid_directions_t &new_grid_direction)
{

    reflection_id.assign(grid_direction.size, -1);
    new_grid_direction.mesh = new interpolation_mesh_direction_t;
    if (read_mesh_direction(glb_files.base_address, *new_grid_direction.mesh))
    {
        printf("mesh wasn't read!\n");
        return e_completion_fail;
    }

    for (size_t i = 0; i < grid_direction.size; i++)
    {
        const Vector3 &dir = grid_direction.directions[i].dir;

        /// \note т.к. трассируем в обратном направлении, знак меняем
        if (k_bound_face_norm.dot(dir) > 0)
        {
            continue; // пропускаем только отражаемые лучи (изнутри)
        }

        // это на самом деле исходное направление
        Vector3 refl_dir = reflection_dir(dir, k_bound_face_norm);

        for (size_t face = 0; face < grid_direction.size; face++)
        {
            Vector3 intersection;
            Type dist = intersection::RayIntersectsTriangle(Vector3::Zero(), refl_dir, surface[face], intersection);
            if (dist > 0)
            {
                reflection_id[i] = face; /// \note для основного расчета нужна ячейка
                new_grid_direction.directions[i].dir = refl_dir;
                new_grid_direction.directions[i].area =
                    grid_direction.directions[face].area / grid_direction.directions[i].area;
                break;
            }
        }
    }

    return e_completion_success;
}

/**
 * @brief Обратное отражение + интерполяция
 *
 * @param grid_direction
 * @param new_grid_direction
 * @param illum_src
 * @param reflection_id
 * @param mode - режим интерполяции (0- cellTopoint, 1-cellTopoint + linear_interpolation)
 * @param reflection_val
 */
static void inverse_inter_reflection(
    const grid_directions_t &grid_direction,
    const grid_directions_t &new_grid_direction,
    const std::vector<Type> &illum_src,
    const std::vector<int> &reflection_id,
    const int mode,
    std::vector<Type> &reflection_val)
{
    // начинаем интерполяцию
    reflection_val.assign(grid_direction.size, -1);
    const int N = grid_direction.size;
    const interpolation_mesh_direction_t *st = new_grid_direction.mesh;

    for (size_t i = 0; i < N; i++)
    {
        int id = reflection_id[i];
        DIE_IF(id >= N);
        if (id < 0)
        {
            continue;
        }
        // вершины треугольника
        Matrix3 vertex;
        vertex.row(0) = Cartesian2Spherical(st->nodes_coord[3 * id]);
        vertex.row(1) = Cartesian2Spherical(st->nodes_coord[3 * id + 1]);
        vertex.row(2) = Cartesian2Spherical(st->nodes_coord[3 * id + 2]);

        Vector3 x = Cartesian2Spherical(new_grid_direction.directions[i].dir); // точка куда отразился луч внутри треугольника

        Vector3 rhs;
        for (size_t k = 0; k < 3; k++)
        {
            const int *id_neigh_ptr = st->node_cells_id.data();
            Type sum = 0;
            int count = 1;
            if (3 * id + k == 0)
            {
                count = st->num_of_neighs[3 * id + k];
            }
            else
            {
                count = st->num_of_neighs[3 * id + k] - st->num_of_neighs[3 * id + k - 1];
                id_neigh_ptr += st->num_of_neighs[3 * id + k - 1];
            }

            for (size_t j = 0; j < count; j++)
            {
                DIE_IF(id_neigh_ptr - st->node_cells_id.data() >= st->node_cells_id.size());
                DIE_IF_ACTION(*id_neigh_ptr >= illum_src.size(),
                              WRITE_LOG("id =%d ", id_neigh_ptr););

                sum += illum_src[*id_neigh_ptr];
                id_neigh_ptr++;
            }
            rhs[k] = sum / count;
        }

        switch (mode)
        {
        case 0:
        {
            reflection_val[i] = (rhs[0] + rhs[1] + rhs[2]) / 3;
        }
        break;
        case 1:
        {
            Vector3 coef = vertex.inverse() * rhs;
            reflection_val[i] = coef[0] * x[0] + coef[1] * x[1] + coef[2];
        }
        break;

        default:
            D_LD;
        }

        //  Type dS = (grid_direction.directions[id].area / grid_direction.directions[i].area);
        // printf("%lf\n", new_grid_direction.directions[i].area);
        // reflection_val[i] *= new_grid_direction.directions[i].area;
    }
}

/// @brief Обратное отражение + расширение шаблона (ячейка + 3 соседа)
static void inverse_avg_cells_reflection(
    const grid_directions_t &grid_direction,
    const std::vector<int> &reflection_id,
    const std::vector<Type> &illum_src,
    const std::vector<IntId> &neighbors,
    std::vector<Type> &reflection_val)
{
    reflection_val.assign(grid_direction.size, 0);
    const int N = grid_direction.size;

    for (int i = 0; i < N; i++)
    {
        int id = reflection_id[i];
        DIE_IF(id >= N);
        if (id < 0)
        {
            continue;
        }

        reflection_val[i] = illum_src[id];
        for (int k = 0; k < 3; k++)
        {
            int idx = neighbors[3 * id + k] / 3;
            reflection_val[i] += illum_src[idx];
        }
        reflection_val[i] /= 4;
    }

    return;
}

int main(int argc, char *argv[])
{
    MPI_START(argc, argv);

    std::string file_config = "config/directories_cfg.json";
    if (argc < 4)
    {
        PRINT_LOG_ERR("Not enough arguments.\n"
                      "Input: \n"
                      "cfg.json\n"
                      "illum_src on sphere.bin\n"
                      "output_reflections_value.bin\n"
                      "reflection mode (int)\n");
        return e_completion_fail;
    }

    file_config = argv[1];
    std::string file_illum = argv[2];
    std::string file_output = argv[3];
    const int mode = std::stoi(argv[4]);

    if (files_sys::json::ReadStartSettings(file_config, glb_files, &_solve_mode, &_hllc_cfg))
        return e_completion_fail;

    std::vector<Type> illum;
    grid_directions_t grid_direction;
    std::vector<Face> surface;
    std::vector<IntId> neighbors;

    uint32_t err = files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid_direction);
    err |= files_sys::bin::ReadSimple(file_illum, illum);
    err |= files_sys::bin::ReadSimple(glb_files.base_address + F_SURFACE_SPHERE_DIRECTION, surface);
    err |= files_sys::bin::ReadSimple(glb_files.base_address + F_NEIGHBOR_SURFACE, neighbors);

    std::vector<int> reflection_id;
    grid_directions_t grid_reflection(grid_direction.size);

    err |= init_reflection(grid_direction, surface, reflection_id, grid_reflection);

    if (err)
    {
        RETURN_ERR("Error reading \n");
    }

    std::vector<Type> reflection_val;

    switch (mode)
    {
    case 0:
        direct_reflection(grid_direction, surface, illum, reflection_val);
        break;
    case 1:
        inverse_reflection(grid_direction, surface, illum, reflection_val);
        break;
    case 2:
        inverse_inter_reflection(grid_direction, grid_reflection, illum, reflection_id, 0, reflection_val);
        break;
    case 3:
        inverse_inter_reflection(grid_direction, grid_reflection, illum, reflection_id, 1, reflection_val);
        break;
    case 4:
        inverse_avg_cells_reflection(grid_direction, reflection_id, illum, neighbors, reflection_val);
        break;

    default:
        break;
    }

    Type sum1 = 0;
    Type sum2 = 0;
    for (size_t i = 0; i < grid_direction.size; i++)
    {
        sum1 += illum[i] * grid_direction.directions[i].area;
        sum2 += reflection_val[i] * grid_reflection.directions[i].area * grid_direction.directions[i].area;
        // sum2 += reflection_val[i] * grid_reflection.directions[i].area;
    }
    sum1 /= grid_direction.full_area;
    sum2 /= grid_direction.full_area;
    printf("Integrate value: %lf => %lf\n", sum1, sum2);

    files_sys::bin::WriteSimple(file_output, reflection_val);

    MPI_END;
    return e_completion_success;
}