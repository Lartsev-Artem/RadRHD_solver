#include "solvers_struct.h"

#include "graph_main.h"
#include "trace_main.h"
#include "illum_calc_gpu_async.h"

#include "global_value.h"
#include "writer_bin.h"
#include "convert_face_to_cell.h"

#include "reader_txt.h"
#include "reader_bin.h"

#include "mpi_shifts.h"

int TracingObserver()
{
    WRITE_LOG("Start TracingObserver\n");

    TracerData trace_data;
    trace_data.Init(glb_files);
        
    grid_directions_t grid_direction;    
    DIE_IF(files_sys::txt::ReadSphereDirectionCartesian(glb_files.name_file_sphere_direction, grid_direction));

    std::vector<IntId> projection; // номера ячеек проецируемых на пиксель соответствующих данному направлению
    DIE_IF(files_sys::bin::ReadSimple(glb_files.trace_address + F_RAY_TRACE + ".bin", projection));

    DIE_IF(projection.size() != grid_direction.size);
    
    int id = get_mpi_id();
    int np =  get_mpi_np();

    std::vector<IdType> disp;
    std::vector<IdType> send;
    GetDisp(np,grid_direction.size, disp);
    GetSend(np,grid_direction.size, send);
    int left  = disp[id];
    int right  = disp[id] + send[id];
    std::vector<Type> illum_projection(send[id], 0);

//    for (int i = id; i < grid_direction.size; i += np)
    for (int i = left; i < right; i ++)
    {    
        Vector3 direction = grid_direction.directions[i].dir;
        WRITE_LOG("Direction #%d: %lf %lf %lf\n",i, direction[0],direction[1],direction[2]);
    
        if(projection[i] < 0)
        {
            continue; //нет пересечения на пиксель, не выполняем здесь расчет
        }        

        graph::RunGraphModule(trace_data, direction);        
        trace::RunTracesModule(trace_data, direction);                    
        illum::separate_gpu::CalculateIllumByDirection(direction, trace_data.X0, trace_data.graph_cell_faces, trace_data.graph_bound_faces, trace_data.geo_grid);        
         
        grid_t& grid = trace_data.geo_grid;
#if 1
        DIE_IF(projection[i]/CELL_SIZE >= grid.size);
        illum_projection[i-left] = grid.Illum[projection[i]/CELL_SIZE];
#else
        std::vector<Type> illum;
        GetCellDataBySelectedDirection(grid.size, grid.size_dir, 0, grid.Illum, illum);
        files_sys::bin::WriteSimple(glb_files.solve_address + std::to_string(i) + F_ILLUM, illum);
#endif
    }

    WRITE_LOG("End TracingObserver\n");
    MPI_BARRIER(MPI_COMM_WORLD);

    files_sys::bin::WriteSimpleMPI(glb_files.trace_address + F_SOLVE + "0.bin", (int)grid_direction.size, illum_projection.data(), left, right);    
    return 0;
}