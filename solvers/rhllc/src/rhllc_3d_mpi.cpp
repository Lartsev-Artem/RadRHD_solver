#if 0 // defined RHLLC && defined SOLVERS && USE_MPI
#include "rhllc_mpi.h"

#include "rhllc_flux_stab.h"
#include "rhllc_utils.h"

using namespace rhllc;

void rhllc_mpi::Hllc3dStab(const Type tau, grid_t &grid) {
  rhllc::max_signal_speed = 0;

#pragma omp parallel default(none) firstprivate(tau) shared(grid, glb_files, rhllc::max_signal_speed)
  {
    const int size_grid = grid.size;
    Type max_speed = 0;
    flux_all_t bound_val;
    const int size_face = grid.faces.size();
// потоки
#pragma omp for
    for (int i = 0; i < size_face; i++) {
      face_t &f = grid.faces[i];
      BoundConditions(f, grid.cells, bound_val);
      max_speed = std::max(max_speed, GetFluxStab(grid.cells[f.geo.id_l].conv_val, bound_val.conv_val, grid.cells[f.geo.id_l].phys_val, bound_val.phys_val, f));
    }

#pragma omp for
    for (int i = 0; i < size_grid; i++) {
      elem_t &el = grid.cells[i];
      flux_t sumF;
      for (int j = 0; j < CELL_SIZE; j++) {
        if (el.geo.sign_n[j]) {
          sumF += grid.faces[el.geo.id_faces[j]].f;
        } else {
          sumF -= grid.faces[el.geo.id_faces[j]].f;
        }
      }
      sumF *= (tau / el.geo.V);
      el.conv_val -= sumF;

      if (GetPhysValueStab(el.conv_val, el.phys_val)) {
        DIE_IF(PhysPressureFix(el.conv_val, el.phys_val));
      }
    }

    if (max_speed > rhllc::max_signal_speed) {
#pragma omp critical
      {
        rhllc::max_signal_speed = std::max(max_signal_speed, max_speed);
      }
    }

  } // omp
}

#endif

#if 0

//#ifdef RHLLC_MPI
namespace mpi_rhllc
{
    std::vector<int> send_cells;                    ///< кол-во отправок для каждого процесса
    std::vector<int> disp_cells;                    ///< смещения по ячейкам для процессов
    std::vector<IntId> id_irregular_faces;          ///< номера граней с границе областей процессов
    std::vector<MPI_request> requests_cast_phys;    ///< запросы на передачу физ. переменных
    std::vector<MPI_request> requests_send_faces;   ///< запросы на передачу потоков грани
    std::vector<MPI_request> requests_rcv_faces;    ///< запросы на приём потоков грани
}

void Init()
{
    int np=3;
    int myid =0;
    int N = grid.size;

    MpiSend(np,N,send_cells);    
    MpiDisp(np,N,disp_cells);
    grid.loc_size = send_cells[myid];//это приведет к правкам на видеокарте(возможно это уже учтено. Надо проверить)
    
    std::vector<int> metis_id={0,0,0,  1,1,1, 2,2,2}; //ReadFileConfig

    for(int id_cell = 0;id_cell<grid.size;id_cell++ )
    {
        grid.cells[id_cell].node = metis_id[id_cell];
    }

    id_irregular_faces.reserve(N/np);
    for(int id_faces = 0;id_faces<grid.size_face;id_faces++ )
    {
        grid.faces[id_faces].is_reg = 0;
        grid.faces[id_faces].node_l = metis_id[idl];
        grid.faces[id_faces].node_r = metis_id[idr];
        if( (metis_id[idr] ==  metis_id[idl]) && ( metis_id[idl] == myid))
        {
            grid.faces[id_faces].is_reg = 1;   //полностью у нас
        }
        else
        {            
            if(metis_id[idl] == myid)
            {
                MPI_request rq_send;
                MPI_request rq_rcv;

                MPI_send_init(&grid.cell[idl],1, myid->metis_id[idr]  &rq_send);
                MPI_rcv_init(&grid.cell[idr],1, metis_id[idr] ->myid  &rq_rcv);
                            
                requests_send_faces.push_back(rq_send);
                requests_rcv_faces.push_back(rq_rcv);
                id_irregular_faces.push_back(id_faces); // только на нашем узле
            }

            if(metis_id[idr] == myid)
            {
                MPI_request rq_send;
                MPI_request rq_rcv;

                MPI_send_init(&grid.cell[idr],1, myid->metis_id[idl]  &rq_send);
                MPI_rcv_init(&grid.cell[idl],1, metis_id[idl] ->myid  &rq_rcv);
                            
                requests_send_faces.push_back(rq_send);
                requests_rcv_faces.push_back(rq_rcv);
                id_irregular_faces.push_back(id_faces); // только на нашем узле
            }         
        }
    }
    id_irregular_faces.shrink_to_fit();
     

     requests_cast_phys.resize(np);
}

void hllc()
{
    int N=9;
    int np=3;
    int metis_id[9]={0,0,0,  1,1,1, 2,2,2};
    int id_irreg_faces[100];

    // (if != 0) wait_send_all(); //завершаем предыдущие вызовы
    // start_all sendrcv
    

    for(auto& f: face)
    {
        if(f.is_reg) // (f.idr==f.idl && f_idr==myid)
        {
            //calc
        }
    }

    //wait_rcv_all

    for(auto id: id_irreg_faces)
    {
        //calc: faces[id]
    }

    for(auto c: cells)
    {
        if(c.node == myid)
        {
            //calc
        }
    }
    
    for(int id=0;id<np;id++)
    {
        IBcast(cells.data()+disp[id],send[id],id,MPI_COMM_WORLD, &requests_cast_phys[id]);
    }

    //здесь можно сделать  так:

    int rcvCast=np-1; //кроме себя
    do
    {
        flag=false;
        for(int id=0; id<np;id++)
        {
            if(myid != id && requests_cast_phys[id]!= NULL)
            if(Mpi_test(requests_cast_phys[id]))
            {
                //calc(alpha,betta,T,logT)[disp[id]:send[id]+disp[id]]
                --rcvCast;
                 requests_cast_phys[id]!= NULL;
            }
        }
    }while(rcvCast);


    Mpiwaitall(requests_cast_phys.size(),requests_cast_phys.data());
}
#endif