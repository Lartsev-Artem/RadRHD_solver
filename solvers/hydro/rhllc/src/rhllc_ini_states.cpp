#ifdef RHLLC
#include "rhllc_ini_states.h"

#include "rhllc_utils.h"

#include "global_consts.h"
#include "gas_state.h"

#include "global_value.h"
#include "reader_bin.h"

using namespace rrhd;

flux_t rhllc::ini::Soda(const Vector3 &x)
{
    flux_t flx;
#if 1
    if (x[0] < 0.5)
    {
        flx.d = 1;
        flx.p = 1;
        flx.v = Vector3(0.9, 0, 0);
    }
    else
    {
        flx.d = 1;
        flx.p = 10;
        flx.v = Vector3(0, 0, 0);
    }
#else
    if (x[0] < 0.5)
    {
        flx.d = 1;
        flx.p = 1;
        flx.v = Vector3(0, 0, 0);
    }
    else
    {
        flx.d = 0.125;
        flx.p = 0.1;
        flx.v = Vector3(0, 0, 0);
    }
#endif
    return flx;
}
flux_t rhllc::ini::Jet(const Vector3 &x)
{
    flux_t flx;
    if (Vector2(x[1], x[2]).norm() < 0.03 && x[0] < 0.1)
    {
        flx.d = 0.1;
        flx.p = 0.01;
        flx.v = Vector3(0.99, 0, 0);
    }
    else
    {
        flx.d = 10;
        flx.p = 0.01;
        flx.v = Vector3(0, 0, 0);
    }
    return flx;
}

flux_t rhllc::ini::Uniform(const Vector3 &x)
{
    flux_t flx;
#if 1
    flx.d = kM_hydrogen * 1e14 / kDensity; //  0.1;
    flx.p = GetPressure(flx.d, kEv);
    flx.v = Vector3(0, 0, 0);
#else
    flx.d = 0.1;
    flx.p = 0.1;
    flx.v = Vector3(0, 0, 0);
#endif
    return flx;
}

int rhllc::Init(const std::string &file_init_value,
                const std::function<flux_t(const Vector3 &)> ini_func,
                grid_t &grid)
{

    if (files_sys::bin::ReadHllcInit(file_init_value, grid.cells) != e_completion_success)
    {
        std::vector<Vector3> centers;
        if (files_sys::bin::ReadSimple(glb_files.base_address + F_CENTERS, centers))
            return e_completion_fail;

#pragma omp parallel for
        for (int i = 0; i < grid.size; i++)
        {
            grid.cells[i].phys_val = ini_func(centers[i]);
        }
        WRITE_LOG("SetDefault hllc value\n");
    }

    HllcPhysToConv(grid);
    rhllc::max_signal_speed = 1;
    _hllc_cfg.tau = GetTimeStep(_hllc_cfg);
    return e_completion_success;
}
#endif