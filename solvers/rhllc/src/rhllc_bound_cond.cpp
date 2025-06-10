#ifdef RHLLC
#include "rhllc_bound_cond.h"

#include "rhllc_flux.h"

#include "global_consts.h"
#include "gas_state.h"
#include "linear_alg.h"

void rhllc::bound::OutSrcBound(const elem_t &cell, flux_all_t &flx)
{
#if 0    
    constexpr Type rho = kM_hydrogen * 1e14 / kDensity; //  0.1;
    constexpr Type prs = GetPressure(rho, 10 * kEv);
    constexpr Type v[3] = {1e3 / kVelocity, 0, 0};

    static const flux_t phys(rho, Vector3(v[0], v[1], v[2]), prs);
    static const flux_t conv = GetConvValue(phys);

    flx.phys_val = phys;
    flx.phys_val = conv;
#endif

    flx.conv_val = cell.conv_val;
    flx.phys_val = cell.phys_val;
}
void rhllc::bound::InSrcBound(const elem_t &cell, flux_all_t &flx)
{
    flx.conv_val = cell.conv_val;
    flx.phys_val = cell.phys_val;
}
void rhllc::bound::OuterSurfBound(const elem_t &cell, flux_all_t &flx)
{
    flx.conv_val = cell.conv_val;
    flx.phys_val = cell.phys_val;
}

void rhllc::BoundConditions(const face_t &f, const std::vector<elem_t> &cells, flux_all_t &bound)
{
    const int id_l = f.geo.id_l;
    const int id_r = f.geo.id_r;
    const elem_t &cell = cells[id_l];
    Matrix3 T;

    switch (id_r) // id соседа она же признак ГУ
    {

    case e_bound_free:
    {
        bound.conv_val = cell.conv_val;
        bound.phys_val = cell.phys_val;
    }
    break;
    case e_bound_lock:
    {
        GetRotationMatrix(f.geo.n, T);

        bound.conv_val.v = T * bound.conv_val.v;
        bound.phys_val.v = T * bound.phys_val.v;

        bound.conv_val.v[0] = -bound.conv_val.v[0];
        bound.phys_val.v[0] = -bound.phys_val.v[0];

        bound.conv_val.v = T.transpose() * bound.conv_val.v;
        bound.phys_val.v = T.transpose() * bound.phys_val.v;
    }
    break;

    case e_bound_out_source:
    {
        bound::OutSrcBound(cell, bound);
    }
    break;

    case e_bound_inner_source:
    {
        bound::InSrcBound(cell, bound);
    }
    break;

    case e_bound_outer_surface:
    {
        bound::OuterSurfBound(cell, bound);
    }
    break;

    default:
        STOP_IF(id_r < 0); // Err bound in RHLLC_3d

        bound.conv_val = cells[id_r].conv_val;
        bound.phys_val = cells[id_r].phys_val;
        break;
    }
}
#endif