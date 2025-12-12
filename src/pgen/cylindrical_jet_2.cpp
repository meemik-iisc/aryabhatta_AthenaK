//========================================================================================
// Athena++ astrophysical MHD code, Kokkos version
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file jet_amb.cpp
//! \brief Problem generator for an AGN Jet in ambient medium
//! \ref (arXiv:2401.00446v1) Dissipation of AGN Jets in clumpy interstellar medium

#include <algorithm>
#include <cmath>
#include <sstream>

#include "parameter_input.hpp"
#include "athena.hpp"
#include "mesh/mesh.hpp"
#include "eos/eos.hpp"
#include "hydro/hydro.hpp"
#include "coordinates/cartesian_ks.hpp"
#include "coordinates/cell_locations.hpp"
#include "pgen.hpp"
#include "srcterms/srcterms.hpp"
#include "globals.hpp"
#include "units/units.hpp"

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem_()
//! \brief Problem Generator for jets in a uniform medium
namespace {
    // made global to share with source terms   
    void AddUserSrcs(Mesh *pm, const Real bdt);
    void AddJets(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                        const DvceArray5D<Real> &w0, const EOS_Data &eos_data);                 
    struct pgen_jet_amb{
        // jet & ambient medium parameters
        Real cs_amb, d_amb, rho_jet, gamma, r_jet, h_jet, l_jet, v_jet;
    };
        pgen_jet_amb* pjet = new pgen_jet_amb();
    // Function to check if a point is inside the jet
    bool InJet(const Real x1v, const Real x2v, const Real x3v, const Real r_jet, const Real h_jet) {
        Real r2 = x2v*x2v+x3v*x3v;
        return (r2<=r_jet*r_jet) && (abs(x1v)<=h_jet);
    }
} // namespace

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
    user_srcs_func = AddUserSrcs;

    MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;
    // capture variables for the kernel
    auto &indcs = pmbp->pmesh->mb_indcs;
    int &is = indcs.is; int &ie = indcs.ie;
    int &js = indcs.js; int &je = indcs.je;
    int &ks = indcs.ks; int &ke = indcs.ke;
    auto &size = pmbp->pmb->mb_size;
    // get initial parameters from input file
    pjet->gamma       = pin->GetReal("hydro", "gamma");
    pjet->cs_amb      = pin->GetReal("problem", "cs_amb");
    pjet->d_amb       = pin->GetOrAddReal("problem", "d_amb",1);
    pjet->l_jet       = pin->GetReal("problem", "l_jet");
    pjet->r_jet       = pin->GetReal("problem", "r_jet");
    pjet->h_jet       = pin->GetReal("problem", "h_jet");
    pjet->v_jet       = pin->GetReal("problem", "v_jet");
    

    Real area_jet       = M_PI*pjet->r_jet*pjet->r_jet;
    pjet->rho_jet      = (2*pjet->l_jet)/(pjet->v_jet*pjet->v_jet*pjet->v_jet*area_jet);
    Real const &gm1    = pjet->gamma - 1;

    if (restart) return;

    // Select either Hydro or MHD
    if (pmbp->phydro!=nullptr){
        auto &u0 = pmbp->phydro->u0;
        int nfluid = pmbp->phydro->nhydro;
        int nscalars = pmbp->phydro->nscalars;
        par_for("jets",DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
        KOKKOS_LAMBDA(int m,int k,int j,int i) {
            Real &xmin = size.d_view(m).x1min;
            Real &xmax = size.d_view(m).x1max;
            int nx1 = indcs.nx1;
            Real x1v = CellCenterX(i-is, nx1, xmin, xmax);

            Real &ymin = size.d_view(m).x2min;
            Real &ymax = size.d_view(m).x2max;
            int nx2 = indcs.nx2;
            Real x2v = CellCenterX(j-js, nx2, ymin, ymax);

            Real &zmin = size.d_view(m).x3min;
            Real &zmax = size.d_view(m).x3max;
            int nx3 = indcs.nx3;
            Real x3v = CellCenterX(k-ks, nx3, zmin, zmax);

            //Initialize pressure and density profile
            Real pres = pjet->d_amb*SQR(pjet->cs_amb); //jet has same pressure as ambient medium
            Real scal = 0.0;

            if (InJet(x1v, x2v, x3v, pjet->r_jet, pjet->h_jet)) {
                scal =1.0;
                u0(m,IDN,k,j,i) = pjet->rho_jet;
                Real        rad = sqrt(SQR(x1v) + SQR(x2v) + SQR(x3v));
                u0(m,IM1,k,j,i) = (rad > 0) ? pjet->rho_jet*pjet->v_jet : 0.0;
                u0(m,IM2,k,j,i) = 0.0;
                u0(m,IM3,k,j,i) = 0.0;
            } else {
                scal = 0.0;
                u0(m,IDN,k,j,i) = pjet->d_amb;
                u0(m,IM1,k,j,i) = 0.0;
                u0(m,IM2,k,j,i) = 0.0;
                u0(m,IM3,k,j,i) = 0.0;
            }
            u0(m,IEN,k,j,i) = pres/gm1 + 0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))+SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
            // add passive scalars
            for (int n=nfluid; n<(nfluid+nscalars); ++n) {
                u0(m,n,k,j,i) = scal;
            }
        });
    return;
    }
}

namespace{
    //----------------------------------------------------------------------------------------
    //! \fn void AddUserSrcs()
    //! \brief Add User Source Terms
    // NOTE source terms must all be computed using primitive (w0) and NOT conserved (u0) vars
    void AddUserSrcs(Mesh *pm, const Real bdt) {
        MeshBlockPack *pmbp = pm->pmb_pack;
        const auto &w0 = pmbp->phydro->w0;
        auto &u0 = pmbp->phydro->u0;
        const EOS_Data &eos_data = pmbp->phydro->peos->eos_data;
        AddJets(pm,bdt,u0,w0,eos_data);
        return;
    }

    void AddJets(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                const DvceArray5D<Real> &w0, const EOS_Data &eos_data) { //Apply jet injection at all time substeps
        MeshBlockPack *pmbp = pm->pmb_pack;
        auto &indcs = pmbp->pmesh->mb_indcs;
        int is = indcs.is, ie = indcs.ie;
        int js = indcs.js, je = indcs.je;
        int ks = indcs.ks, ke = indcs.ke;
        int nmb1 = pmbp->nmb_thispack - 1;
        auto size = pmbp->pmb->mb_size;
        Real const &gm1 = pjet->gamma - 1;
        int nfluid = pmbp->phydro->nhydro;
        int nscalars = pmbp->phydro->nscalars;

        par_for("jet_inject", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

            Real &x1min = size.d_view(m).x1min;
            Real &x1max = size.d_view(m).x1max;
            int nx1 = indcs.nx1;
            Real x1v = CellCenterX(i-is, nx1, x1min, x1max);

            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            int nx2 = indcs.nx2;
            Real x2v = CellCenterX(j-js, nx2, x2min, x2max);

            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            int nx3 = indcs.nx3;
            Real x3v = CellCenterX(k-ks, nx3, x3min, x3max);

            if (InJet(x1v, x2v, x3v, pjet->r_jet, pjet->h_jet)) {
                Real scal = 1.0;
                Real pres = pjet->d_amb*SQR(pjet->cs_amb); //jet has same pressure as ambient medium
                u0(m,IDN,k,j,i) = pjet->rho_jet;
                Real rad = sqrt(SQR(x1v) + SQR(x2v) + SQR(x3v));
                u0(m,IM1,k,j,i) = (rad > 0) ? pjet->rho_jet*pjet->v_jet : 0.0;
                u0(m,IM2,k,j,i) = 0.0;
                u0(m,IM3,k,j,i) = 0.0;
                u0(m,IEN,k,j,i) = pres/gm1 + 0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))+SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
                // add passive scalars
                for (int n=nfluid; n<(nfluid+nscalars); ++n) {
                    u0(m,n,k,j,i) = scal;
                }
            }
        });
        return;
    }

}; //namespace
