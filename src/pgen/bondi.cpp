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
  void AddBHGrav(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                      const DvceArray5D<Real> &w0, const EOS_Data &eos_data);
                                       
  struct pgen_bh{
      // jet & ambient medium parameters
      Real cs_amb, d_amb, M_bh, gamma, CONST_G, epsilon;
  };
      pgen_bh* pbh = new pgen_bh();
  // Function to check if a point is inside the jet
//   bool InJet(const Real x1v, const Real x2v, const Real x3v, const Real r_jet, const Real alpha_jet, const Real theta_jet) {
//     Real r = sqrt(SQR(x1v) + SQR(x2v) + SQR(x3v));
//     Real nhat_dot_rhat = (r > 0) ? (x1v*cos(theta_jet) + x2v*sin(theta_jet))/r : 0.0;
//     return (r < r_jet) && (abs(nhat_dot_rhat) > cos(alpha_jet));
//   }
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
    pbh->cs_amb        = pin->GetReal("problem", "cs_amb");
    pbh->d_amb         = pin->GetReal("problem", "d_amb");
    pbh->M_bh          = pin->GetReal("problem", "M_bh");
    pbh->epsilon       = pin->GetReal("problem","epsilon");
    pbh->gamma         = pin->GetReal("hydro", "gamma");
    pbh->CONST_G       = pin->GetReal("problem", "CONST_G");

    Real const &gm1     = pbh->gamma - 1;

    if (restart) return;
    // Select either Hydro or MHD
    if (pmbp->phydro!=nullptr){

        auto &u0 = pmbp->phydro->u0;
        par_for("bondi",DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
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
            Real pres = pbh->d_amb*SQR(pbh->cs_amb); //jet has same pressure as ambient medium
            // std::cout << "pres= " << pres << "\n";
            u0(m,IDN,k,j,i) = pbh->d_amb;
            u0(m,IM1,k,j,i) = 0.0;
            u0(m,IM2,k,j,i) = 0.0;
            u0(m,IM3,k,j,i) = 0.0;
            u0(m,IEN,k,j,i) = pres/gm1 + 0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))+SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
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
        AddBHGrav(pm,bdt,u0,w0,eos_data);
        return;
    }

    void AddBHGrav(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                const DvceArray5D<Real> &w0, const EOS_Data &eos_data) { //Apply BH Grav at all timesteps
        MeshBlockPack *pmbp = pm->pmb_pack;
        auto &indcs = pmbp->pmesh->mb_indcs;
        int is = indcs.is, ie = indcs.ie;
        int js = indcs.js, je = indcs.je;
        int ks = indcs.ks, ke = indcs.ke;
        int nmb1 = pmbp->nmb_thispack - 1;
        auto size = pmbp->pmb->mb_size;

        Real const &gm1 = pbh->gamma - 1;

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

            Real rad = sqrt(SQR(x1v)+SQR(x2v)+SQR(x3v));

            Real grad_phi_by_r = (pbh->CONST_G*pbh->M_bh)/(pow((SQR(rad)+SQR(pbh->epsilon)),1.5));

            u0(m,IM1,k,j,i) -= (pbh->d_amb*grad_phi_by_r*bdt)*x1v;
            u0(m,IM2,k,j,i) -= (pbh->d_amb*grad_phi_by_r*bdt)*x2v;
            u0(m,IM3,k,j,i) -= (pbh->d_amb*grad_phi_by_r*bdt)*x3v;
            Real p_dot_r = (u0(m,IM1,k,j,i)*x1v)+(u0(m,IM2,k,j,i)*x2v)+(u0(m,IM3,k,j,i)*x3v);
            u0(m,IEN,k,j,i) -=p_dot_r*grad_phi_by_r*bdt;

        });
        return;
    }
}; //namespace
