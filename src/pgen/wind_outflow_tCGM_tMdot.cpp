//========================================================================================
// Athena++ astrophysical MHD code, Kokkos version
// Copyright(C) 2020 James M. Stone <jmstone@ias.edu> and the Athena code team
// Licensed under the 3-clause BSD License (the "LICENSE")
//========================================================================================
//! \file wind_outflow_time_dependent_cgm.cpp
//! \brief Problem generator for time-dependent CGM winds blowing in a stationary BH
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
#include "srcterms/ismcooling.hpp"
#include "globals.hpp"
#include "units/units.hpp"

//----------------------------------------------------------------------------------------
//! \fn ProblemGenerator::UserProblem_()
//! \brief Problem Generator for jets in a uniform medium
namespace {
    // made global to share with source terms   
    // made global to share with source terms   
    struct PgenBh{
        Real CONST_G,
        CONST_kB_cgs, CONST_mp, CONST_mu, 
        M_bh, v_bh, epsilon,
        M_dot_0,tau_w, v_w, r_inj,
        rho0_cgm, cs_cgm, tau, gamma_gas, 
        length_cgs, mass_cgs, time_cgs, temp_floor;
        KOKKOS_INLINE_FUNCTION PgenBh() = default;
    };
    PgenBh pbh;  // Host global
    void AddUserSrcs(Mesh *pm, const Real bdt);
    void AddISMCooling(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                      const DvceArray5D<Real> &w0, const EOS_Data &eos_data);
    void AddBHGrav(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                      const DvceArray5D<Real> &w0, const EOS_Data &eos_data);
    void AddOutflow(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                      const DvceArray5D<Real> &w0, const EOS_Data &eos_data);
    void WindX1BC(Mesh *pm);            
    // Functions Time dependent CGM
    KOKKOS_INLINE_FUNCTION
    static Real Rho_CGM(const Real t,const Real tau, const Real rho0_CGM) {
        Real rho_cgm = rho0_CGM*exp((-1*t)/tau);
        return rho_cgm;
  }
  KOKKOS_INLINE_FUNCTION
    static Real M_dot_w(const Real t,const Real tau_w, const Real M_dot_w) {
        Real Mdot_w = M_dot_w*exp((-1*t)/tau_w);
        return Mdot_w;
  }
} // namespace
  

void ProblemGenerator::UserProblem(ParameterInput *pin, const bool restart) {
    user_srcs_func = AddUserSrcs;
     // User boundary function
    user_bcs_func = WindX1BC;
    MeshBlockPack *pmbp = pmy_mesh_->pmb_pack;

    // capture variables for the kernel
    auto &indcs = pmbp->pmesh->mb_indcs;
    int &is = indcs.is; int &ie = indcs.ie;
    int &js = indcs.js; int &je = indcs.je;
    int &ks = indcs.ks; int &ke = indcs.ke;
    auto &size = pmbp->pmb->mb_size;

    int nfluid   = pmbp->phydro->nhydro;   
    int nscalars = pmbp->phydro->nscalars;

    // get initial parameters from input file
    pbh.CONST_G        = pin->GetReal("problem","CONST_G");    
    pbh.CONST_kB_cgs   = pin->GetReal("problem","CONST_kB_cgs"); 
    pbh.CONST_mp       = pin->GetReal("problem","CONST_mp");    
    pbh.CONST_mu       = pin->GetReal("problem","CONST_mu"); 
    pbh.gamma_gas      = pin->GetReal("hydro", "gamma");

    pbh.rho0_cgm       = pin->GetReal("problem","rho0_cgm");    
    pbh.cs_cgm         = pin->GetReal("problem","cs_cgm");
    pbh.tau            = pin->GetReal("problem", "tau");
    // pbh.alpha_cgm      = pin->GetReal("problem", "alpha_cgm");

    pbh.M_bh           = pin->GetReal("problem","M_bh");
    pbh.v_bh           = pin->GetReal("problem","v_bh");
    pbh.epsilon        = pin->GetReal("problem","epsilon");

    pbh.M_dot_0        = pin->GetReal("problem","M_dot_0");
    pbh.tau_w          = pin->GetReal("problem","tau_w");
    pbh.v_w            = pin->GetReal("problem","v_w");
    pbh.r_inj          = pin->GetReal("problem","r_inj");

    pbh.temp_floor     = pin->GetOrAddReal("problem", "temp_floor", 1.0e4 );
    pbh.length_cgs     = pin->GetReal("units", "length_cgs");
    pbh.mass_cgs       = pin->GetReal("units", "mass_cgs");
    pbh.time_cgs       = pin->GetReal("units", "time_cgs");

    Real const &gm1     = pbh.gamma_gas - 1;

    Real time           = pmbp->pmesh->time;

    if (restart) return;

    // Select either Hydro or MHD
    if (pmbp->phydro!=nullptr){

        auto &u0 = pmbp->phydro->u0;
        par_for("rad_bondi",DevExeSpace(),0,(pmbp->nmb_thispack-1),ks,ke,js,je,is,ie,
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

            Real rad = sqrt(SQR(x1v)+SQR(x2v)+SQR(x3v));

            if (rad<pbh.r_inj){
                Real v_inj      = pbh.v_w/std::sqrt(2);
                Real Mdot_w     = M_dot_w(time,pbh.tau_w,pbh.M_dot_0);
                Real rho_w      = (Mdot_w)/(4*M_PI*SQR(pbh.r_inj)*v_inj);
                Real pres_w     = (rho_w*SQR(v_inj))/pbh.gamma_gas;
                u0(m,IDN,k,j,i) = rho_w;
                u0(m,IM1,k,j,i) = (rho_w*v_inj*x1v)/rad;
                u0(m,IM2,k,j,i) = (rho_w*v_inj*x2v)/rad;
                u0(m,IM3,k,j,i) = (rho_w*v_inj*x3v)/rad;
                u0(m,IEN,k,j,i) = pres_w/gm1 + 0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))+SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
                //Add outflow tracer
                u0(m, nfluid, k, j, i) = rho_w*1.0;
            }
            else{
                Real rho_cgm = Rho_CGM(time, pbh.tau, pbh.rho0_cgm);
                Real pres = rho_cgm*SQR(pbh.cs_cgm);
                u0(m,IDN,k,j,i) = rho_cgm;
                u0(m,IM1,k,j,i) = -1*rho_cgm*pbh.v_bh;
                u0(m,IM2,k,j,i) = 0.0;
                u0(m,IM3,k,j,i) = 0.0;
                u0(m,IEN,k,j,i) = pres/gm1 + 0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))+SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
                //Add tracer = 0 outside
                u0(m, nfluid, k, j, i) = 0.0;
            }

            // //Initialize pressure and density profile
            // Real pres = pbh->d_amb*SQR(pbh->cs_amb); //jet has same pressure as ambient medium
            // // std::cout << "pres= " << pres << "\n";
            // u0(m,IDN,k,j,i) = pbh->d_amb;
            // u0(m,IM1,k,j,i) = 0.0;
            // u0(m,IM2,k,j,i) = 0.0;
            // u0(m,IM3,k,j,i) = 0.0;
            // u0(m,IEN,k,j,i) = pres/gm1 + 0.5*(SQR(u0(m,IM1,k,j,i))+SQR(u0(m,IM2,k,j,i))+SQR(u0(m,IM3,k,j,i)))/u0(m,IDN,k,j,i);
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
        AddOutflow(pm,bdt,u0,w0,eos_data);
        AddISMCooling(pm,bdt,u0,w0,eos_data);
        AddBHGrav(pm,bdt,u0,w0,eos_data);
        return;
    }
    void AddISMCooling(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                const DvceArray5D<Real> &w0, const EOS_Data &eos_data) { //Apply BH Grav at all timesteps
        MeshBlockPack *pmbp = pm->pmb_pack;
        auto &indcs = pmbp->pmesh->mb_indcs;
        int is = indcs.is, ie = indcs.ie;
        int js = indcs.js, je = indcs.je;
        int ks = indcs.ks, ke = indcs.ke;
        int nmb1 = pmbp->nmb_thispack - 1;
        auto size = pmbp->pmb->mb_size;

        Real const &gm1         = pbh.gamma_gas - 1;
        Real rho_unit           = pbh.mass_cgs/(std::pow(pbh.length_cgs,3.0));
        Real v_unit             = pbh.length_cgs/pbh.time_cgs;
        Real temp_unit          = pbh.CONST_mu*pbh.CONST_mp/(pbh.CONST_kB_cgs);
        Real energy_unit        = pbh.mass_cgs * SQR(v_unit);
        Real vol_unit           = pow(pbh.length_cgs, 3.0);
        Real cooling_rate_unit  = energy_unit / (vol_unit * pbh.time_cgs);
        Real pres_unit          = rho_unit*v_unit*v_unit;
        Real time               = pmbp->pmesh->time;

        // Real rho_cgm            = Rho_CGM(time, pbh.tau, pbh.rho0_cgm);
        Real T_cgm              = 1.0e6;
        // Real T_thresh_low       = 0.9*T_cgm;
        // Real T_thresh_high      = 1.1*T_cgm;
        Real T_thresh_low       = 0.5*T_cgm;
        Real T_thresh_high      = 1.5*T_cgm;        
        // std::cout<<"rho_unit = "<<rho_unit<<"v_unit = "<<v_unit<<"temp_unit = "<<temp_unit<<"pres_unit = "<<pres_unit<<" cool_rate_unit = "<<cooling_rate_unit<<std::endl;

        par_for("ism_cooling", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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

            //Compute Pressure Unit
            
            //Calculate temperature in CGS
            Real rho_code   = w0(m, IDN, k, j, i);
            
            Real dens_cgs   = rho_code*rho_unit;
            Real pres_cgs   = w0(m, IPR, k, j, i)*pres_unit;
            Real temp_cgs   = (pres_cgs/dens_cgs)*temp_unit;
            Real n_cgs      = dens_cgs/(pbh.CONST_mu*pbh.CONST_mp);
            //Put a hard temperature floor
            if (temp_cgs<pbh.temp_floor){
                Real pres_floor_cgs     = (dens_cgs*pbh.CONST_kB_cgs*pbh.temp_floor)/(pbh.CONST_mu*pbh.CONST_mp);
                Real pres_floor_code    = pres_floor_cgs/pres_unit;
                w0(m,IPR,k,j,i)         = pres_floor_code;
                // std::cout<<"Temp less than 1e4K, temp_cgs = "<<temp_cgs<<" pres_cgs = "<<pres_cgs<<" dens_cgs = "<<dens_cgs<<" pres_floor = "<<pres_floor_cgs<<std::endl;
                Real vx     = w0(m,IVX,k,j,i);
                Real vy     = w0(m,IVY,k,j,i);
                Real vz     = w0(m,IVZ,k,j,i);
                Real rho    = w0(m,IDN,k,j,i);

                Real kin    = 0.5 * rho * (SQR(vx) + SQR(vy) + SQR(vz));
                u0(m,IEN,k,j,i) = pres_floor_code/gm1 + kin;
            }else if (temp_cgs<T_thresh_low or temp_cgs>T_thresh_high){
                //Add Radiative cooling if temp>temp_floor
                Real lambda_cgs         = ISMCoolFn(temp_cgs);
                Real cooling_rate_cgs   = SQR(n_cgs)*lambda_cgs;
                // std::cout<<"T>1e4, temp_cgs = "<<temp_cgs<<" lambda_cgs = "<<lambda_cgs<<" cool_rate_cgs = "<<cooling_rate_cgs<<std::endl;
                //Convert to code units
                Real cooling_rate_code  = cooling_rate_cgs/cooling_rate_unit;
                //Update internal Energy
                u0(m, IEN, k, j, i)     -= cooling_rate_code*bdt;
            }
            

        });
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

        Real const &gm1 = pbh.gamma_gas - 1;

        par_for("bh_gravity", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
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

            Real grad_phi_by_r = (pbh.CONST_G*pbh.M_bh)/(pow((SQR(rad)+SQR(pbh.epsilon)),1.5));

            u0(m,IM1,k,j,i) -= (u0(m,IDN,k,j,i)*grad_phi_by_r*bdt)*x1v;
            u0(m,IM2,k,j,i) -= (u0(m,IDN,k,j,i)*grad_phi_by_r*bdt)*x2v;
            u0(m,IM3,k,j,i) -= (u0(m,IDN,k,j,i)*grad_phi_by_r*bdt)*x3v;
            Real p_dot_r = (u0(m,IM1,k,j,i)*x1v)+(u0(m,IM2,k,j,i)*x2v)+(u0(m,IM3,k,j,i)*x3v);
            u0(m,IEN,k,j,i) -=p_dot_r*grad_phi_by_r*bdt;

        });
        return;
    }
    void AddOutflow(Mesh *pm, const Real bdt, DvceArray5D<Real> &u0,
                const DvceArray5D<Real> &w0, const EOS_Data &eos_data) { //Apply Radial Outflow at all timesteps
        MeshBlockPack *pmbp = pm->pmb_pack;
        auto &indcs = pmbp->pmesh->mb_indcs;
        int is = indcs.is, ie = indcs.ie;
        int js = indcs.js, je = indcs.je;
        int ks = indcs.ks, ke = indcs.ke;
        int nmb1 = pmbp->nmb_thispack - 1;
        auto size = pmbp->pmb->mb_size;

        // number of hydro vars and passive scalars
        int nfluid   = pmbp->phydro->nhydro;
        int nscalars = pmbp->phydro->nscalars;
        Real time       = pm->pmb_pack->pmesh->time;

        Real const &gm1 = pbh.gamma_gas - 1; 

        par_for("outflow", DevExeSpace(), 0, nmb1, ks, ke, js, je, is, ie,
        KOKKOS_LAMBDA(const int m, const int k, const int j, const int i) {

            Real &x1min = size.d_view(m).x1min;
            Real &x1max = size.d_view(m).x1max;
            int nx1     = indcs.nx1;
            Real x1v    = CellCenterX(i-is, nx1, x1min, x1max);

            Real &x2min = size.d_view(m).x2min;
            Real &x2max = size.d_view(m).x2max;
            int nx2     = indcs.nx2;
            Real x2v    = CellCenterX(j-js, nx2, x2min, x2max);

            Real &x3min = size.d_view(m).x3min;
            Real &x3max = size.d_view(m).x3max;
            int nx3     = indcs.nx3;
            Real x3v    = CellCenterX(k-ks, nx3, x3min, x3max);

            Real rad    = sqrt(SQR(x1v)+SQR(x2v)+SQR(x3v));

            Real v_inj  = pbh.v_w/std::sqrt(2);
            Real Vol_inj= (4*M_PI*std::pow(pbh.r_inj,3))/3;
            Real Mdot_w = M_dot_w(time, pbh.tau_w, pbh.M_dot_0);    
            Real S_rho  = Mdot_w/Vol_inj;
            Real S_e    = (S_rho*SQR(v_inj))*(0.5+(1/(pbh.gamma_gas*gm1))); 

            if (rad<pbh.r_inj){
                Real dens       = w0(m,IDN,k,j,i);     
                u0(m,IDN,k,j,i) += S_rho*bdt;
                u0(m,IM1,k,j,i) += (S_rho*v_inj*x1v/rad)*bdt;
                u0(m,IM2,k,j,i) += (S_rho*v_inj*x2v/rad)*bdt;
                u0(m,IM3,k,j,i) += (S_rho*v_inj*x3v/rad)*bdt;
                u0(m,IEN,k,j,i) += S_e*bdt;
                //Add tracer
                u0(m,nfluid,k,j,i) = dens*1.0;
            }
        });
        return;
    }
    void WindX1BC(Mesh *pm) {
        auto &indcs     = pm->mb_indcs;
        int &ng         = indcs.ng;
        int n1          = indcs.nx1 + 2*ng;
        int n2          = (indcs.nx2 > 1)? (indcs.nx2 + 2*ng) : 1;
        int n3          = (indcs.nx3 > 1)? (indcs.nx3 + 2*ng) : 1;
        int &is         = indcs.is;  int &ie  = indcs.ie;
        int &js         = indcs.js;  int &je  = indcs.je;
        int &ks         = indcs.ks;  int &ke  = indcs.ke;
        auto &mb_bcs    = pm->pmb_pack->pmb->mb_bcs;

        DvceArray5D<Real> u0_, w0_;
        u0_             = pm->pmb_pack->phydro->u0;
        w0_             = pm->pmb_pack->phydro->w0;
        int nmb         = pm->pmb_pack->nmb_thispack;
        int nvar        = u0_.extent_int(1);

        Real time       = pm->pmb_pack->pmesh->time;
        int nfluid_bc   = pm->pmb_pack->phydro->nhydro;
        int nscalars_bc = pm->pmb_pack->phydro->nscalars;
        Real rho_cgm    = Rho_CGM(time, pbh.tau, pbh.rho0_cgm);
        Real pres_cgm   = rho_cgm*SQR(pbh.cs_cgm);
        // std::cout << "time= " << time <<"\t rho = "<<rho_cgm_t<<"\t pres = "<<pres_cgm_t<<"\n"; //Code Units
        // Real v_bh     = pbh->v_bh;  // Set wind speed here

        // // Copy to device-visible scalars for the lambda
        // Real rho_cgm_t_loc  = rho_cgm_t;
        // Real pres_cgm_t_loc = pres_cgm_t;
        // Copy to device-visible scalars for the lambda
        // Real rho_cgm_t_loc  = pbh->rho0_cgm;
        // Real pres_cgm_t_loc = pbh->rho0_cgm*SQR(pbh->cs_cgm);
        // Real v_bh_loc       = pbh->v_bh;


        // ConsToPrim over all X1 ghost zones
        pm->pmb_pack->phydro->peos->ConsToPrim(u0_, w0_, false, is-ng, ie+ng, 0, n2-1, 0, n3-1);

        // Set outer X1 ghost zones to fixed CGM inflow
        par_for("wind_x1_bc", DevExeSpace(), 0, nmb-1, ks, ke, js, je,
        KOKKOS_LAMBDA(int m, int k, int j) {
            if (mb_bcs.d_view(m, BoundaryFace::outer_x1) == BoundaryFlag::user) {
                for (int i = 0; i < ng; ++i) {
                    int ig = ie + i + 1;

                    // Density
                    w0_(m, IDN, k, j, ig) = rho_cgm;

                    // Velocity: inflow with -v_bh along x
                    w0_(m, IVX, k, j, ig) = -1*pbh.v_bh;
                    w0_(m, IVY, k, j, ig) = 0.0;
                    w0_(m, IVZ, k, j, ig) = 0.0;

                    // Pressure = rho_cgm * cs_cgm^2
                    w0_(m, IPR, k, j, ig) = rho_cgm*SQR(pbh.cs_cgm);
                    //Add tracer = 0
                    if (nscalars_bc > 0) {
                        w0_(m, nfluid_bc, k, j, ig) = 0.0;  // tracer=0 for CGM inflow
                    }
                }
            }
        });

        // PrimToCons on X1 ghost zones
        pm->pmb_pack->phydro->peos->PrimToCons(w0_,u0_,ie+1,ie+ng,0,(n2-1),0,(n3-1));

        return;
    }

}; //namespace
