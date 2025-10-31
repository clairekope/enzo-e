
/// See LICENSE_CELLO file for license and copyright information

/// @file	  enzo_EnzoMethodFeedbackMechanical.cpp
/// @author Claire Kopenhafer (kopenhaf@msu.edu, clairekope@gmail.com)
//             
/// @date
/// @brief  Implements the mechanical feedback method in Enzo that was
///         written by Cassi Lochhaas for the FOGGIE collaboration.
///         This method was based on Kimm & Cen 2014. Since this method
///         is similar to the STARSS algorithm, this Enzo-E adaptation
///         is built using the STARSS implementation as a starting point.


#include "Cello/cello.hpp"
#include "Enzo/enzo.hpp"
#include "Enzo/particle/particle.hpp"

#include <time.h>
#include <random>


// =============================================================================
/* 
  Following Azton & Andrew's comments in STARSS & Distributed FB methods,
  I should create a helper class that performs table lookups
  and potentially other methods as needed. 
  An example is EnzoMethodM1Closure and the M1Tables class.
*/
// =============================================================================

void EnzoMethodFeedbackMechanical::compute_molecular_weight_(EnzoBlock * enzo_block) {

  Field field = enzo_block->data()->field();

  int mx, my, mz, gx, gy, gz, nx, ny, nz;
  field.size(&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;

  CelloView<enzo_float,3> d  = field.view<enzo_float>("density");

  // this temp field must be allocated before this function is called
  CelloView<enzo_float,3> mu = field.view<enzo_float>(i_mu);

  CelloView<enzo_float,3> dHI, dHII, dHeI, dHeII, dHeIII, d_el,
    dH2I, dH2II, dHM, dDI, dDII, dHDI;

  if (field.is_field("HI_density"))  dHI   = field.view<enzo_float>("HI_density");
  if (field.is_field("HII_density")) dHII  = field.view<enzo_float>("HII_density");
  if (field.is_field("HeI_density")) dHeI  = field.view<enzo_float>("HeI_density");
  if (field.is_field("HeII_density")) dHeII = field.view<enzo_float>("HeII_density");
  if (field.is_field("HeIII_density")) dHeIII = field.view<enzo_float>("HeIII_density");
  if (field.is_field("e_density"))    d_el  = field.view<enzo_float>("e_density");
  if (field.is_field("H2I_density"))  dH2I  = field.view<enzo_float>("H2I_density");
  if (field.is_field("H2II_density")) dH2II = field.view<enzo_float>("H2II_density");
  if (field.is_field("HM_density"))  dHM   = field.view<enzo_float>("HM_density");
  if (field.is_field("DI_density"))  dDI   = field.view<enzo_float>("DI_density");
  if (field.is_field("DII_density")) dDII  = field.view<enzo_float>("DII_density");
  if (field.is_field("HDI_density")) dHDI  = field.view<enzo_float>("HDI_density");

  const GrackleChemistryData * grackle_chem = enzo::grackle_chemistry();
  const int primordial_chemistry = (grackle_chem == nullptr) ?
    0 : grackle_chem->get<int>("primordial_chemistry");

  const double dflt_mu = static_cast<double>(enzo::fluid_props()->mol_weight());

  for (int iz=0; iz<mz; iz++) {
    for (int iy=0; iy<my; iy++) {
      for (int ix=0; ix<my; ix++) {
        if (primordial_chemistry > 0) {
          mu(ix, iy, iz) = d_el(ix, iy, iz) + dHI(ix, iy, iz) + dHII(ix, iy, iz) 
             + 0.25*(dHeI(ix, iy, iz)+dHeII(ix, iy, iz)+dHeIII(ix, iy, iz));

          if (primordial_chemistry > 1) {
            mu(ix, iy, iz) += dHM(ix, iy, iz) + 0.5*(dH2I(ix, iy, iz)+dH2II(ix, iy, iz));
          }
          if (primordial_chemistry > 2) {
            mu(ix, iy, iz) += 0.5*(dDI(ix, iy, iz) + dDII(ix, iy, iz)) + dHDI(ix, iy, iz)/3.0;
          }
          mu(ix, iy, iz) /= d(ix, iy, iz);
        } else {
          mu(ix, iy, iz) = dflt_mu;
          // in an older version, mu = dflt_mu/d(ix, iy, iz), but I think that was a typo
        }
      }
    }
  }
}

void draw_stochastic(double &nsn_ii, double &nsn_ia) {
  // Draw Poisson-distributed integers with means nsn_ii and nsn_ia
  // and replace the input with the drawn values.
  // Input values are doubles because the feedback tables return doubles.
  
  // Use a thread-local RNG so repeated calls in the same thread are
  // not re-seeded each time.
  static thread_local std::mt19937_64 gen((std::random_device())());
  
  // Ensure non-negative means
  double mean_ii = (nsn_ii < 0.0) ? 0.0 : nsn_ii;
  double mean_ia = (nsn_ia < 0.0) ? 0.0 : nsn_ia;
  
  // The original Fortran version uses Knuth's algorithm.
  // Here we use std::poisson_distribution
  // which is typically more efficient.
  std::poisson_distribution<int> dist_ii(mean_ii);
  std::poisson_distribution<int> dist_ia(mean_ia);

  int k_ii = dist_ii(gen);
  int k_ia = dist_ia(gen);

  nsn_ii = static_cast<double>(k_ii);
  nsn_ia = static_cast<double>(k_ia);

}

void compute_snii_energy_momentum(
  double nsn,
  int distcells,
  const CelloView<const enzo_float,3> &d,  // pointer & values are const
  double n_avg, double Z_floor,
  int ic, int jc, int kc,
  double vol_cell, double vol_cell_oct,
  double mass_per_cell,
  double mom_per_cell[3][3][3],
  double energy_per_cell[3][3][3],
  double mom_mult
) {
  EnzoUnits * enzo_units = enzo::units();
  
  double dunit = enzo_units->density();
  double vunit = enzo_units->velocity();
  double lunit = enzo_units->length();
  double tunit = enzo_units->time();

  // Convert several Fortran expressions to C++.
  const double ergs_51_sqr = std::pow(10.0, 51.0/2.0); // matches Fortran 10**(51/2)
  const double Nnbors = distcells - 1; // 26 for 3^3
  const double betaSN = 1.0 / static_cast<double>(distcells); // 1/27
  const double Mej = mass_per_cell * static_cast<double>(distcells) * dunit * vol_cell; // cgs mass
  const double SN_energy_unit = std::pow(10.0,51.0); // 1e51 (cgs)
  const double SN_energy = (Nnbors > 0.0) ? (nsn * SN_energy_unit / Nnbors) : 0.0;

  double chi_th = 69.58 * std::pow(nsn, -2.0/17.0) *
                          std::pow(n_avg, -4.0/17.0) *
                          std::pow(Z_floor, -0.28);
  double d_fbck_cell = d(ic,jc,kc); // code density units

  // Loop neighbors
  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {

        const int ii = ic + i;
        const int jj = jc + j;
        const int kk = kc + k;

        // assume caller avoids boundaries (as in original Fortran)
        double den_cell = d(ii,jj,kk) * dunit; // cgs density
        double dMej = ((1.0 - betaSN) * Mej) / Nnbors;

        double dMswept = den_cell * vol_cell_oct
                       + ((1.0 - betaSN) * d_fbck_cell * dunit * vol_cell) / Nnbors
                       + dMej;
        double chi = dMswept / dMej;

        double pSN = 0.0;
        double energy_term = 0.0;

        if (chi < chi_th) {
          // resolved case
          double fe = 1.0 - ((chi - 1.0) / (3.0 * (chi_th - 1.0)));
          fe = std::max(0.0, fe);
          
          double tmp = 2.0 * chi * Mej * fe * nsn;
          if (tmp < 0.0) tmp = 0.0;  // Keep the safety copilot added during conversion
          pSN = mom_mult * std::sqrt(tmp) * ergs_51_sqr;  // cgs

          // kinetic energy from momentum (compare to SN_energy),
          double denom = den_cell * vol_cell + dMej;
          double ke_from_p = (denom > 0.0) ? 0.5 * std::pow(pSN / Nnbors, 2) / denom : 0.0;

          if (ke_from_p < SN_energy) {
            // add thermal energy (code units) to get up to 10^51 erg
            energy_term = (SN_energy - ke_from_p) / (dunit * vol_cell * vunit * vunit);
          } else {
            energy_term = 0.0;
          }
        } else {
          // unresolved / terminal momentum case
          double base = 0.0;
          if (nsn > 0.0 && n_avg > 0.0 && Z_floor > 0.0) {
            base = mom_mult * 3.0e10 * enzo_constants::mass_solar
                   * std::pow(nsn, 16.0/17.0)
                   * std::pow(n_avg, -2.0/17.0)
                   * std::pow(Z_floor, -0.14);
          }
          pSN = base;  // cgs
          energy_term = 0.0;
        }

        // Convert pSN from cgs (as constructed above) to code units:
        if (Nnbors > 0.0) {  // again, copilot is being very careful but I appreciate it
          pSN = pSN / (dunit * vol_cell * vunit * Nnbors);
        } else {
          pSN = 0.0;
        }

        // If center cell, zero injection (Fortran: dist == 0)
        int dist = i*i + j*j + k*k;
        if (dist == 0) {
          pSN = 0.0;
          energy_term = 0.0;
        }

        // store results at indices offset by +1 (range 0..2)
        const int ai = i + 1;
        const int aj = j + 1;
        const int ak = k + 1;

        mom_per_cell[ai][aj][ak] = pSN;             // code-mass * code-velocity / code-volume
        energy_per_cell[ai][aj][ak] = energy_term;  // code specific energy (velocity^2)
      }
    }
  }

  return;
}

void transform_momentum(
  const CelloView<const enzo_float,3> &d, // pointer & values are const
  const CelloView<enzo_float,3> &u,       // just pointer is const
  const CelloView<enzo_float,3> &v,
  const CelloView<enzo_float,3> &w,
  const int ic, const int jc, const int kc,
  double up, double vp, double wp,
  const int idir
) {
  // Converts velocity in the grid frame to momentum in the explosion frame
  // and back again:
  // idir == +1 : convert vel -> mom  : (u - up) * rho
  // idir == -1 : convert mom -> vel  : (u / rho) + up
  if (idir != +1 && idir != -1) {
    ERROR("transform_momentum","invalid idir (must be +1 or -1)");
  }

  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {
        int ii = ic + i;
        int jj = jc + j;
        int kk = kc + k;

        if (idir == +1) {
          // convert velocity -> momentum (explosion frame)
          u(ii,jj,kk) = (u(ii,jj,kk) - up) * d(ii,jj,kk);
          v(ii,jj,kk) = (v(ii,jj,kk) - vp) * d(ii,jj,kk);
          w(ii,jj,kk) = (w(ii,jj,kk) - wp) * d(ii,jj,kk);
        } else {
          // convert momentum -> velocity (grid frame)
          u(ii,jj,kk) = u(ii,jj,kk) / d(ii,jj,kk) + up;
          v(ii,jj,kk) = v(ii,jj,kk) / d(ii,jj,kk) + vp;
          w(ii,jj,kk) = w(ii,jj,kk) / d(ii,jj,kk) + wp;
        }
      }
    }
  }
  return;
}

void sum_energy(
  const CelloView<const enzo_float,3> &d,
  const CelloView<const enzo_float,3> &pu,  // pointer & values are const
  const CelloView<const enzo_float,3> &pv,
  const CelloView<const enzo_float,3> &pw,
  const int ic, const int jc, const int kc,
  double &kin_energy_sum
) {
  // Sum mass, kinetic energy and momentum magnitude over 3x3x3 cube
  kin_energy_sum = 0.0;

  for (int k = -1; k <= 1; ++k) {
    for (int j = -1; j <= 1; ++j) {
      for (int i = -1; i <= 1; ++i) {

        const int ii = ic + i;
        const int jj = jc + j;
        const int kk = kc + k;

        double mass_term = static_cast<double>(d(ii,jj,kk)); // mass/vol

        double pu2 = static_cast<double>(pu(ii,jj,kk));
        double pv2 = static_cast<double>(pv(ii,jj,kk));
        double pw2 = static_cast<double>(pw(ii,jj,kk));
        double mom_term_sq = pu2*pu2 + pv2*pv2 + pw2*pw2; // mass^2*velocity^2/volume^2

        double kin_energy = 0.0;
        if (mass_term > 0.0) {
          kin_energy = mom_term_sq / (2.0 * mass_term);
        } else {
          // avoid divide by zero; treat cell with zero mass as contributing nothing
          kin_energy = 0.0;
        }

        kin_energy_sum += kin_energy;
      }
    }
  }
  return;
}

void add_feedback_SNe(
  const CelloView<enzo_float,3> &pu,  // pointer can't change; values can
  const CelloView<enzo_float,3> &pv,
  const CelloView<enzo_float,3> &pw,
  const CelloView<enzo_float,3> &d,
  const CelloView<enzo_float,3> &ge,
  const CelloView<enzo_float,3> &te,
  const CelloView<enzo_float,3> &metals,
  const CelloView<enzo_float,3> &metalSNII,
  const CelloView<enzo_float,3> &metalSNIa,
  const int ic, const int jc, const int kc,
  const double mass_per_cell,
  const double mom_per_cell[3][3][3],
  const double mzeject, const double mzeject_ia, const double mzeject_ii,
  const double cell_volume,
  const bool track_metal_sources_,
  const bool cap_velocity_kick
) {
  // pu,pv,pw are momentum-like quantities (mass*vel/volume)
  // mass_per_cell in code density units (mass/volume)
  // mom_per_cell indexed [0..2] for i=-1..1
  // mzeject variables are metal density in code units (mass/volume)

  EnzoUnits * enzo_units = enzo::units();

  double dunit = enzo_units->density();
  double vunit = enzo_units->velocity();
  double lunit = enzo_units->length();

  const double to_km_per_s = vunit / 1e5;  // 1000 km/s in code units

  for (int i = -1; i <= 1; ++i) {
    for (int j = -1; j <= 1; ++j) {
      for (int k = -1; k <= 1; ++k) {

        int ii = ic + i;
        int jj = jc + j;
        int kk = kc + k;

        int dist = i*i + j*j + k*k;
        double mult = 0.0;
        if (dist == 0) {
          mult = 0.0;
        } else if (dist == 1) {
          mult = 1.0;
        } else if (dist == 2) {
          mult = 1.0 / std::sqrt(2.0);
        } else if (dist == 3) {
          mult = 1.0 / std::sqrt(3.0);
        }

        int ai = i + 1;  // 0..2
        int aj = j + 1;
        int ak = k + 1;

        double delta_pu = i * mult * mom_per_cell[ai][aj][ak];
        double delta_pv = j * mult * mom_per_cell[ai][aj][ak];
        double delta_pw = k * mult * mom_per_cell[ai][aj][ak];
        double delta_p = std::sqrt(
          delta_pu*delta_pu + delta_pv*delta_pv + delta_pw*delta_pw
        );

        double dens = static_cast<double>(d(ii,jj,kk));
        double new_dens = (dens + mass_per_cell);
        double vel_kick = 0.0;
        if (new_dens > 0.0) {
          vel_kick = delta_p / new_dens * to_km_per_s;
        }

        if (cap_velocity_kick && vel_kick > 1000.0) {
          delta_p = 1000.0 / to_km_per_s * new_dens;
          delta_pu = i * mult * delta_p;
          delta_pv = j * mult * delta_p;
          delta_pw = k * mult * delta_p;
        }

        // add momentum (pu,pv,pw are momentum-like quantities)
        double pu_start = static_cast<double>(pu(ii,jj,kk));
        double pv_start = static_cast<double>(pv(ii,jj,kk));
        double pw_start = static_cast<double>(pw(ii,jj,kk));

        pu(ii,jj,kk) = static_cast<enzo_float>(pu_start + delta_pu);
        pv(ii,jj,kk) = static_cast<enzo_float>(pv_start + delta_pv);
        pw(ii,jj,kk) = static_cast<enzo_float>(pw_start + delta_pw);

        // adjust thermal (specific) energies to account for added mass.
        // energies will be updated outside this function at the very end.
        double dratio = 1.0;
        dratio = dens / new_dens;

        te(ii,jj,kk) = static_cast<enzo_float>( static_cast<double>(te(ii,jj,kk)) * dratio );
        ge(ii,jj,kk) = static_cast<enzo_float>( static_cast<double>(ge(ii,jj,kk)) * dratio );

        if (!metals.is_null()) {
          double metal_old = static_cast<double>(metals(ii,jj,kk)); // metal density
          double metal_new = metal_old + mzeject;
          double mf_new = std::max(metal_new / new_dens, 0.95); // Cap metal fraction at 0.95
          metals(ii,jj,kk) = static_cast<enzo_float>(mf_new) * new_dens;
        }

        if (track_metal_sources_ && !metalSNII.is_null() && !metalSNIa.is_null()) {
          double mzii_old = static_cast<double>(metalSNII(ii,jj,kk));
          double mzia_old = static_cast<double>(metalSNIa(ii,jj,kk));

          double mzii_new = mzii_old + mzeject_ii;
          double mzia_new = mzia_old + mzeject_ia;

          double mfii_new = std::max(mzii_new / new_dens, 0.95);
          double mfia_new = std::max(mzia_new / new_dens, 0.95);

          metalSNII(ii,jj,kk) = static_cast<enzo_float>(mfii_new) * new_dens;
          metalSNIa(ii,jj,kk) = static_cast<enzo_float>(mfia_new) * new_dens;
        }

        // Finally add mass to density field
        d(ii,jj,kk) = static_cast<enzo_float>(dens + mass_per_cell);

      }
    }
  }
  return;
}

// =============================================================================

EnzoMethodFeedbackMechanical::EnzoMethodFeedbackMechanical(ParameterGroup p)
  : Method()
  , ir_feedback_(-1)
{
  cello::particle_descr()->check_particle_attribute("star", "initial_mass");

  FieldDescr * field_descr = cello::field_descr();
  EnzoUnits * enzo_units = enzo::units();

  ASSERT("EnzoMethodFeedbackMechanical::EnzoMethodFeedbackMechanical",
         "untested without dual-energy formalism",
         ! enzo::fluid_props()->dual_energy_config().is_disabled());
         
  // parameters
  stochastic_            = p.value_logical("stochastic_supernovae",false);
  pre_sne_               = p.value_logical("pre_sne",true);
  ejecta_mass_fraction_  = p.value_float("ejecta_mass_fraction",0.25);
  ejecta_metal_fraction_ = p.value_float("ejecta_metal_fraction",0.02);
  min_nsn_per_timestep_  = p.value_integer("nsn_per_timestep",1000);
  momentum_mult_         = p.value_float("momentum_multiplier",1.0);
  cap_velocity_kick_     = p.value_logical("cap_velocity_kick",true);
  track_metal_sources_   = p.value_logical("track_metal_sources",true);

  // required fields
  cello::define_field("density");
  cello::define_field("pressure");
  cello::define_field("total_energy");
  cello::define_field("internal_energy");
  cello::define_field("velocity_x");
  cello::define_field("velocity_y");
  cello::define_field("velocity_z");
  cello::define_field("metal_density");
  
  cello::define_field_in_group("metal_density","color");
  
  if (track_metal_sources_) {
    cello::define_field("metal_snII_density");
    cello::define_field("metal_snIa_density");

    cello::define_field_in_group("metal_snII_density","color");
    cello::define_field_in_group("metal_snIa_density","color");
  }

  // Initialize refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  // Initialize temporary fields
  i_yld_nsn = cello::field_descr()->insert_temporary();
  i_yld_mass = cello::field_descr()->insert_temporary();
  i_yld_metl = cello::field_descr()->insert_temporary();
  i_yld_vxp = cello::field_descr()->insert_temporary();
  i_yld_vyp = cello::field_descr()->insert_temporary();
  i_yld_vzp = cello::field_descr()->insert_temporary();
  
  i_d_dep  = cello::field_descr()->insert_temporary();
  i_te_dep = cello::field_descr()->insert_temporary();
  i_ge_dep = cello::field_descr()->insert_temporary();
  i_md_dep = cello::field_descr()->insert_temporary();
  i_vx_dep = cello::field_descr()->insert_temporary();
  i_vy_dep = cello::field_descr()->insert_temporary();
  i_vz_dep = cello::field_descr()->insert_temporary();
  
  i_d_dep_a  = cello::field_descr()->insert_temporary();
  i_te_dep_a = cello::field_descr()->insert_temporary();
  i_ge_dep_a = cello::field_descr()->insert_temporary();
  i_md_dep_a = cello::field_descr()->insert_temporary();
  i_vx_dep_a = cello::field_descr()->insert_temporary();
  i_vy_dep_a = cello::field_descr()->insert_temporary();
  i_vz_dep_a = cello::field_descr()->insert_temporary();
  
  if (track_metal_sources_) {
    i_yld_snii = cello::field_descr()->insert_temporary();
    i_yld_snia = cello::field_descr()->insert_temporary();

    i_mdii_dep = cello::field_descr()->insert_temporary();
    i_mdia_dep = cello::field_descr()->insert_temporary();
    i_mdii_dep_a = cello::field_descr()->insert_temporary();
    i_mdia_dep_a = cello::field_descr()->insert_temporary();
  }

  // Deposition across grid boundaries is handled using refresh with set_accumulate=true.
  // The set_accumulate flag tells Cello to include ghost zones in the refresh operation,
  // and adds the ghost zone values from the "src" field to the corresponding active zone
  // values in the "dst" field.
 
  // I'm currently using a set of two initially empty fields for each field that SNe directly affect.
  // FB deposition directly modifies the *_dep_* fields (including ghost zones). The ghost zone
  // values are then sent to the *_dep_*_a fields during the refresh operation. Values are then
  // copied back to the original field. 
  
  ir_feedback_ = add_refresh_();
  cello::simulation()->refresh_set_name(ir_feedback_,name()+":add");
  Refresh * refresh_fb = cello::refresh(ir_feedback_); 

  refresh_fb->set_accumulate(true);

  refresh_fb->add_field_src_dst(i_d_dep,  i_d_dep_a);
  refresh_fb->add_field_src_dst(i_te_dep, i_te_dep_a);
  refresh_fb->add_field_src_dst(i_ge_dep, i_ge_dep_a);
  refresh_fb->add_field_src_dst(i_md_dep, i_md_dep_a);
  refresh_fb->add_field_src_dst(i_vx_dep, i_vx_dep_a);
  refresh_fb->add_field_src_dst(i_vy_dep, i_vy_dep_a);
  refresh_fb->add_field_src_dst(i_vz_dep, i_vz_dep_a);
  
  if (track_metal_sources_) {
    refresh_fb->add_field_src_dst(i_mdii_dep, i_mdii_dep_a);
    refresh_fb->add_field_src_dst(i_mdia_dep, i_mdia_dep_a);
  } 

  // TODO add SNII and SNIa fields to p_method_feedback_mech_end callback
  refresh_fb->set_callback(CkIndex_EnzoBlock::p_method_feedback_mech_end());

  return;
}

void EnzoMethodFeedbackMechanical::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  // parameters
  p | stochastic_;
  p | pre_sne_;
  p | ejecta_mass_fraction_;
  p | ejecta_metal_fraction_;
  p | min_nsn_per_timestep_;
  p | momentum_mult_;
  p | cap_velocity_kick_;
  p | track_metal_sources_;

  // temporary yield fields
  p | i_yld_nsn;
  p | i_yld_mass;
  p | i_yld_metl;
  p | i_yld_snii;
  p | i_yld_snia;
  p | i_yld_vxp;
  p | i_yld_vyp;
  p | i_yld_vzp;

  // accumulation fields
  p | i_d_dep;
  p | i_te_dep;
  p | i_ge_dep;
  p | i_md_dep;
  p | i_mdii_dep;
  p | i_mdia_dep;
  p | i_vx_dep;
  p | i_vy_dep;
  p | i_vz_dep;

  p | i_d_dep_a;
  p | i_te_dep_a;
  p | i_ge_dep_a;
  p | i_md_dep_a;
  p | i_mdii_dep_a;
  p | i_mdia_dep_a;
  p | i_vx_dep_a;
  p | i_vy_dep_a;
  p | i_vz_dep_a;

  return;
}

void EnzoMethodFeedbackMechanical::compute (Block * block) throw()
{

  if (block->is_leaf()){
    this->compute_(block);
  }

  else {
    block->compute_done();
  }
  return;
}

void EnzoMethodFeedbackMechanical::compute_ (Block * block)
{
  
  //----------------------------------------------------
  // some constants here that might get moved to parameters or something else
  const float z_solar = enzo_constants::metallicity_solar; //Solar metal fraction (0.012)
  //-----------------------------------------------------

  EnzoBlock * enzo_block = enzo::block(block);
  Particle particle = enzo_block->data()->particle();
  EnzoUnits * enzo_units = enzo::units();

  double munit = enzo_units->mass();
  double lunit = enzo_units->length();
  double tunit = enzo_units->time();
  double dunit = enzo_units->density();

  double current_time  = block->time();

  Field field = enzo_block->data()->field();

  // Obtain grid sizes and ghost sizes
  int mx, my, mz, gx, gy, gz, nx, ny, nz;
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  field.size(&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);
  block->data()->lower(&xm,&ym,&zm);
  block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;

  double cell_volume = hx*hy*hz;  // TODO does this need to account for cosmology?
  double cell_volume_octant = cell_volume/8;
  int fb_cells = 27; // number of cells in feedback region (3^3 cube)

  const int rank = cello::rank();

  // apply feedback depending on particle type
  // for now, just do this for all star particles

  int it = particle.type_index("star");

  if (particle.num_particles(it) > 0){

    // get current field data
    CelloView<enzo_float,3> d  = field.view<enzo_float>("density");
    CelloView<enzo_float,3> vx = field.view<enzo_float>("velocity_x");
    CelloView<enzo_float,3> vy = field.view<enzo_float>("velocity_y");
    CelloView<enzo_float,3> vz = field.view<enzo_float>("velocity_z");
    CelloView<enzo_float,3> te = field.view<enzo_float>("total_energy");
    CelloView<enzo_float,3> ge = field.view<enzo_float>("internal_energy");
    CelloView<enzo_float,3> md = field.view<enzo_float>("metal_density");

    CelloView<enzo_float,3> mdii, mdia;  // uninitialized acts as nullptr
    if (track_metal_sources_) {
      mdii = field.view<enzo_float>("metal_snII_density");
      mdia = field.view<enzo_float>("metal_snIa_density");
    }

    // allocate temporary fields
    allocate_temporary_yields_(enzo_block);
    allocate_temporary_fluids_(enzo_block);

    // setup temporary mu field (allocated w/ fluids)
    compute_molecular_weight_(enzo_block);
    CelloView<enzo_float,3> mu = field.view<enzo_float>(i_mu);

    // initialize temporary fields as zero
    CelloView<enzo_float,3> nsn_yld  = field.view<enzo_float>(i_yld_nsn);
    CelloView<enzo_float,3> m_yld  = field.view<enzo_float>(i_yld_mass);
    CelloView<enzo_float,3> mz_yld = field.view<enzo_float>(i_yld_metl);
    CelloView<enzo_float,3> vxp_yld = field.view<enzo_float>(i_yld_vxp);
    CelloView<enzo_float,3> vyp_yld = field.view<enzo_float>(i_yld_vyp);
    CelloView<enzo_float,3> vzp_yld = field.view<enzo_float>(i_yld_vzp);
    
    CelloView<enzo_float,3> d_dep  = field.view<enzo_float>(i_d_dep);
    CelloView<enzo_float,3> te_dep = field.view<enzo_float>(i_te_dep);
    CelloView<enzo_float,3> ge_dep = field.view<enzo_float>(i_ge_dep);
    CelloView<enzo_float,3> md_dep = field.view<enzo_float>(i_md_dep);
    CelloView<enzo_float,3> vx_dep = field.view<enzo_float>(i_vx_dep);
    CelloView<enzo_float,3> vy_dep = field.view<enzo_float>(i_vy_dep);
    CelloView<enzo_float,3> vz_dep = field.view<enzo_float>(i_vz_dep);
    
    CelloView<enzo_float,3> d_dep_a  = field.view<enzo_float>(i_d_dep_a);
    CelloView<enzo_float,3> te_dep_a = field.view<enzo_float>(i_te_dep_a);
    CelloView<enzo_float,3> ge_dep_a = field.view<enzo_float>(i_ge_dep_a);
    CelloView<enzo_float,3> md_dep_a = field.view<enzo_float>(i_md_dep_a);
    CelloView<enzo_float,3> vx_dep_a = field.view<enzo_float>(i_vx_dep_a);
    CelloView<enzo_float,3> vy_dep_a = field.view<enzo_float>(i_vy_dep_a);
    CelloView<enzo_float,3> vz_dep_a = field.view<enzo_float>(i_vz_dep_a);
    
    CelloView<enzo_float,3> mzii_yld, mzia_yld;  // uninitialized acts as nullptr
    CelloView<enzo_float,3> mdii_dep, mdia_dep;
    CelloView<enzo_float,3> mdii_dep_a, mdia_dep_a;
    if (track_metal_sources_) {
      mzii_yld = field.view<enzo_float>(i_yld_snii);
      mzia_yld = field.view<enzo_float>(i_yld_snia);

      mdii_dep = field.view<enzo_float>(i_mdii_dep);
      mdia_dep = field.view<enzo_float>(i_mdia_dep);

      mdii_dep_a = field.view<enzo_float>(i_mdii_dep_a);
      mdia_dep_a = field.view<enzo_float>(i_mdia_dep_a);
    }

    // also initialize ghost zones
    for (int iz=0; iz<mz; iz++) {
      for (int iy=0; iy<my; iy++) {
        for (int ix=0; ix<my; ix++) {
          nsn_yld (ix, iy, iz) = 0.0;
          m_yld (ix, iy, iz) = 0.0;
          mz_yld(ix, iy, iz) = 0.0;
          vxp_yld(ix, iy, iz) = 0.0;
          vyp_yld(ix, iy, iz) = 0.0;
          vzp_yld(ix, iy, iz) = 0.0;
          
          d_dep (ix, iy, iz) = 0.0;
          te_dep(ix, iy, iz) = 0.0;
          ge_dep(ix, iy, iz) = 0.0;
          md_dep(ix, iy, iz) = 0.0;
          vx_dep(ix, iy, iz) = 0.0;
          vy_dep(ix, iy, iz) = 0.0;
          vz_dep(ix, iy, iz) = 0.0;
          
          d_dep_a (ix, iy, iz) = 0.0;
          te_dep_a(ix, iy, iz) = 0.0;
          ge_dep_a(ix, iy, iz) = 0.0;
          md_dep_a(ix, iy, iz) = 0.0;
          vx_dep_a(ix, iy, iz) = 0.0;
          vy_dep_a(ix, iy, iz) = 0.0;
          vz_dep_a(ix, iy, iz) = 0.0;

          if (track_metal_sources_) {
            mzii_yld(ix, iy, iz) = 0.0;
            mzia_yld(ix, iy, iz) = 0.0;

            mdii_dep(ix, iy, iz) = 0.0;
            mdia_dep(ix, iy, iz) = 0.0;

            mdii_dep_a(ix, iy, iz) = 0.0;
            mdia_dep_a(ix, iy, iz) = 0.0;
          }
        }
      }
    }

    // information for accessing particle data
    const int ia_m = particle.attribute_index (it, "mass");
    const int ia_x  = particle.attribute_index (it, "x");
    const int ia_y  = particle.attribute_index (it, "y");
    const int ia_z  = particle.attribute_index (it, "z");
    const int ia_vx = particle.attribute_index (it, "vx");
    const int ia_vy = particle.attribute_index (it, "vy");
    const int ia_vz = particle.attribute_index (it, "vz");
    const int ia_l = particle.attribute_index (it, "lifetime");
    const int ia_c = particle.attribute_index (it, "creation_time");
    const int ia_im = particle.attribute_index (it, "initial_mass");
    const int ia_mf = particle.attribute_index (it, "metal_fraction");

    const int dm = particle.stride(it, ia_m);
    const int dp = particle.stride(it, ia_x);
    const int dv = particle.stride(it, ia_vx);
    const int dl = particle.stride(it, ia_l);
    const int dc = particle.stride(it, ia_c);
    const int dim = particle.stride(it, ia_im);
    const int dmf = particle.stride(it, ia_mf);

    const int nb = particle.num_batches(it);

    // Iterate over particles to deposit their yield to a grid
    for (int ib=0; ib<nb; ib++){
      enzo_float *px=0, *py=0, *pz=0, *pvx=0, *pvy=0, *pvz=0;
      enzo_float *plifetime=0, *pcreation=0, *pmass=0, *pmetal=0, *pimass=0;

      // get current particle data
      pmass = (enzo_float *) particle.attribute_array(it, ia_m, ib);
      pmetal = (enzo_float *) particle.attribute_array(it, ia_mf, ib);
      pimass = (enzo_float *) particle.attribute_array(it, ia_im, ib);

      px  = (enzo_float *) particle.attribute_array(it, ia_x, ib);
      py  = (enzo_float *) particle.attribute_array(it, ia_y, ib);
      pz  = (enzo_float *) particle.attribute_array(it, ia_z, ib);
      pvx = (enzo_float *) particle.attribute_array(it, ia_vx, ib);
      pvy = (enzo_float *) particle.attribute_array(it, ia_vy, ib);
      pvz = (enzo_float *) particle.attribute_array(it, ia_vz, ib);

      plifetime = (enzo_float *) particle.attribute_array(it, ia_l, ib);
      pcreation = (enzo_float *) particle.attribute_array(it, ia_c, ib);

      int np = particle.num_particles(it,ib);

      // process each particle, calculating yields and adding that data to tmp grids
      for (int ip=0; ip<np; ip++){
        int ip_p = ip*dp; // pos
        int ip_m = ip*dm; // mass
        int ip_v = ip*dv; // velocity
        int ip_l = ip*dl; // lifetime
        int ip_c = ip*dc; // creation time
        int ip_mf = ip*dmf; // metallicity
        int ip_im = ip*dim; // initial mass

        if (pmass[ip_m] > 0.0 && plifetime[ip_l] > 0.0){
          double pmass_solar = pmass[ip_m] * munit/enzo_constants::mass_solar;
          const double age = (current_time - pcreation[ip_c]) * enzo_units->time() / enzo_constants::Myr_s;

          // Get yields from tables
          double nsn, nsn_ii, nsn_ia, mej, mej_ii, mej_ia,
            mzej, mzej_ii, mzej_ia;
          // TODO ensure metal yields are returned as densities

          // Stochasitcally sample a Poisson distribution using the expected
          // number of SNe as the mean.
          if (stochastic_){
            double nsn_ii_sto, nsn_ia_sto;
            nsn_ii_sto = nsn_ii;
            nsn_ia_sto = nsn_ia;
            draw_stochastic(nsn_ii_sto, nsn_ia_sto);

            // Adjust mass & metal yields based on the drawn number of SNe.
            // Convert to an ejection fraction of the initial particle mass
            // and then scale to the current particle mass.
            // User ejecta fraction parameters serve as a limit
            if (nsn_ii_sto > 0.0) {
              mej_ii = std::min(
                (mej_ii/nsn_ii * nsn_ii_sto) * pmass[ip_m]/pimass[ip_im],
                ejecta_mass_fraction_ * pmass[ip_m]
              );
              mzej_ii = std::min(
                (mzej_ii/nsn_ii * nsn_ii_sto) * pmass[ip_m]/pimass[ip_im],
                ejecta_metal_fraction_ * pmass[ip_m]  // TODO ensure particle mass is a density
              );
            } else {
              mej_ii = 0.0;
              mzej_ii = 0.0;
            }

            if (nsn_ia_sto > 0.0) {
              mej_ia = std::min(
                (mej_ia/nsn_ia * nsn_ia_sto) * pmass[ip_m]/pimass[ip_im],
                ejecta_mass_fraction_ * pmass[ip_m]
              );
              mzej_ia = std::min(
                (mzej_ia/nsn_ia * nsn_ia_sto) * pmass[ip_m]/pimass[ip_im],
                ejecta_metal_fraction_ * pmass[ip_m]
              );
            } else {
              mej_ia = 0.0;
              mzej_ia = 0.0;
            }

            // reset SNe numbers to those that were drawn
            nsn_ii = nsn_ii_sto;
            nsn_ia = nsn_ia_sto;
          } // end stochastic

          // sum yields
          nsn = nsn_ii + nsn_ia;
          if (nsn < min_nsn_per_timestep_)
            nsn = 0.0;

          mej = mej_ii + mej_ia;
          mzej = mzej_ii + mzej_ia;

          // Calculate energy yield
          double energy = nsn * 1.0e51; // in ergs
          energy = energy / (munit * lunit*lunit / (tunit*tunit)); // code units
          energy = energy / cell_volume; // energy per unit volume

          // Subtract mass from particle
          pmass[ip_m] -= mej * enzo_constants::mass_solar / munit;
          
          // Get position of particle in space
          double x = px[ip_p];
          double y = py[ip_p];
          double z = pz[ip_p];

          // compute index of central feedback cell
          // this must account for ghost zones 
          // as xm, ym, zm reference the active zone
          int ix = (int) floor((x - xm) / hx + gx);
          int iy = (int) floor((y - ym) / hy + gy);
          int iz = (int) floor((z - zm) / hz + gz);

          // Store yields in temporary grids
          // so we can sum over all particles a given cell
          nsn_yld(ix, iy, iz) += (enzo_float) nsn;
          m_yld(ix, iy, iz) += (enzo_float) mej;
          mz_yld(ix, iy, iz) += (enzo_float) mzej;
          if (track_metal_sources_) {
            mzii_yld(ix, iy, iz) += (enzo_float) mzej_ii;
            mzia_yld(ix, iy, iz) += (enzo_float) mzej_ia;
          }

          // Weight avg velocity by number of SNe going off
          vxp_yld(ix, iy, iz) += (enzo_float) (pvx[ip_v] * nsn);
          vyp_yld(ix, iy, iz) += (enzo_float) (pvy[ip_v] * nsn);
          vzp_yld(ix, iy, iz) += (enzo_float) (pvz[ip_v] * nsn);
          
        } // if mass and lifetime > 0
      } // end particle loop
    } // end batch loop

    // Now, iterate over yield grids and apply to actual fields
    // using 3x3x3 dummy grids to capture momentum cancellation
    double m_per_cell, mz_per_cell;
    double mzii_per_cell=0.0, mzia_per_cell=0.0;
    double mom_per_cell[3][3][3], eng_per_cell[3][3][3], ke_old_grid[3][3][3];

    CelloView<enzo_float, 3> px1(3,3,3), py1(3,3,3), pz1(3,3,3);
    CelloView<enzo_float, 3> d1(3,3,3), ge1(3,3,3), te1(3,3,3);
    CelloView<enzo_float, 3> null; // uninitialized CelloView analogous to nullptr

    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
        for (int ix=gx; ix<my-gx; ix++) {
          if (nsn_yld(ix, iy, iz) == 0.0) continue; // no SNe from this cell

          // compute average density & metallicity around this cell
          double avg_Z = 0.0, avg_n = 0.0;
          double cell_mf, cell_Z, cell_n;
          for (int k=-1; k<=1; k++){
            for (int j=-1; j<=1; j++){
              for (int i=-1; i<=1; i++){
                cell_mf = md(ix+i, iy+j, iz+k) / d(ix+i, iy+j, iz+k);
                cell_Z += cell_mf / z_solar;
                cell_n += d(ix+i, iy+j, iz+k)*dunit 
                        * mu(ix+i, iy+j, iz+k)/enzo_constants::mass_hydrogen;
              }
            }
          }
          // The Kimm & Cen (2014) scheme (eq A5) requires a floor on Z_avg
          avg_Z = std::max(avg_Z/fb_cells, 0.01);
          avg_n /= fb_cells;

          // Compute mass, momentum, & energy this cell is responsible for injecting
          m_per_cell = m_yld(ix, iy, iz) / fb_cells;
          mz_per_cell = mz_yld(ix, iy, iz) / fb_cells;
          if (track_metal_sources_) {  // these are initialized to zero
            mzii_per_cell = mzii_yld(ix, iy, iz) / fb_cells;
            mzia_per_cell = mzia_yld(ix, iy, iz) / fb_cells;
          }

          compute_snii_energy_momentum(
            nsn_yld(ix, iy, iz), fb_cells,
            d, avg_n, avg_Z,
            ix, iy, iz,
            cell_volume, cell_volume_octant,
            m_per_cell, mom_per_cell, eng_per_cell,
            momentum_mult_
          );

          // check that momentum & energy injection are nonzero
          // this is the only time I'll use ijk insead of kji ordering
          // because I'm only working with an array of doubles
          double mom_inj = 0.0, eng_inj = 0.0;
          for (int i=0; i<=2; i++){
            for (int j=0; j<=2; j++){
              for (int k=0; k<=2; k++){
                mom_inj += mom_per_cell[i][j][k];
                eng_inj += eng_per_cell[i][j][k];
              }
            }
          }
          if ((mom_inj <= 0.0) || (eng_inj <= 0.0))
            continue;

          // find initial kinetic energy in grid frame before injection
          // units of mass*velocity^2/volume
          for (int k=-1; k<=1; k++)
            for (int j=-1; j<=1; j++)
              for (int i=-1; i<=1; i++)
                ke_old_grid[i+1][j+1][k+1] = d(ix+i, iy+j, iz+k) * 0.5 * (
                  vx(ix+i, iy+j, iz+k)*vx(ix+i, iy+j, iz+k) +
                  vy(ix+i, iy+j, iz+k)*vy(ix+i, iy+j, iz+k) +
                  vz(ix+i, iy+j, iz+k)*vz(ix+i, iy+j, iz+k)
                );

          // convert current velocities to momenta and transform into the frame
          // comoving with explosion-weighted average velocity in this cell
          double vxp_avg = vxp_yld(ix, iy, iz)/nsn_yld(ix, iy, iz);
          double vyp_avg = vyp_yld(ix, iy, iz)/nsn_yld(ix, iy, iz);
          double vzp_avg = vzp_yld(ix, iy, iz)/nsn_yld(ix, iy, iz);
          transform_momentum(
            d, vx, vy, vz,
            ix, iy, iz,
            vxp_avg, vyp_avg, vzp_avg,
            1
          );

          // sum kintetic energy before injection.
          double ke_old_explosion;
          sum_energy(
            d, vx, vy, vz,
            ix, iy, iz,
            ke_old_explosion);

          // Zero dummy grids
          // use kji for faster CelloView access
          for (int k=0; k<=2; k++){
            for (int j=0; j<=2; j++){
              for (int i=0; i<=2; i++){

                px1(i,j,k) = 0.0;
                py1(i,j,k) = 0.0;
                pz1(i,j,k) = 0.0;

                // copy host-grid density into dummy density array
                d1(i,j,k) = d(ix + i-1, iy + j-1, iz + k-1);

                ge1(i,j,k) = 0.0;
                te1(i,j,k) = 0.0;
              }
            }
          }

          // Add FB to dummy grids
          add_feedback_SNe(
            px1, py1, pz1,
            d1, ge1, te1,
            null, null, null, // no metals in dummy
            1, 1, 1,          // center of dummy grid (runs 0..2)
            m_per_cell,
            mom_per_cell,
            0.0, 0.0, 0.0,    // no metals in dummy
            cell_volume,
            track_metal_sources_,
            cap_velocity_kick_
          );

          // sum mass, KE, and momentum in particle frame (dummy grids)
          double ke_dummy_explosion;
          sum_energy(
            d1, px1, py1, pz1,
            1, 1, 1,  // center of dummy grid
            ke_dummy_explosion
          );

          // Add FB to actual grids
          add_feedback_SNe(
            vx_dep, vy_dep, vz_dep,
            d_dep, ge_dep, te_dep,
            md_dep, mdii_dep, mdia_dep,
            ix, iy, iz,
            m_per_cell,
            mom_per_cell,
            mz_per_cell,
            mzii_per_cell,
            mzia_per_cell,
            cell_volume,
            track_metal_sources_,
            cap_velocity_kick_
          );

          // sum mass, KE, and momentum after injection
          double ke_new_explosion;
          sum_energy(
            d, vx, vy, vz,
            ix, iy, iz,
            ke_new_explosion
          );

          // The KE added to real grid can be less than added to the dummy
          // because of momentum cancellation
          // If this is true, add the missing energy as thermal energy
          double ke_injected = ke_new_explosion - ke_old_explosion;
          double ke_deficit = ke_dummy_explosion - ke_injected;
          ke_deficit = std::max(ke_deficit, 0.0); // must be positive
          double ke_deficit_per_cell = ke_deficit / fb_cells;

          // Convert momenta back to velocities & transform to simulation frame
          transform_momentum(
            d, vx, vy, vz,
            ix, iy, iz,
            vxp_avg, vyp_avg, vzp_avg,
            -1
          );

          // Track change in KE in the total energy field
          // and inject both energy lost due to momentum cancellation (ke_deficit_per_cell)
          // and energy from resolved Sedov-Taylor phase (eng_per_cell) as thermal energy.
          // Note that we have to recalculate the KE because of the frame change.
          for (int k=-1; k<=1; k++){
            for (int j=-1; j<=1; j++){
              for (int i=-1; i<=1; i++){

                double ke_new_grid = 0.5 * d(ix+i, iy+j, iz+k) * (
                  vx(ix+i, iy+j, iz+k)*vx(ix+i, iy+j, iz+k) +
                  vy(ix+i, iy+j, iz+k)*vy(ix+i, iy+j, iz+k) +
                  vz(ix+i, iy+j, iz+k)*vz(ix+i, iy+j, iz+k)
                );
                double ke_added = ke_new_grid - ke_old_grid[i+1][j+1][k+1];
                double ge_added = ke_deficit_per_cell + eng_per_cell[i+1][j+1][k+1];

                // add, converting from energy density to specific energy
                te_dep(ix+i, iy+j, iz+k) += (ke_added+ge_added) / d(ix+i, iy+j, iz+k);
                ge_dep(ix+i, iy+j, iz+k) += ge_added / d(ix+i, iy+j, iz+k);
              } // end loop over local z
            }
          }
        } // end loop over grid z
      }
    }

    // Deallocate yield fields; keep fluid fields for the refresh
    deallocate_temporary_yields_(enzo_block);

    // refresh
    cello::refresh(ir_feedback_)->set_active(enzo_block->is_leaf());
    enzo_block->refresh_start(ir_feedback_, CkIndex_EnzoBlock::p_method_feedback_mech_end());
    return;
  } // end check for particles
}

void EnzoBlock::p_method_feedback_mech_end() 
{  
  EnzoMethodFeedbackMechanical * method = static_cast<EnzoMethodFeedbackMechanical*> (this->method());
  method->add_accumulate_fields(this);
  
  compute_done();
  return;
}

void EnzoMethodFeedbackMechanical::add_accumulate_fields(EnzoBlock * enzo_block) throw()
{
  Field field = enzo_block->data()->field();

  int mx, my, mz, gx, gy, gz, nx, ny, nz;
  double xm, ym, zm, xp, yp, zp, hx, hy, hz;
  field.size(&nx,&ny,&nz);
  field.ghost_depth(0,&gx,&gy,&gz);
  enzo_block->data()->lower(&xm,&ym,&zm);
  enzo_block->data()->upper(&xp,&yp,&zp);
  field.cell_width(xm,xp,&hx,ym,yp,&hy,zm,zp,&hz);

  mx = nx + 2*gx;
  my = ny + 2*gy;
  mz = nz + 2*gz;

  // add accumulated values over and reset them to zero

  CelloView<enzo_float,3> d  = field.view<enzo_float>("density");
  CelloView<enzo_float,3> vx = field.view<enzo_float>("velocity_x");
  CelloView<enzo_float,3> vy = field.view<enzo_float>("velocity_y");
  CelloView<enzo_float,3> vz = field.view<enzo_float>("velocity_z");
  CelloView<enzo_float,3> te = field.view<enzo_float>("total_energy");
  CelloView<enzo_float,3> ge = field.view<enzo_float>("internal_energy");
  CelloView<enzo_float,3> md = field.view<enzo_float>("metal_density");

  CelloView<enzo_float,3> mdii, mdia;  // uninitialized acts as nullptr
  if (track_metal_sources_) {
    mdii = field.view<enzo_float>("metal_snII_density");
    mdia = field.view<enzo_float>("metal_snIa_density");
  }

  CelloView<enzo_float,3> d_dep_a  = field.view<enzo_float>(i_d_dep_a);
  CelloView<enzo_float,3> vx_dep_a = field.view<enzo_float>(i_vx_dep_a);
  CelloView<enzo_float,3> vy_dep_a = field.view<enzo_float>(i_vy_dep_a);
  CelloView<enzo_float,3> vz_dep_a = field.view<enzo_float>(i_vz_dep_a);
  CelloView<enzo_float,3> te_dep_a = field.view<enzo_float>(i_te_dep_a);
  CelloView<enzo_float,3> ge_dep_a = field.view<enzo_float>(i_ge_dep_a);
  CelloView<enzo_float,3> md_dep_a = field.view<enzo_float>(i_md_dep_a);

  CelloView<enzo_float,3> mdii_dep_a, mdia_dep_a;  // uninitialized acts as nullptr
  if (track_metal_sources_) {
    mdii_dep_a = field.view<enzo_float>(i_mdii_dep_a);
    mdia_dep_a = field.view<enzo_float>(i_mdia_dep_a);
  }

  EnzoUnits * enzo_units = enzo::units();
  double cell_volume_code = hx*hy*hz;
  double cell_volume_cgs = cell_volume_code * enzo_units->volume();
  double rhounit = enzo_units->density();

  double maxEvacFraction = 0.75; // TODO: make this a parameter

  // multiply by this value to convert cell density (in code units)
  // to cell mass in Msun
  double rho_to_m = rhounit*cell_volume_cgs / enzo_constants::mass_solar;

  for (int iz=gz; iz<nz+gz; iz++){
    for (int iy=gy; iy<ny+gy; iy++){
      for (int ix=gx; ix<nx+gx; ix++){

        if (d_dep_a(ix,iy,iz) != 0) { // if any deposition

          double d_old = d(ix,iy,iz);
          
          // Could have a race condition here where if one particle updates the density of a cell,
          // that update won't get communicated to the other block until the end of the cycle.
          // In rare cases, this could result in a cell having a negative density because "centralMass"
          // and "centralMetals" in deposit_feedback() will be overpredicted going into the negative-mass
          // CiC for clearing out gas from the central cell. Catch this case by setting density
          // to (1-maxEvacFraction) * density[i] and metal_density to (1-maxEvacFraction) * metal_density[i] if
          // either go negative.

          if (d(ix,iy,iz) + d_dep_a(ix,iy,iz) < 0) {
            d(ix,iy,iz) *= 1-maxEvacFraction;
          }
          else {
            d(ix,iy,iz) += d_dep_a(ix,iy,iz);
          }

          if (md(ix,iy,iz) + md_dep_a(ix,iy,iz) < 0) {
            md(ix,iy,iz) *= 1-maxEvacFraction;
          }
          else {
            md(ix,iy,iz) += md_dep_a(ix,iy,iz);
          }

          vx(ix,iy,iz) += vx_dep_a(ix,iy,iz);
          vy(ix,iy,iz) += vy_dep_a(ix,iy,iz);
          vz(ix,iy,iz) += vz_dep_a(ix,iy,iz);
          
          double d_new = d(ix,iy,iz);
          double d_ratio = d_old/d_new;
          // energies need to be rescaled by the new density after adding (at the old density)
          ge(ix,iy,iz) = (ge(ix,iy,iz) + std::max( (double) ge_dep_a(ix,iy,iz), 0.0)) * d_ratio;
          te(ix,iy,iz) = (te(ix,iy,iz) + std::max( (double) te_dep_a(ix,iy,iz), 0.0)) * d_ratio; 
          
          // rescale color fields to account for new densities
          //EnzoMethodStarMaker::rescale_densities(enzo_block, i, M_scale_tot);
          // undo rescaling of metal_density field
          md(ix,iy,iz) *= d_ratio;

          

         }
      }
    }
  }

  deallocate_temporary_fluids_(enzo_block);
  return;
}

// ----------------------------------------------------------------------------

double EnzoMethodFeedbackMechanical::timestep (Block * block) throw()
{
  // In general this is not needed, but could imagine putting timestep
  // limiters in situations where, for example, one would want
  // dt < star_lifetime (or something like that), especially if
  // important things happen throughout the star's lifetime.
  EnzoUnits * enzo_units = enzo::units();
  
  // return 1000.0 * enzo_constants::yr_s / enzo_units->time();
  return std::numeric_limits<double>::max();
}