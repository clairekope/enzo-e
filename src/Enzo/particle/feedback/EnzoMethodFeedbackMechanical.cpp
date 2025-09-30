
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

  // Initialize refresh object
  cello::simulation()->refresh_set_name(ir_post_,name());
  Refresh * refresh = cello::refresh(ir_post_);
  refresh->add_all_fields();

  stochastic_            = p.value_logical("supernovae",true);
  ejecta_mass_fraction_  = p.value_float("ejecta_mass_fraction",0.25);
  ejecta_metal_fraction_ = p.value_float("ejecta_metal_fraction",0.02);
  min_nsn_per_timestep_  = p.value_integer("nsn_per_timestep",1000);

  // Initialize temporary fields
  i_dep_mass  = cello::field_descr()->insert_temporary();
  i_dep_metl = cello::field_descr()->insert_temporary();
  i_dep_snii = cello::field_descr()->insert_temporary();
  i_dep_snia = cello::field_descr()->insert_temporary();
  i_dep_px = cello::field_descr()->insert_temporary();
  i_dep_py = cello::field_descr()->insert_temporary();
  i_dep_pz = cello::field_descr()->insert_temporary();
  i_nsn = cello::field_descr()->insert_temporary();

  i_nsn_a  = cello::field_descr()->insert_temporary();
  i_dep_mass_a = cello::field_descr()->insert_temporary();
  i_dep_metl_a = cello::field_descr()->insert_temporary();
  i_dep_snii_a = cello::field_descr()->insert_temporary();
  i_dep_snia_a = cello::field_descr()->insert_temporary();
  i_dep_px_a = cello::field_descr()->insert_temporary();
  i_dep_py_a = cello::field_descr()->insert_temporary();
  i_dep_pz_a= cello::field_descr()->insert_temporary();

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

  refresh_fb->add_field_src_dst(i_dep_mass, i_dep_mass_a);
  refresh_fb->add_field_src_dst(i_dep_metl, i_dep_metl_a);
  refresh_fb->add_field_src_dst(i_dep_snii, i_dep_snii_a);
  refresh_fb->add_field_src_dst(i_dep_snia, i_dep_snia_a);
  refresh_fb->add_field_src_dst(i_dep_px, i_dep_px_a);
  refresh_fb->add_field_src_dst(i_dep_py, i_dep_py_a);
  refresh_fb->add_field_src_dst(i_dep_pz, i_dep_pz_a);
  refresh_fb->add_field_src_dst(i_nsn, i_nsn_a);
 
  refresh_fb->set_callback(CkIndex_EnzoBlock::p_method_feedback_mech_end());

  return;
}

void EnzoMethodFeedbackMechanical::pup (PUP::er &p)
{
  /// NOTE: Change this function whenever attributes change

  TRACEPUP;

  Method::pup(p);

  p | stochastic_;
  p | ejecta_mass_fraction_;
  p | ejecta_metal_fraction_;
  p | min_nsn_per_timestep_;

  p | i_nsn;
  p | i_nsn_a;
  p | i_dep_mass;
  p | i_dep_mass_a;
  p | i_dep_metl;
  p | i_dep_metl_a;
  p | i_dep_snii;
  p | i_dep_snii_a; 
  p | i_dep_snia;
  p | i_dep_snia_a;
  p | i_dep_px;
  p | i_dep_px_a;
  p | i_dep_py;
  p | i_dep_py_a;
  p | i_dep_pz;
  p | i_dep_pz_a;

  return;
}

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

void EnzoBlock::p_method_feedback_mech_end() 
{  
  EnzoMethodFeedbackMechanical * method = static_cast<EnzoMethodFeedbackMechanical*> (this->method());
  //method->add_accumulate_fields(this);
  
  compute_done();
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

  const int rank = cello::rank();

  // apply feedback depending on particle type
  // for now, just do this for all star particles

  int it = particle.type_index("star");

  if (particle.num_particles(it) > 0){

    // get current field data
    CelloView<enzo_float,3> d  = field.view<enzo_float>("density");
    CelloView<enzo_float,3> te = field.view<enzo_float>("total_energy");
    CelloView<enzo_float,3> ge = field.view<enzo_float>("internal_energy");
    CelloView<enzo_float,3> mf = field.view<enzo_float>("metal_density");

    allocate_temporary_(enzo_block);

    // initialize temporary fields as zero
    CelloView<enzo_float,3> n_dep   = field.view<enzo_float>(i_nsn);
    CelloView<enzo_float,3> n_dep_a = field.view<enzo_float>(i_nsn_a);

    CelloView<enzo_float,3> m_dep  = field.view<enzo_float>(i_dep_mass);
    CelloView<enzo_float,3> mz_dep = field.view<enzo_float>(i_dep_metl);
    CelloView<enzo_float,3> mzii_dep = field.view<enzo_float>(i_dep_snii);
    CelloView<enzo_float,3> mzia_dep = field.view<enzo_float>(i_dep_snia);
    CelloView<enzo_float,3> px_dep = field.view<enzo_float>(i_dep_px);
    CelloView<enzo_float,3> py_dep = field.view<enzo_float>(i_dep_py);
    CelloView<enzo_float,3> pz_dep = field.view<enzo_float>(i_dep_pz);

    CelloView<enzo_float,3> m_dep_a  = field.view<enzo_float>(i_dep_mass_a);
    CelloView<enzo_float,3> mz_dep_a = field.view<enzo_float>(i_dep_metl_a);
    CelloView<enzo_float,3> mzii_dep_a = field.view<enzo_float>(i_dep_snii_a);
    CelloView<enzo_float,3> mzia_dep_a = field.view<enzo_float>(i_dep_snia_a);
    CelloView<enzo_float,3> px_dep_a = field.view<enzo_float>(i_dep_px_a);
    CelloView<enzo_float,3> py_dep_a = field.view<enzo_float>(i_dep_py_a);
    CelloView<enzo_float,3> pz_dep_a = field.view<enzo_float>(i_dep_pz_a);

    for (int iz=gz; iz<mz-gz; iz++) {
      for (int iy=gy; iy<my-gy; iy++) {
        for (int ix=gx; ix<my-gx; ix++) {
          n_dep(ix, iy, iz) = 0.0;
          m_dep (ix, iy, iz) = 0.0;
          mz_dep(ix, iy, iz) = 0.0;
          mzii_dep(ix, iy, iz) = 0.0;
          mzia_dep(ix, iy, iz) = 0.0;
          px_dep(ix, iy, iz) = 0.0;
          py_dep(ix, iy, iz) = 0.0;
          pz_dep(ix, iy, iz) = 0.0;
          
          n_dep_a(ix, iy, iz) = 0.0;
          m_dep_a(ix, iy, iz) = 0.0;
          mz_dep_a(ix, iy, iz) = 0.0;
          mzii_dep_a(ix, iy, iz) = 0.0;
          mzia_dep_a(ix, iy, iz) = 0.0;
          px_dep_a(ix, iy, iz) = 0.0;
          py_dep_a(ix, iy, iz) = 0.0;
          pz_dep_a(ix, iy, iz) = 0.0;
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
                ejecta_metal_fraction_ * pmass[ip_m]
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

            nsn_ii = nsn_ii_sto;
            nsn_ia = nsn_ia_sto;
          } // end stochastic

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
          n_dep(ix, iy, iz) += (enzo_float) nsn;
          m_dep(ix, iy, iz) += (enzo_float) mej;
          mz_dep(ix, iy, iz) += (enzo_float) mzej;
          mzii_dep(ix, iy, iz) += (enzo_float) mzej_ii;
          mzia_dep(ix, iy, iz) += (enzo_float) mzej_ia;

          // Weight avg velocity by number of SNe going off
          px_dep(ix, iy, iz) += (enzo_float) (pvx[ip_v] * nsn);
          py_dep(ix, iy, iz) += (enzo_float) (pvy[ip_v] * nsn);
          pz_dep(ix, iy, iz) += (enzo_float) (pvz[ip_v] * nsn);
          
        } // if mass and lifetime > 0
      } // end particle loop
    } // end batch loop

    // refresh
    cello::refresh(ir_feedback_)->set_active(enzo_block->is_leaf());
    enzo_block->refresh_start(ir_feedback_, CkIndex_EnzoBlock::p_method_feedback_mech_end());
    return;
  } // end check for particles
}


// ----------------------------------------------------------------------------

