// See LICENSE_CELLO file for license and copyright information

/// @file     enzo_EnzoConfig.cpp
/// @author   James Bordner (jobordner@ucsd.edu)
/// @date     2012-10-03
/// @brief    Implementation of the EnzoConfig class

#include "cello.hpp"
#include "enzo.hpp"

extern CProxy_EnzoSimulation proxy_enzo_simulation;

//----------------------------------------------------------------------

EnzoConfig g_enzo_config;

EnzoConfig::EnzoConfig() throw ()
  :
  adapt_mass_type(0),
  ppm_diffusion(false),
  ppm_dual_energy(false),
  ppm_dual_energy_eta_1(0.0),
  ppm_dual_energy_eta_2(0.0),
  ppm_flattening(0),
  ppm_minimum_pressure_support_parameter(0),
  ppm_number_density_floor(0.0),
  ppm_density_floor(0.0),
  ppm_pressure_floor(0.0),
  ppm_pressure_free(false),
  ppm_temperature_floor(0.0),
  ppm_steepening(false),
  ppm_use_minimum_pressure_support(false),
  ppm_mol_weight(0.0),
  field_gamma(0.0),
  physics_cosmology(false),
  physics_cosmology_hubble_constant_now(0.0),
  physics_cosmology_omega_matter_now(0.0),
  physics_cosmology_omega_lamda_now(0.0),
  physics_cosmology_omega_baryon_now(1.0),
  physics_cosmology_omega_cdm_now(0.0),
  physics_cosmology_comoving_box_size(0.0),
  physics_cosmology_max_expansion_rate(0.0),
  physics_cosmology_initial_redshift(0.0),
  physics_cosmology_final_redshift(0.0),
  physics_gravity(false),
  // EnzoInitialBCenter
  initial_bcenter_update_etot(false),
  // EnzoInitialCloud
  initial_cloud_subsample_n(0),
  initial_cloud_radius(0.),
  initial_cloud_center_x(0.0),
  initial_cloud_center_y(0.0),
  initial_cloud_center_z(0.0),
  initial_cloud_density_cloud(0.0),
  initial_cloud_density_wind(0.0),
  initial_cloud_velocity_wind(0.0),
  initial_cloud_etot_wind(0.0),
  initial_cloud_eint_wind(0.0),
  initial_cloud_metal_mass_frac(0.0),
  initial_cloud_initialize_uniform_bfield(false),
  initial_cloud_perturb_stddev(0.0),
  initial_cloud_trunc_dev(0.0),
  initial_cloud_perturb_seed(0),
  // EnzoInitialCosmology
  initial_cosmology_temperature(0.0),
  // EnzoInitialCollapse
  initial_collapse_rank(0),
  initial_collapse_radius_relative(0.0),
  initial_collapse_particle_ratio(0.0),
  initial_collapse_mass(0.0),
  initial_collapse_temperature(0.0),
  // EnzoInitialGrackleTest
#ifdef CONFIG_USE_GRACKLE
  initial_grackle_test_maximum_H_number_density(1000.0),
  initial_grackle_test_maximum_metallicity(1.0),
  initial_grackle_test_maximum_temperature(1.0E8),
  initial_grackle_test_minimum_H_number_density(0.1),
  initial_grackle_test_minimum_metallicity(1.0E-4),
  initial_grackle_test_minimum_temperature(10.0),
  initial_grackle_test_reset_energies(0),
#endif /* CONFIG_USE_GRACKLE */
  // EnzoInitialInclinedWave
  initial_inclinedwave_alpha(0.0),
  initial_inclinedwave_beta(0.0),
  initial_inclinedwave_amplitude(0.0),
  initial_inclinedwave_lambda(0.0),
  initial_inclinedwave_parallel_vel(std::numeric_limits<double>::min()),
  initial_inclinedwave_positive_vel(true),
  initial_inclinedwave_wave_type(""),
  // EnzoInitialMusic
  initial_music_field_files(),
  initial_music_field_datasets(),
  initial_music_field_names(),
  initial_music_field_coords(),
  initial_music_particle_files(),
  initial_music_particle_datasets(),
  initial_music_particle_coords(),
  initial_music_particle_types(),
  initial_music_particle_attributes(),
  initial_music_throttle_internode(),
  initial_music_throttle_intranode(),
  initial_music_throttle_node_files(),
  initial_music_throttle_close_count(),
  initial_music_throttle_group_size(),
  initial_music_throttle_seconds_stagger(),
  initial_music_throttle_seconds_delay(),
  // EnzoInitialPm
  initial_pm_field(""),
  initial_pm_mpp(0.0),
  initial_pm_level(0),
  // EnzoInitialSedov[23]
  initial_sedov_rank(0),
  initial_sedov_radius_relative(0.0),
  initial_sedov_pressure_in(0.0),
  initial_sedov_pressure_out(0.0),
  initial_sedov_density(0.0),
  // EnzoInitialSedovRandom
  initial_sedov_random_half_empty(false),
  initial_sedov_random_grackle_cooling(false),
  initial_sedov_random_max_blasts(0),
  initial_sedov_random_radius_relative(0.0),
  initial_sedov_random_pressure_in(0.0),
  initial_sedov_random_pressure_out(0.0),
  initial_sedov_random_density(0.0),
  initial_sedov_random_te_multiplier(0),
  // EnzoInitialShockTube
  initial_shock_tube_setup_name(""),
  initial_shock_tube_aligned_ax(""),
  initial_shock_tube_axis_velocity(0.0),
  initial_shock_tube_trans_velocity(0.0),
  initial_shock_tube_flip_initialize(false),
  // EnzoInitialSoup
  initial_soup_rank(0),
  initial_soup_file(""),
  initial_soup_rotate(false),
  initial_soup_pressure_in(0.0),
  initial_soup_pressure_out(0.0),
  initial_soup_density(0.0),
  // EnzoInitialTurbulence
  initial_turbulence_density(0.0),
  initial_turbulence_pressure(0.0),
  initial_turbulence_temperature(0.0),
  // EnzoProlong
  interpolation_method(""),
  // EnzoMethodCheckGravity
  method_check_gravity_particle_type(),
  // EnzoMethodHeat
  method_heat_alpha(0.0),
  // EnzoMethodHydro
  method_hydro_method(""),
  method_hydro_dual_energy(false),
  method_hydro_dual_energy_eta_1(0.0),
  method_hydro_dual_energy_eta_2(0.0),
  method_hydro_reconstruct_method(""),
  method_hydro_reconstruct_conservative(0),
  method_hydro_reconstruct_positive(0),
  method_hydro_riemann_solver(""),
  // EnzoMethodTurbulence
  method_turbulence_edot(0.0),
  method_turbulence_mach_number(0.0),
  method_grackle_use_grackle(false),
#ifdef CONFIG_USE_GRACKLE
  method_grackle_chemistry(),
  method_grackle_use_cooling_timestep(false),
  method_grackle_radiation_redshift(-1.0),
#endif
  // EnzoMethodGravity
  method_gravity_grav_const(0.0),
  method_gravity_solver(""),
  method_gravity_order(4),
  method_gravity_accumulate(false),
  /// EnzoMethodPmDeposit
  method_pm_deposit_alpha(0.5),
  /// EnzoMethodPmUpdate
  method_pm_update_max_dt(std::numeric_limits<double>::max()),
  /// EnzoMethodMHDVlct
  method_vlct_riemann_solver(""),
  method_vlct_half_dt_reconstruct_method(""),
  method_vlct_full_dt_reconstruct_method(""),
  method_vlct_theta_limiter(0.0),
  method_vlct_density_floor(0.0),
  method_vlct_pressure_floor(0.0),
  method_vlct_mhd_choice(""),
  method_vlct_dual_energy(false),
  method_vlct_dual_energy_eta(0.0),
  /// EnzoProlong
  prolong_enzo_type(),
  prolong_enzo_positive(true),
  /// EnzoSolverMg0
  solver_pre_smooth(),
  solver_post_smooth(),
  solver_last_smooth(),
  solver_coarse_solve(),
  solver_domain_solve(),
  solver_weight(),
  solver_restart_cycle(),
  /// EnzoSolver<Krylov>
  solver_precondition(),
  solver_coarse_level(),
  solver_is_unigrid(),
  stopping_redshift()

{
  for (int i=0; i<3; i++) {
    initial_cloud_uniform_bfield[i] = 0;
    initial_sedov_array[i] = 0;
    initial_soup_array[i]  = 0;
    initial_soup_d_pos[i]  = 0.0;
    initial_soup_d_size[i] = 0.0;
    initial_collapse_array[i] = 0;
  }

#ifdef CONFIG_USE_GRACKLE
    method_grackle_chemistry = NULL;
#endif
}

//----------------------------------------------------------------------

EnzoConfig::~EnzoConfig() throw ()
{
#ifdef CONFIG_USE_GRACKLE
  if (method_grackle_chemistry){
    delete[] method_grackle_chemistry->grackle_data_file;
    delete method_grackle_chemistry;
  }
#endif // CONFIG_USE_GRACKLE
}

//----------------------------------------------------------------------

void EnzoConfig::pup (PUP::er &p)
{

  Config::pup(p);
  TRACEPUP;

  // NOTE: change this function whenever attributes change

  p | adapt_mass_type;

  p | ppm_diffusion;
  p | ppm_dual_energy;
  p | ppm_dual_energy_eta_1;
  p | ppm_dual_energy_eta_2;
  p | ppm_flattening;
  p | ppm_minimum_pressure_support_parameter;
  p | ppm_number_density_floor;
  p | ppm_density_floor;
  p | ppm_pressure_floor;
  p | ppm_pressure_free;
  p | ppm_temperature_floor;
  p | ppm_steepening;
  p | ppm_use_minimum_pressure_support;
  p | ppm_mol_weight;

  p | field_gamma;

  p | physics_cosmology;
  p | physics_cosmology_hubble_constant_now;
  p | physics_cosmology_omega_lamda_now;
  p | physics_cosmology_omega_matter_now;
  p | physics_cosmology_omega_baryon_now;
  p | physics_cosmology_omega_cdm_now;
  p | physics_cosmology_comoving_box_size;
  p | physics_cosmology_max_expansion_rate;
  p | physics_cosmology_initial_redshift;
  p | physics_cosmology_final_redshift;

  p | physics_gravity;

  p | initial_bcenter_update_etot;

  p | initial_cloud_subsample_n;
  p | initial_cloud_radius;
  p | initial_cloud_center_x;
  p | initial_cloud_center_y;
  p | initial_cloud_center_z;
  p | initial_cloud_density_cloud;
  p | initial_cloud_density_wind;
  p | initial_cloud_velocity_wind;
  p | initial_cloud_etot_wind;
  p | initial_cloud_eint_wind;
  p | initial_cloud_metal_mass_frac;
  p | initial_cloud_initialize_uniform_bfield;
  PUParray(p,initial_cloud_uniform_bfield,3);
  p | initial_cloud_perturb_stddev;
  p | initial_cloud_trunc_dev;
  p | initial_cloud_perturb_seed;

  p | initial_cosmology_temperature;

  p | initial_collapse_rank;
  PUParray(p,initial_collapse_array,3);
  p | initial_collapse_radius_relative;
  p | initial_collapse_particle_ratio;
  p | initial_collapse_mass;
  p | initial_collapse_temperature;

#ifdef CONFIG_USE_GRACKLE
  p | initial_grackle_test_minimum_H_number_density;
  p | initial_grackle_test_maximum_H_number_density;
  p | initial_grackle_test_minimum_temperature;
  p | initial_grackle_test_maximum_temperature;
  p | initial_grackle_test_minimum_metallicity;
  p | initial_grackle_test_maximum_metallicity;
  p | initial_grackle_test_reset_energies;
#endif /* CONFIG_USE_GRACKLE */

  p | initial_inclinedwave_alpha;
  p | initial_inclinedwave_beta;
  p | initial_inclinedwave_amplitude;
  p | initial_inclinedwave_lambda;
  p | initial_inclinedwave_parallel_vel;
  p | initial_inclinedwave_positive_vel;
  p | initial_inclinedwave_wave_type;

  p | initial_sedov_rank;
  PUParray(p,initial_sedov_array,3);
  p | initial_sedov_radius_relative;
  p | initial_sedov_pressure_in;
  p | initial_sedov_pressure_out;
  p | initial_sedov_density;

  PUParray(p,initial_sedov_random_array,3);
  p | initial_sedov_random_half_empty;
  p | initial_sedov_random_grackle_cooling;
  p | initial_sedov_random_max_blasts;
  p | initial_sedov_random_radius_relative;
  p | initial_sedov_random_pressure_in;
  p | initial_sedov_random_pressure_out;
  p | initial_sedov_random_density;
  p | initial_sedov_random_te_multiplier;

  p | initial_turbulence_density;
  p | initial_turbulence_pressure;
  p | initial_turbulence_temperature;

  p | initial_music_field_files;
  p | initial_music_field_datasets;
  p | initial_music_field_names;
  p | initial_music_field_coords;

  p | initial_music_particle_files;
  p | initial_music_particle_datasets;
  p | initial_music_particle_coords;
  p | initial_music_particle_types;
  p | initial_music_particle_attributes;
  p | initial_music_throttle_internode;
  p | initial_music_throttle_intranode;
  p | initial_music_throttle_node_files;
  p | initial_music_throttle_close_count;
  p | initial_music_throttle_group_size;
  p | initial_music_throttle_seconds_stagger;
  p | initial_music_throttle_seconds_delay;

  p | initial_pm_field;
  p | initial_pm_mpp;
  p | initial_pm_level;

  p | initial_shock_tube_setup_name;
  p | initial_shock_tube_aligned_ax;
  p | initial_shock_tube_axis_velocity;
  p | initial_shock_tube_trans_velocity;
  p | initial_shock_tube_flip_initialize;

  p | initial_soup_rank;
  p | initial_soup_file;
  p | initial_soup_rotate;
  PUParray(p,initial_soup_array,3);
  PUParray(p,initial_soup_d_pos,3);
  PUParray(p,initial_soup_d_size,3);
  p | initial_soup_pressure_in;
  p | initial_soup_pressure_out;
  p | initial_soup_density;

  p | interpolation_method;

  p | method_check_gravity_particle_type;

  p | method_heat_alpha;

  p | method_hydro_method;
  p | method_hydro_dual_energy;
  p | method_hydro_dual_energy_eta_1;
  p | method_hydro_dual_energy_eta_2;
  p | method_hydro_reconstruct_method;
  p | method_hydro_reconstruct_conservative;
  p | method_hydro_reconstruct_positive;
  p | method_hydro_riemann_solver;

  p | method_turbulence_edot;

  p | method_gravity_grav_const;
  p | method_gravity_solver;
  p | method_gravity_order;
  p | method_gravity_accumulate;

  p | method_pm_deposit_alpha;
  p | method_pm_update_max_dt;

  p | method_vlct_riemann_solver;
  p | method_vlct_half_dt_reconstruct_method;
  p | method_vlct_full_dt_reconstruct_method;
  p | method_vlct_theta_limiter;
  p | method_vlct_density_floor;
  p | method_vlct_pressure_floor;
  p | method_vlct_mhd_choice;
  p | method_vlct_dual_energy;
  p | method_vlct_dual_energy_eta;

  p | prolong_enzo_type;
  p | prolong_enzo_positive;

  p | solver_pre_smooth;
  p | solver_post_smooth;
  p | solver_last_smooth;
  p | solver_coarse_solve;
  p | solver_domain_solve;
  p | solver_weight;
  p | solver_restart_cycle;
  p | solver_precondition;
  p | solver_coarse_level;
  p | solver_is_unigrid;

  p | stopping_redshift;

  p | units_mass;
  p | units_density;
  p | units_length;
  p | units_time;

  p  | method_grackle_use_grackle;

#ifdef CONFIG_USE_GRACKLE
  if (method_grackle_use_grackle) {
    p  | method_grackle_use_cooling_timestep;
    p  | method_grackle_radiation_redshift;
    if (p.isUnpacking()) { method_grackle_chemistry = new chemistry_data; }
    p | *method_grackle_chemistry;
  } else {
    method_grackle_chemistry = nullptr;
  }
#endif /* CONFIG_USE_GRACKLE */

}

//----------------------------------------------------------------------

void EnzoConfig::read(Parameters * p) throw()
{
  TRACE("BEGIN EnzoConfig::read()");

  // Read Cello parameters


  TRACE("EnzoCharm::read calling Config::read()");

  ((Config*)this) -> read (p);

  adapt_mass_type.resize(num_adapt);

  for (int ia=0; ia<num_adapt; ia++) {

    std::string prefix = "Adapt:" + adapt_list[ia] + ":";
    adapt_mass_type[ia] = p->value_string(prefix+"mass_type","unknown");
    ASSERT2("EnzoConfig::read()",
	    "Unknown mass_type %s for parameter %s",
	    adapt_mass_type[ia].c_str(),(prefix+"mass_type").c_str(),
	    (adapt_type[ia] != "mass" ||
	     (adapt_mass_type[ia]=="dark" ||
	      adapt_mass_type[ia]=="baryon")));
  }

  double floor_default = 1e-6;

  ppm_diffusion = p->value_logical
    ("Method:ppm:diffusion", false);
  ppm_dual_energy = p->value_logical
    ("Method:ppm:dual_energy",false);
  ppm_dual_energy_eta_1 = p->value_float
    ("Method:ppm:dual_energy_eta_1", 0.001);
  ppm_dual_energy_eta_2 = p->value_float
    ("Method:ppm:dual_energy_eta_2", 0.1);
  ppm_flattening = p->value_integer
    ("Method:ppm:flattening", 3);
  ppm_minimum_pressure_support_parameter = p->value_integer
    ("Method:ppm:minimum_pressure_support_parameter",100);
  ppm_number_density_floor = p->value_float
    ("Method:ppm:number_density_floor", floor_default);
  ppm_density_floor = p->value_float
    ("Method:ppm:density_floor", floor_default);
  ppm_pressure_floor = p->value_float
    ("Method:ppm:pressure_floor", floor_default);
  ppm_pressure_free = p->value_logical
    ("Method:ppm:pressure_free",false);
  ppm_temperature_floor = p->value_float
    ("Method:ppm:temperature_floor", floor_default);
  ppm_steepening = p->value_logical
    ("Method:ppm:steepening", false);
  ppm_use_minimum_pressure_support = p->value_logical
    ("Method:ppm:use_minimum_pressure_support",false);
  ppm_mol_weight = p->value_float
    ("Method:ppm:mol_weight",0.6);

  // InitialMusic

  std::string name_initial = "Initial:music:";
  int num_files = p->list_length (name_initial + "file_list");
  for (int index_file=0; index_file<num_files; index_file++) {
    std::string file_id = name_initial +
      p->list_value_string (index_file,name_initial+"file_list") + ":";

    std::string type    = p->value_string (file_id+"type","");
    std::string name    = p->value_string (file_id+"name","");
    std::string file    = p->value_string (file_id+"file","");
    std::string dataset = p->value_string (file_id+"dataset","");
    std::string coords  = p->value_string (file_id+"coords","xyz");

    if (type == "particle") {
      std::string attribute = p->value_string (file_id+"attribute","");
      //      if (name != "") {
      initial_music_particle_files.     push_back(file);
      initial_music_particle_datasets.  push_back(dataset);
      initial_music_particle_coords.    push_back(coords);
      initial_music_particle_types.     push_back(name);
      initial_music_particle_attributes.push_back(attribute);
      //      }
    } else if (type == "field") {

      initial_music_field_files.        push_back(file);
      initial_music_field_datasets.     push_back(dataset);
      initial_music_field_names.        push_back(name);
      initial_music_field_coords.       push_back(coords);
    } else {
      ERROR2 ("EnzoConfig::read",
	      "Unknown particle type %s for parameter %s",
	      type.c_str(),(file_id+"type").c_str());
    }
  }
  // "sleep_by_process", "limit_per_node"
  initial_music_throttle_internode = p->value_logical
    ("Initial:music:throttle_internode",false);
  initial_music_throttle_intranode = p->value_logical
    ("Initial:music:throttle_intranode",false);
  initial_music_throttle_node_files = p->value_logical
    ("Initial:music:throttle_node_files",false);
  initial_music_throttle_close_count = p->value_integer
    ("Initial:music:throttle_close_count",0);
  initial_music_throttle_group_size = p->value_integer
    ("Initial:music:throttle_group_size",std::numeric_limits<int>::max());
  initial_music_throttle_seconds_stagger = p->value_float
    ("Initial:music:throttle_seconds_stagger",0.0);
  initial_music_throttle_seconds_delay = p->value_float
    ("Initial:music:throttle_seconds_delay",0.0);

  // PM method and initialization

  method_pm_deposit_alpha = p->value_float ("Method:pm_deposit:alpha",0.5);

  method_pm_update_max_dt = p->value_float
    ("Method:pm_update:max_dt", std::numeric_limits<double>::max());

  // ENZO interpolation
  prolong_enzo_type     = p->value_logical ("Prolong:enzo:type","2A");
  prolong_enzo_positive = p->value_logical ("Prolong:enzo:positive",true);

  // Particle method initialization
  initial_pm_field        = p->value_string  ("Initial:pm:field","density");
  initial_pm_mpp          = p->value_float   ("Initial:pm:mpp",-1.0);
  initial_pm_level        = p->value_integer ("Initial:pm:level",-1);

  field_gamma = p->value_float ("Field:gamma",5.0/3.0);

  // InitialInclinedWave initialization
  
  initial_inclinedwave_alpha          = p->value_float
    ("Initial:inclined_wave:alpha",0.0);
  initial_inclinedwave_beta           = p->value_float
    ("Initial:inclined_wave:beta",0.0);
  initial_inclinedwave_amplitude      = p->value_float
    ("Initial:inclined_wave:amplitude",1.e-6);
  initial_inclinedwave_lambda         = p->value_float
    ("Initial:inclined_wave:lambda",1.0);
  // The default vaue for parallel_vel is known by EnzoInitialInclinedWave
  // to mean that a value was not specified
  initial_inclinedwave_parallel_vel   = p->value_float
    ("Initial:inclined_wave:parallel_vel", std::numeric_limits<double>::min());
  initial_inclinedwave_positive_vel   = p->value_logical
    ("Initial:inclined_wave:positive_vel",true);
  initial_inclinedwave_wave_type      = p->value_string
    ("Initial:inclined_wave:wave_type","alfven");

  // InitialSoup initialization

  initial_soup_rank      = p->value_integer ("Initial:soup:rank",0);
  initial_soup_file      = p->value_string ("Initial:soup:file","soup.png");
  initial_soup_rotate    = p->value_logical ("Initial:soup:rotate",false);
  for (int axis=0; axis<3; axis++) {
    initial_soup_array[axis]  = p->list_value_integer
      (axis,"Initial:soup:array",1);
    initial_soup_d_pos[axis]  = p->list_value_float
      (axis,"Initial:soup:d_pos",0.0);
    initial_soup_d_size[axis] = p->list_value_float
      (axis,"Initial:soup:d_size",0.0);
  }
  initial_soup_pressure_in =
    p->value_float("Initial:soup:pressure_in",1.0);
  initial_soup_pressure_out =
    p->value_float("Initial:soup:pressure_out",1e-5);
  initial_soup_density =
    p->value_float("Initial:soup:density",1.0);

  // Sedov initialization

  TRACE1("field_gamma = %f",field_gamma);

  initial_sedov_rank = p->value_integer ("Initial:sedov:rank",0);

  initial_sedov_array[0] = p->list_value_integer (0,"Initial:sedov:array",1);
  initial_sedov_array[1] = p->list_value_integer (1,"Initial:sedov:array",1);
  initial_sedov_array[2] = p->list_value_integer (2,"Initial:sedov:array",1);

  initial_sedov_radius_relative =
    p->value_float("Initial:sedov:radius_relative",0.1);
  initial_sedov_pressure_in =
    p->value_float("Initial:sedov:pressure_in",1.0);
  initial_sedov_pressure_out =
    p->value_float("Initial:sedov:pressure_out",1e-5);
  initial_sedov_density =
    p->value_float("Initial:sedov:density",1.0);

  // Sedov Random Initialization

  initial_sedov_random_array[0] =
    p->list_value_integer (0,"Initial:sedov_random:array",1);
  initial_sedov_random_array[1] =
    p->list_value_integer (1,"Initial:sedov_random:array",1);
  initial_sedov_random_array[2] =
    p->list_value_integer (2,"Initial:sedov_random:array",1);

  initial_sedov_random_half_empty =
    p->value_logical ("Initial:sedov_random:half_empty",false);
  initial_sedov_random_grackle_cooling =
    p->value_logical ("Initial:sedov_random:grackle_cooling",false);
  initial_sedov_random_max_blasts =
    p->value_integer ("Initial:sedov_random:max_blasts",1);
  initial_sedov_random_radius_relative =
    p->value_float   ("Initial:sedov_random:radius_relative",0.1);
  initial_sedov_random_pressure_in =
    p->value_float   ("Initial:sedov_random:pressure_in",1.0);
  initial_sedov_random_pressure_out =
    p->value_float   ("Initial:sedov_random:pressure_out",1e-5);
  initial_sedov_random_density =
    p->value_float   ("Initial:sedov_random:density",1.0);
  initial_sedov_random_te_multiplier =
    p->value_integer  ("Initial:sedov_random:te_multiplier",1);

  // Shock Tube Initialization
  initial_shock_tube_setup_name = p->value_string
    ("Initial:shock_tube:setup_name","");
  initial_shock_tube_aligned_ax = p->value_string
    ("Initial:shock_tube:aligned_ax","x");
  initial_shock_tube_axis_velocity = p->value_float
    ("Initial:shock_tube:axis_velocity",0.0);
  initial_shock_tube_trans_velocity = p->value_float
    ("Initial:shock_tube:transverse_velocity",0.0);
  initial_shock_tube_flip_initialize = p -> value_logical
    ("Initial:shock_tube:flip_initialize", false);

  // VL+CT b-field initialization
  initial_bcenter_update_etot = p->value_logical
    ("Initial:vlct_bfield:update_etot",false);
  
  // Cloud Crush Initialization
  initial_cloud_subsample_n     = p->value_integer
    ("Initial:cloud:subsample_n",0);
  initial_cloud_radius          = p->value_float
    ("Initial:cloud:cloud_radius",0.0);
  initial_cloud_center_x        = p->value_float
    ("Initial:cloud:cloud_center_x",0.0);
  initial_cloud_center_y        = p->value_float
    ("Initial:cloud:cloud_center_y",0.0);
  initial_cloud_center_z        = p->value_float
    ("Initial:cloud:cloud_center_z",0.0);
  initial_cloud_density_cloud   = p->value_float
    ("Initial:cloud:cloud_density",0.0);
  initial_cloud_density_wind    = p->value_float
    ("Initial:cloud:wind_density",0.0);
  initial_cloud_velocity_wind   = p->value_float
    ("Initial:cloud:wind_velocity",0.0);
  initial_cloud_etot_wind       = p->value_float
    ("Initial:cloud:wind_total_energy",0.0);
  initial_cloud_eint_wind       = p->value_float
    ("Initial:cloud:wind_internal_energy",0.0);
  initial_cloud_metal_mass_frac = p->value_float
    ("Initial:cloud:metal_mass_fraction",0.0);
  initial_cloud_perturb_stddev  = p->value_float
    ("Initial:cloud:perturb_standard_deviation",0.0);
  initial_cloud_trunc_dev       = p->value_float
    ("Initial:cloud:perturb_truncation_deviation",0.0);
  int init_cloud_perturb_seed_  = p->value_integer
    ("Initial:cloud:perturb_seed",0);
  ASSERT("EnzoConfig::read()", "Initial:cloud:perturb_seed must be >=0",
	 init_cloud_perturb_seed_ >= 0);
  initial_cloud_perturb_seed = (unsigned int) init_cloud_perturb_seed_;

  int initial_cloud_uniform_bfield_length = p->list_length
    ("Initial:cloud:uniform_bfield");
  if (initial_cloud_uniform_bfield_length == 0){
    initial_cloud_initialize_uniform_bfield = false;
  } else if (initial_cloud_uniform_bfield_length == 3){
    initial_cloud_initialize_uniform_bfield = true;
    for (int i = 0; i <3; i++){
      initial_cloud_uniform_bfield[i] = p->list_value_float
	(i,"Initial:cloud:uniform_bfield");
    }
  } else {
    ERROR("EnzoConfig::read",
	  "Initial:cloud:uniform_bfield must contain 0 or 3 entries.");
  }

  // Cosmology initialization
  initial_cosmology_temperature = p->value_float("Initial:cosmology:temperature",0.0);

  // Collapse initialization

  initial_collapse_rank =  p->value_integer("Initial:collapse:rank",0);
  for (int i=0; i<initial_collapse_rank; i++) {
    initial_collapse_array[i] =
      p->list_value_integer (i,"Initial:collapse:array",1);
  }
  for (int i=initial_collapse_rank; i<3; i++) {
    initial_collapse_array[i] = 1;
  }
  initial_collapse_radius_relative =
    p->value_float("Initial:collapse:radius_relative",0.1);
  initial_collapse_particle_ratio =
    p->value_float("Initial:collapse:particle_ratio",0.0);
  initial_collapse_mass =
    p->value_float("Initial:collapse:mass",cello::mass_solar);
  initial_collapse_temperature =
    p->value_float("Initial:collapse:temperature",10.0);

  // Grackle test initialization
#ifdef CONFIG_USE_GRACKLE
  initial_grackle_test_minimum_H_number_density =
    p->value_float("Initial:grackle_test:minimum_H_number_density",0.1);
  initial_grackle_test_maximum_H_number_density =
    p->value_float("Initial:grackle_test:maximum_H_number_density",1000.0);
  initial_grackle_test_minimum_temperature =
    p->value_float("Initial:grackle_test:minimum_temperature",10.0);
  initial_grackle_test_maximum_temperature =
    p->value_float("Initial:grackle_test:maximum_temperature",1.0E8);
  initial_grackle_test_minimum_metallicity =
    p->value_float("Initial:grackle_test:minimum_metallicity", 1.0E-4);
  initial_grackle_test_maximum_metallicity =
    p->value_float("Initial:grackle_test:maximum_metallicity", 1.0);
  initial_grackle_test_reset_energies =
    p->value_integer("Initial:grackle_test:reset_energies",0);
#endif /* CONFIG_USE_GRACKLE */

  // Turbulence method and initialization

  initial_turbulence_density = p->value_float
    ("Initial:turbulence:density",1.0);

  // Must specify pressure or temperature
  initial_turbulence_pressure =    p->value_float
    ("Initial:turbulence:pressure",   0.0);
  initial_turbulence_temperature = p->value_float
    ("Initial:turbulence:temperature",0.0);

  bool uses_turbulence = false;
  for (size_t i=0; i<method_list.size(); i++) {
    if (method_list[i] == "turbulence") uses_turbulence=true;
  }

  if (uses_turbulence) {
    ASSERT ("EnzoConfig::read",
  	    "Either initial turbulence pressure or temperature must be defined",
  	    ! ((initial_turbulence_pressure == 0.0) &&
  	       (initial_turbulence_temperature == 0.0)));
    ASSERT ("EnzoConfig::read",
  	    "Initial turbulence pressure and temperature cannot "
	    "both be defined",
  	    ! ((initial_turbulence_pressure != 0.0) &&
  	       (initial_turbulence_temperature != 0.0)));
  }

  method_turbulence_edot = p->value_float
    ("Method:turbulence:edot",-1.0);
  method_turbulence_mach_number = p->value_float
    ("Method:turbulence:mach_number",0.0);

  interpolation_method = p->value_string
    ("Field:interpolation_method","SecondOrderA");

  method_check_gravity_particle_type = p->value_string
    ("Method:check_gravity:particle_type","dark");

  method_heat_alpha = p->value_float
    ("Method:heat:alpha",1.0);

  method_hydro_method = p->value_string
    ("Method:hydro:method","ppm");

  method_hydro_dual_energy = p->value_logical
    ("Method:hydro:dual_energy",false);
  method_hydro_dual_energy_eta_1 = p->value_float
    ("Method:hydro:dual_energy_eta_1",0.001);
  method_hydro_dual_energy_eta_2 = p->value_float
    ("Method:hydro:dual_energy_eta_2",0.1);

  method_hydro_reconstruct_method = p->value_string
    ("Method:hydro:reconstruct_method","ppm");

  method_hydro_reconstruct_conservative = p->value_logical
    ("Method:hydro:reconstruct_conservative",false);

  method_hydro_reconstruct_positive = p->value_logical
    ("Method:hydro:reconstruct_positive",false);

  method_hydro_riemann_solver = p->value_string
    ("Method:hydro:riemann_solver","ppm");

  method_gravity_grav_const = p->value_float
    ("Method:gravity:grav_const",6.67384e-8);

  method_gravity_solver = p->value_string
    ("Method:gravity:solver","unknown");

  method_gravity_order = p->value_integer
    ("Method:gravity:order",4);

  method_gravity_accumulate = p->value_logical
    ("Method:gravity:accumulate",true);

  method_vlct_riemann_solver = p->value_string
    ("Method:mhd_vlct:riemann_solver","hlld");
  method_vlct_half_dt_reconstruct_method = p->value_string
    ("Method:mhd_vlct:half_dt_reconstruct_method","nn");
  method_vlct_full_dt_reconstruct_method = p->value_string
    ("Method:mhd_vlct:full_dt_reconstruct_method","plm");
  method_vlct_theta_limiter = p->value_float
    ("Method:mhd_vlct:theta_limiter", 1.5);
  method_vlct_density_floor = p->value_float
    ("Method:mhd_vlct:density_floor", 0.0);
  method_vlct_pressure_floor = p->value_float
    ("Method:mhd_vlct:pressure_floor", 0.0);
  method_vlct_dual_energy = p->value_logical
    ("Method:mhd_vlct:dual_energy", false);
  method_vlct_dual_energy_eta = p->value_float
    ("Method:mhd_vlct:dual_energy_eta", 0.001);

  // we should raise an error if mhd_choice is not specified
  bool uses_vlct = false;
  for (size_t i=0; i<method_list.size(); i++) {
    if (method_list[i] == "mhd_vlct") uses_vlct=true;
  }
  method_vlct_mhd_choice = p->value_string
    ("Method:mhd_vlct:mhd_choice", "");
  if (uses_vlct && (method_vlct_mhd_choice == "")){
    ERROR("EnzoConfig::read", "Method:mhd_vlct:mhd_choice was not specified");
  }

  //--------------------------------------------------
  // Physics
  //--------------------------------------------------

  num_physics = p->list_length("Physics:list");

  for (int index_physics=0; index_physics<num_physics; index_physics++) {

    std::string name =
      p->list_value_string(index_physics,"Physics:list");

    std::string full_name = std::string("Physics:") + name;

    if (physics_list[index_physics] == "cosmology") {

      physics_cosmology = true;

      physics_cosmology_hubble_constant_now = p->value_float
	(full_name + ":hubble_constant_now",0.701);

      physics_cosmology_omega_matter_now = p->value_float
	(full_name + ":omega_matter_now",   0.279);

      physics_cosmology_omega_baryon_now = p->value_float
	(full_name + ":omega_baryon_now",   1.0);

      physics_cosmology_omega_cdm_now = p->value_float
	(full_name + ":omega_cdm_now",   0.0);

      physics_cosmology_omega_lamda_now = p->value_float
	(full_name + ":omega_lambda_now",   0.721);


      physics_cosmology_comoving_box_size = p->value_float
	(full_name + ":comoving_box_size", 64.0);

      physics_cosmology_max_expansion_rate = p->value_float
	(full_name + ":max_expansion_rate", 0.01);

      physics_cosmology_initial_redshift = p->value_float
	(full_name + ":initial_redshift",  20.0);;

      physics_cosmology_final_redshift = p->value_float
	(full_name + ":final_redshift",  0.0);;

    }

    if (physics_list[index_physics] == "gravity") {

      physics_gravity = true;

    }
  }

  //======================================================================
  // SOLVER
  //======================================================================

  num_solvers = p->list_length("Solver:list");

  solver_pre_smooth.  resize(num_solvers);
  solver_coarse_solve.resize(num_solvers);
  solver_domain_solve.resize(num_solvers);
  solver_post_smooth. resize(num_solvers);
  solver_last_smooth. resize(num_solvers);
  solver_weight.      resize(num_solvers);
  solver_restart_cycle.resize(num_solvers);
  solver_precondition.resize(num_solvers);
  solver_coarse_level.resize(num_solvers);
  solver_is_unigrid.resize(num_solvers);

  for (int index_solver=0; index_solver<num_solvers; index_solver++) {

    std::string solver_name =
      std::string("Solver:") + p->list_value_string(index_solver,"Solver:list");

    std::string solver;

    solver = p->value_string (solver_name + ":precondition","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_precondition[index_solver] = solver_index[solver];
    } else {
      solver_precondition[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":pre_smooth","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_pre_smooth[index_solver] = solver_index[solver];
    } else {
      solver_pre_smooth[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":coarse_solve","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_coarse_solve[index_solver] = solver_index[solver];
    } else {
      solver_coarse_solve[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":domain_solve","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_domain_solve[index_solver] = solver_index[solver];
    } else {
      solver_domain_solve[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":post_smooth","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_post_smooth[index_solver] = solver_index[solver];
    } else {
      solver_post_smooth[index_solver] = -1;
    }

    solver = p->value_string (solver_name + ":last_smooth","unknown");
    if (solver_index.find(solver) != solver_index.end()) {
      solver_last_smooth[index_solver] = solver_index[solver];
    } else {
      solver_last_smooth[index_solver] = -1;
    }

    solver_weight[index_solver] =
      p->value_float(solver_name + ":weight",1.0);

    solver_restart_cycle[index_solver] =
      p->value_integer(solver_name + ":restart_cycle",1);

    solver_coarse_level[index_solver] = 
      p->value_integer (solver_name + ":coarse_level",
			solver_min_level[index_solver]);

    solver_is_unigrid[index_solver] =
      p->value_logical (solver_name + ":is_unigrid",false);

  }

  //======================================================================
  // STOPPING
  //======================================================================

  stopping_redshift = p->value_float ("Stopping:redshift",0.0);

  //======================================================================
  // GRACKLE 3.0
  //======================================================================

  this->method_grackle_use_grackle = false;
#ifdef CONFIG_USE_GRACKLE


  /// Grackle parameters

  for (size_t i=0; i<method_list.size(); i++) {
    if (method_list[i] == "grackle") method_grackle_use_grackle=true;
  }

  // Defaults alert PUP::er() to ignore
  if (method_grackle_use_grackle) {

    method_grackle_chemistry = new chemistry_data;
    *method_grackle_chemistry = _set_default_chemistry_parameters();

    /* this must be set AFTER default values are set */
    method_grackle_chemistry->use_grackle = method_grackle_use_grackle;

    // Copy over parameters from Enzo-P to Grackle
    method_grackle_chemistry->Gamma = field_gamma;

    //
    method_grackle_use_cooling_timestep = p->value_logical
      ("Method:grackle:use_cooling_timestep", false);

    // for when not using cosmology - redshift of UVB
    method_grackle_radiation_redshift = p->value_float
      ("Method:grackle:radiation_redshift", -1.0);

    // Set Grackle parameters from parameter file
    method_grackle_chemistry->with_radiative_cooling = p->value_integer
      ("Method:grackle:with_radiative_cooling",
        method_grackle_chemistry->with_radiative_cooling);

    method_grackle_chemistry->primordial_chemistry = p->value_integer
      ("Method:grackle:primordial_chemistry",
        method_grackle_chemistry->primordial_chemistry);

    method_grackle_chemistry->metal_cooling = p->value_integer
      ("Method:grackle:metal_cooling",
        method_grackle_chemistry->metal_cooling);

    method_grackle_chemistry->h2_on_dust = p->value_integer
      ("Method:grackle:h2_on_dust",
        method_grackle_chemistry->h2_on_dust);

    method_grackle_chemistry->three_body_rate = p->value_integer
      ("Method:grackle:three_body_rate",
        method_grackle_chemistry->three_body_rate);

    method_grackle_chemistry->cmb_temperature_floor = p->value_integer
      ("Method:grackle:cmb_temperature_floor",
        method_grackle_chemistry->cmb_temperature_floor);

    std::string grackle_data_file_ = p->value_string
      ("Method:grackle:data_file", "");
    ASSERT("EnzoConfig::read",
	   "no value specified for \"Method:grackle:data_file\"",
	   grackle_data_file_.length() > 0);

    method_grackle_chemistry->grackle_data_file
      = new char[grackle_data_file_.length() + 1];
    strcpy(method_grackle_chemistry->grackle_data_file,
	   grackle_data_file_.c_str());

    method_grackle_chemistry->cie_cooling = p->value_integer
      ("Method:grackle:cie_cooling",
        method_grackle_chemistry->cie_cooling);

    method_grackle_chemistry->h2_optical_depth_approximation = p->value_integer
      ("Method:grackle:h2_optical_depth_approximation",
        method_grackle_chemistry->h2_optical_depth_approximation);

    method_grackle_chemistry->photoelectric_heating = p->value_integer
      ("Method:grackle:photoelectric_heating",
        method_grackle_chemistry->photoelectric_heating);

    method_grackle_chemistry->photoelectric_heating_rate = p->value_float
      ("Method:grackle:photoelectric_heating_rate",
        method_grackle_chemistry->photoelectric_heating_rate);

    method_grackle_chemistry->CaseBRecombination = p->value_integer
      ("Method:grackle:CaseBRecombination",
       method_grackle_chemistry->CaseBRecombination);

    method_grackle_chemistry->UVbackground = p->value_integer
      ("Method:grackle:UVbackground",
        method_grackle_chemistry->UVbackground);

    method_grackle_chemistry->use_volumetric_heating_rate = p->value_integer
      ("Method:grackle:use_volumetric_heating_rate",
      method_grackle_chemistry->use_volumetric_heating_rate);

    method_grackle_chemistry->use_specific_heating_rate = p->value_integer
      ("Method:grackle:use_specific_heating_rate",
       method_grackle_chemistry->use_specific_heating_rate);

    method_grackle_chemistry->self_shielding_method = p->value_integer
      ("Method:grackle:self_shielding_method",
       method_grackle_chemistry->self_shielding_method);

    method_grackle_chemistry->H2_self_shielding = p->value_integer
      ("Method:grackle:H2_self_shielding",
       method_grackle_chemistry->H2_self_shielding);

    method_grackle_chemistry->HydrogenFractionByMass = p->value_float
      ("Method:grackle:HydrogenFractionByMass",
        method_grackle_chemistry->HydrogenFractionByMass);

    method_grackle_chemistry->DeuteriumToHydrogenRatio = p->value_float
      ("Method:grackle:DeuteriumToHydrogenRatio",
       method_grackle_chemistry->DeuteriumToHydrogenRatio);

    method_grackle_chemistry->SolarMetalFractionByMass = p->value_float
      ("Method:grackle:SolarMetalFractionByMass",
       method_grackle_chemistry->SolarMetalFractionByMass);

    method_grackle_chemistry->Compton_xray_heating = p->value_integer
      ("Method:grackle:Compton_xray_heating",
       method_grackle_chemistry->Compton_xray_heating);

    method_grackle_chemistry->LWbackground_sawtooth_suppression = p->value_integer
      ("Method:grackle:LWbackground_sawtooth_suppression",
       method_grackle_chemistry->LWbackground_sawtooth_suppression);

    method_grackle_chemistry->LWbackground_intensity = p->value_float
      ("Method:grackle:LWbackground_intensity",
       method_grackle_chemistry->LWbackground_intensity);

    method_grackle_chemistry->UVbackground_redshift_on = p->value_float
      ("Method:grackle:UVbackground_redshift_on",
       method_grackle_chemistry->UVbackground_redshift_on);

    method_grackle_chemistry->UVbackground_redshift_off = p->value_float
      ("Method:grackle:UVbackground_redshift_off",
       method_grackle_chemistry->UVbackground_redshift_off);

    method_grackle_chemistry->UVbackground_redshift_fullon = p->value_float
      ("Method:grackle:UVbackground_redshift_fullon",
        method_grackle_chemistry->UVbackground_redshift_fullon);

    method_grackle_chemistry->UVbackground_redshift_drop = p->value_float
     ("Method:grackle:UVbackground_redshift_drop",
      method_grackle_chemistry->UVbackground_redshift_drop);

    // When radiative transfer is eventually included, make
    // sure to set the below parameter to match the Enzo-P
    // parameter for turning RT on / off:
    //   method_grackle_chemistry->use_radiative_transfer = ENZO_P_PARAMETER_NAME;

  }
#endif /* CONFIG_USE_GRACKLE */

  TRACE("END   EnzoConfig::read()");
}

//======================================================================
