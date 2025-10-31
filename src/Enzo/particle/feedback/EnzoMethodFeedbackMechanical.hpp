/// See LICENSE_CELLO file for license and copyright information

/// @file	  enzo_EnzoMethodFeedbackMechanical.hpp
/// @author Claire Kopenhafer (kopenhaf@msu.edu, clairekope@gmail.com)
//             
/// @date
/// @brief  Implements the mechanical feedback method in Enzo that was
///         written by Cassi Lochhaas for the FOGGIE collaboration.
///         This method was based on Kimm & Cen 2014. Since this method
///         is similar to the STARSS algorithm, this Enzo-E adaptation
///         is built using the STARSS implementation as a starting point.


#ifndef ENZO_ENZO_METHOD_FEEDBACK_MECHANICAL
#define ENZO_ENZO_METHOD_FEEDBACK_MECHANICAL

class EnzoMethodFeedbackMechanical : public Method {

  /// @class   EnzoMethodFeedbackMechanical 
  /// @ingroup Enzo
  /// @btief   [\ref Enzo] Encapsulate Feedback Routines

public:

  EnzoMethodFeedbackMechanical(ParameterGroup p);

  /// Destructor
  virtual ~EnzoMethodFeedbackMechanical() throw() {};

  /// Charm++ Pup::able declarations
  PUPable_decl(EnzoMethodFeedbackMechanical);

  /// Charm++ Pup::able migration Constructor
  EnzoMethodFeedbackMechanical (CkMigrateMessage *m)
   : Method (m)
   , ir_feedback_(-1)  // ? used for dummy fields
  {  }

  /// Charm++ Pack / Unpack function
  void pup(PUP::er &p);

  /// Apply the method
  virtual void compute (Block * block) throw();
  
  void compute_ (Block * block);  // no throw()?
  
  /// name
  virtual std::string name() throw()
  { return "feedback"; }

  // Compute the maximum timestep for this method
  virtual double timestep (Block * block) throw();
  
  void add_accumulate_fields(EnzoBlock * enzo_block) throw();
   
protected:

  void allocate_temporary_yields_(EnzoBlock * enzo_block)
  {
    Field field = enzo_block->data()->field();
    field.allocate_temporary(i_yld_nsn);
    field.allocate_temporary(i_yld_mass);
    field.allocate_temporary(i_yld_metl);
    field.allocate_temporary(i_yld_vxp);
    field.allocate_temporary(i_yld_vyp);
    field.allocate_temporary(i_yld_vzp);

    if (track_metal_sources_) {
      field.allocate_temporary(i_yld_snii);
      field.allocate_temporary(i_yld_snia);
    }
  }

  void deallocate_temporary_yields_(EnzoBlock * enzo_block)
  {
    Field field = enzo_block->data()->field();
    field.deallocate_temporary(i_yld_nsn);
    field.deallocate_temporary(i_yld_mass);
    field.deallocate_temporary(i_yld_metl);
    field.deallocate_temporary(i_yld_vxp);
    field.deallocate_temporary(i_yld_vyp);
    field.deallocate_temporary(i_yld_vzp);
    
    if (track_metal_sources_) {
      field.deallocate_temporary(i_yld_snii);
      field.deallocate_temporary(i_yld_snia);
    }
  }

  void allocate_temporary_fluids_(EnzoBlock * enzo_block)
  {
    Field field = enzo_block->data()->field();

    field.allocate_temporary(i_mu);

    field.allocate_temporary(i_d_dep);
    field.allocate_temporary(i_te_dep);
    field.allocate_temporary(i_ge_dep);
    field.allocate_temporary(i_md_dep);
    field.allocate_temporary(i_vx_dep);
    field.allocate_temporary(i_vy_dep);
    field.allocate_temporary(i_vz_dep);

    field.allocate_temporary(i_d_dep_a);
    field.allocate_temporary(i_te_dep_a);
    field.allocate_temporary(i_ge_dep_a);
    field.allocate_temporary(i_md_dep_a);
    field.allocate_temporary(i_vx_dep_a);
    field.allocate_temporary(i_vy_dep_a);
    field.allocate_temporary(i_vz_dep_a);

    if (track_metal_sources_) {
      field.allocate_temporary(i_mdii_dep);
      field.allocate_temporary(i_mdia_dep);
      field.allocate_temporary(i_mdii_dep_a);
      field.allocate_temporary(i_mdia_dep_a);
    }
  }

  void deallocate_temporary_fluids_(EnzoBlock * enzo_block)
  {
    Field field = enzo_block->data()->field();

    field.deallocate_temporary(i_mu);

    field.deallocate_temporary(i_d_dep);
    field.deallocate_temporary(i_te_dep);
    field.deallocate_temporary(i_ge_dep);
    field.deallocate_temporary(i_md_dep);
    field.deallocate_temporary(i_vx_dep);
    field.deallocate_temporary(i_vy_dep);
    field.deallocate_temporary(i_vz_dep);

    field.deallocate_temporary(i_d_dep_a);
    field.deallocate_temporary(i_te_dep_a);
    field.deallocate_temporary(i_ge_dep_a);
    field.deallocate_temporary(i_md_dep_a);
    field.deallocate_temporary(i_vx_dep_a);
    field.deallocate_temporary(i_vy_dep_a);
    field.deallocate_temporary(i_vz_dep_a);

    if (track_metal_sources_) {
      field.deallocate_temporary(i_mdii_dep);
      field.deallocate_temporary(i_mdia_dep);
      field.deallocate_temporary(i_mdii_dep_a);
      field.deallocate_temporary(i_mdia_dep_a);
    }
  }

  void compute_molecular_weight_(EnzoBlock * enzo_block);

  // configuration parameters (these are directly set by the user)
  bool stochastic_;
  bool pre_sne_;
  double ejecta_mass_fraction_;
  double ejecta_metal_fraction_;
  int min_nsn_per_timestep_;
  double momentum_mult_;
  bool cap_velocity_kick_;
  bool track_metal_sources_;

  // Refresh ID
  int ir_feedback_;

  // a temporary field for mean molecular weight
  int i_mu;

  // deposit field id's
  // these are used for figuring out feedback
  int i_yld_nsn;
  int i_yld_mass;
  int i_yld_metl;
  int i_yld_snii;
  int i_yld_snia;
  int i_yld_vxp;
  int i_yld_vyp;
  int i_yld_vzp;
  
  // these are used for accumulating across ghost zones
  // _a fields are for accumulation
  int i_d_dep , i_d_dep_a;
  int i_te_dep, i_te_dep_a;
  int i_ge_dep, i_ge_dep_a;
  int i_md_dep, i_md_dep_a;
  int i_mdii_dep, i_mdii_dep_a;
  int i_mdia_dep, i_mdia_dep_a;
  int i_vx_dep, i_vx_dep_a;
  int i_vy_dep, i_vy_dep_a;
  int i_vz_dep, i_vz_dep_a;

};

#endif
