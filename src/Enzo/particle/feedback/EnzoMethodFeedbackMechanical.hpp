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

protected: // methods
   void allocate_temporary_(EnzoBlock * enzo_block)
   {
     Field field = enzo_block->data()->field();
     field.allocate_temporary(i_nsn);
     field.allocate_temporary(i_dep_mass);
     field.allocate_temporary(i_dep_metl);
     field.allocate_temporary(i_dep_snii);
     field.allocate_temporary(i_dep_snia);
     field.allocate_temporary(i_dep_px);
     field.allocate_temporary(i_dep_py);
     field.allocate_temporary(i_dep_pz);

     field.allocate_temporary(i_nsn_a);
     field.allocate_temporary(i_dep_mass_a);
     field.allocate_temporary(i_dep_metl_a);
     field.allocate_temporary(i_dep_snii_a);
     field.allocate_temporary(i_dep_snia_a);
     field.allocate_temporary(i_dep_px_a);
     field.allocate_temporary(i_dep_py_a);
     field.allocate_temporary(i_dep_pz_a);
   }

   void deallocate_temporary_(EnzoBlock * enzo_block)
   {
     Field field = enzo_block->data()->field();
     field.deallocate_temporary(i_nsn);
     field.deallocate_temporary(i_dep_mass);
     field.deallocate_temporary(i_dep_metl);
     field.deallocate_temporary(i_dep_snii);
     field.deallocate_temporary(i_dep_snia);
     field.deallocate_temporary(i_dep_px);
     field.deallocate_temporary(i_dep_py);
     field.deallocate_temporary(i_dep_pz);

     field.deallocate_temporary(i_nsn_a);
     field.deallocate_temporary(i_dep_mass_a);
     field.deallocate_temporary(i_dep_metl_a);
     field.deallocate_temporary(i_dep_snii_a);
     field.deallocate_temporary(i_dep_snia_a);
     field.deallocate_temporary(i_dep_px_a);
     field.deallocate_temporary(i_dep_py_a);
     field.deallocate_temporary(i_dep_pz_a);
   }

protected:

  // configuration parameters (these are directly set by the user)
  bool stochastic_;
  double ejecta_mass_fraction_;
  double ejecta_metal_fraction_;
  int min_nsn_per_timestep_;

  // Refresh ID
  int ir_feedback_;

  // deposit field id's
  // _a fields are for accumulation
  int i_nsn, i_nsn_a;
  int i_dep_mass, i_dep_mass_a;
  int i_dep_metl, i_dep_metl_a;
  int i_dep_snii, i_dep_snii_a; 
  int i_dep_snia, i_dep_snia_a;
  int i_dep_px, i_dep_px_a;
  int i_dep_py, i_dep_py_a;
  int i_dep_pz, i_dep_pz_a;

};

#endif
