////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Tianhong Yu.
//
// This is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This software is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
//! \file
//! \brief Example code for running a network calculation for a single type
//!        of particle.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>

#include "nnt/write_output_xml.h"
#include "user/remove_duplicate.h"
#include "user/my_network_limiter.h"
#include "user/my_rate_modifiers.h"
#include "user/my_user_rate_functions.h"

#include "carbon_evolve.h"
#include "carbon_rate_functions.h"
#include "carbon_hydro.h"
#include "my_bin_utilities.h"
#include "carbon_molecule_utilities.h"
#include "carbon_reaction_utilities.h"

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_DT0          1.e-20       // Initial time step
#define D_REG_T        0.15         // Time step change regulator for dt update
#define D_REG_Y        0.15         // Abundance change regulator for dt update 
#define D_Y_MIN_DT     1.e-10       // Smallest y for dt update
#define S_SOLVER       nnt::s_ARROW // Solver type: ARROW or GSL

#define S_NETWORK_END  "network end"
#define S_PHOTON_END   "photon end"

int
create_generic_net( nnt::Zone & );

void
my_print_abundances( nnt::Zone & );

int
print_abundance( Libnucnet__Species *, Libnucnet__Zone * );

int
my_species_sort_function(
  const Libnucnet__Species *, const Libnucnet__Species *
);

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step, i_label = 0;
  double d_t, d_dt;
  Libnucnet *p_my_nucnet = NULL, *p_my_output;
  nnt::Zone zone;

  //============================================================================
  // Get the nucnet.
  //============================================================================

  p_my_nucnet = get_nucnet( argc, argv );

  //============================================================================
  // Get the zone.
  //============================================================================

  if( !set_zone( p_my_nucnet, zone, argv ) )
  {
    std::cerr << "Couldn't set zone."  << std::endl;
    return EXIT_FAILURE;
  }

  //============================================================================
  // Create generic network.
  //============================================================================

  create_generic_net( zone );

  //============================================================================
  // Register rate functions.
  //============================================================================

  register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Update user rate data. 
  //============================================================================

  update_my_rate_functions_data( zone );

  //============================================================================
  // Initialize time.
  //============================================================================

  if( zone.hasProperty( nnt::s_DTIME ) )
    d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );
  else
  {
    d_dt = D_DT0;
    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );
  }

  if( zone.hasProperty( nnt::s_TIME ) )
    d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );
  else
  {
    d_t = 0;
    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );
  }

  //============================================================================
  // Initialize zone.
  //============================================================================

  initialize_zone( zone, argv );

  if( 
    zone.hasProperty( "run bin" ) && zone.getProperty( "run bin" ) == "yes" 
  )
    initialize_bin( zone );

  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if( S_SOLVER == nnt::s_ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function) my_species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    zone.updateProperty( nnt::s_SOLVER, nnt::s_ARROW );

    zone.updateProperty( nnt::s_ARROW_WIDTH, "1" );

  }

  //============================================================================
  // Create output.
  //============================================================================

  p_my_output = nnt::create_output( p_my_nucnet );

  Libnucnet__setZoneCompareFunction(
    p_my_output,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  //============================================================================
  // Toggle off the reverse rates for reactions from detailed balance. 
  //============================================================================

  Libnucnet__Zone__toggleReverseRateDetailedBalance(
    zone.getNucnetZone(), "off"
  );

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  i_step = 0;

    std::cout << "species number = " <<
      Libnucnet__Nuc__getNumberOfSpecies(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        )
      ) +
      boost::lexical_cast<int>( zone.getProperty( "bin number" ) )
      << std::endl;

  while ( 
    d_t < boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) )
  )
  {

    d_t += d_dt;

    //==========================================================================
    // Set dt and t.
    //==========================================================================

    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

    //==========================================================================
    // Update temperature and density.
    //==========================================================================

    update_zone_properties( zone );

    if( 
      zone.hasProperty( "run bin" ) && zone.getProperty( "run bin" ) == "yes" 
    )
      update_bin_properties( zone );

    double d_rho = 
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

    if( d_rho < 0. ) {
      std::cout << "rho is " << d_rho << std::endl;
      exit(EXIT_FAILURE);
    }

    //==========================================================================
    // Evolve abundances.
    //==========================================================================

    evolve( zone );

    //==========================================================================
    // Print out abundances.
    //==========================================================================

    if(
      i_step % 
        boost::lexical_cast<int>( zone.getProperty( nnt::s_STEPS ) ) == 0 ||
      d_t >= boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) )
    )
    {

      Libnucnet__relabelZone(
        p_my_nucnet,
        zone.getNucnetZone(),
        ( boost::lexical_cast<std::string>( ++i_label ) ).c_str(),
        NULL,
        NULL
      );

      if( 
        zone.hasProperty( "idl print" ) &&
        zone.getProperty( "idl print" ) == "yes"
      ) {
        my_print_abundances( zone );
      } else {
        zone.printAbundances();
      }

      zone.updateProperty(
        nnt::s_YE,
        boost::lexical_cast<std::string>(
          Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
        )
      );

      nnt::write_xml( p_my_output, zone.getNucnetZone() );

    }

    //==========================================================================
    // Update timestep.
    //==========================================================================

    d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if( 
      zone.hasProperty( "run bin" ) && zone.getProperty( "run bin" ) == "yes" 
    )
      bin_update_timestep( zone, d_dt, D_REG_T, D_REG_Y, D_Y_MIN_DT );

    if( boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) > 10. )
      zone.normalizeAbundances();

    if(
      d_t + d_dt >
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) )
    )
    {
      d_dt =
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) ) - d_t;
    }

    i_step++;

  }  

  //============================================================================
  // Write output.
  //============================================================================

  Libnucnet__updateZoneXmlMassFractionFormat(
    p_my_output,
    "%.15e"
  );

  Libnucnet__writeToXmlFile(
    p_my_output,
    zone.getProperty( S_OUTPUT_XML_FILE ).c_str()
  );

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );
  Libnucnet__free( p_my_output );

  return EXIT_SUCCESS;

}

//############################################################################
// my_print_abundances()
//############################################################################

void
my_print_abundances(
  nnt::Zone & zone
)
{

  double d_t, d_dt, d_t9, d_rho;

  d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

  d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  d_t9 = boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) );

  d_rho = boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

  printf(
    "t = %10.4e\ndt = %10.4e\nt9 = %10.4e\nrho (g/cc) = %10.4e\n\n",
    d_t, d_dt, d_t9, d_rho
  );

  Libnucnet__Nuc__setSpeciesCompareFunction(
    Libnucnet__Net__getNuc( 
      Libnucnet__Zone__getNet( zone.getNucnetZone() ) 
    ),
    (Libnucnet__Species__compare_function) nnt::species_sort_by_z_then_a
  );

  Libnucnet__Nuc__sortSpecies(
    Libnucnet__Net__getNuc( 
      Libnucnet__Zone__getNet( zone.getNucnetZone() ) 
    ) 
  );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc(
      Libnucnet__Zone__getNet( zone.getNucnetZone() )
    ),
    (Libnucnet__Species__iterateFunction) print_abundance,
    zone.getNucnetZone()
  );

  if( 
    zone.hasProperty( "run bin" ) && zone.getProperty( "run bin" ) == "yes" 
  ) {

    print_bin_abundances( zone );

    std::cout << "1 - Xsum = " <<
      1. - 
      Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) -
      compute_bin_sum( zone ) <<
      std::endl << std::endl;

  } else {

    std::cout << "1 - Xsum = " <<
      1. - 
      Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) <<
      std::endl << std::endl;

  }

  Libnucnet__Nuc__setSpeciesCompareFunction(
    Libnucnet__Net__getNuc( 
      Libnucnet__Zone__getNet( zone.getNucnetZone() ) 
    ),
    (Libnucnet__Species__compare_function) my_species_sort_function
  );

  Libnucnet__Nuc__sortSpecies(
    Libnucnet__Net__getNuc( 
      Libnucnet__Zone__getNet( zone.getNucnetZone() ) 
    ) 
  );

}

//############################################################################
// print_abundance()
//############################################################################

int
print_abundance(
  Libnucnet__Species *p_species,
  Libnucnet__Zone *p_zone
)
{

  double d_abund;

  d_abund =
    Libnucnet__Zone__getSpeciesAbundance( p_zone, p_species );

  printf( "%20lu%20lu%16.6e%16.6e\n",
    (unsigned long) Libnucnet__Species__getZ( p_species ),
    (unsigned long) Libnucnet__Species__getA( p_species ),
    d_abund,
    d_abund *
      boost::lexical_cast<double>( Libnucnet__Species__getA( p_species ) )
  );

  return 1;

}

//############################################################################
// my_species_sort_function()
//############################################################################

int
my_species_sort_function(
  const Libnucnet__Species *p_species1,
  const Libnucnet__Species *p_species2
)
{

  int i;

  if( !strcmp( Libnucnet__Species__getName( p_species1 ), "h1" ) )
    return 1;
  else if( !strcmp( Libnucnet__Species__getName( p_species2 ), "h1" ) )
    return -1;

  if( !strcmp( Libnucnet__Species__getName( p_species1 ), "h1g" ) )
    return 1;
  else if( !strcmp( Libnucnet__Species__getName( p_species2 ), "h1g" ) )
    return -1;

  if( 
      Libnucnet__Species__getZ( p_species1 ) <
      Libnucnet__Species__getZ( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getZ( p_species1 ) >
      Libnucnet__Species__getZ( p_species2 )
  )
    return 1;

  if( 
      Libnucnet__Species__getA( p_species1 ) <
      Libnucnet__Species__getA( p_species2 )
  )
    return -1;
  else if(
      Libnucnet__Species__getA( p_species1 ) >
      Libnucnet__Species__getA( p_species2 )
  )
    return 1;

  i =
    strcmp(
      Libnucnet__Species__getName( p_species1 ),
      Libnucnet__Species__getName( p_species2 )
    );

  if( i == 0 ) {
    return 0;
  } else {
    return GSL_SIGN( i );
  }

}

//##############################################################################
// create_generic_net(). 
//##############################################################################

int
create_generic_net(
  nnt::Zone & zone
)
{

  unsigned i_network_end;
  unsigned i_photon_end;

  i_network_end = 
    boost::lexical_cast<unsigned>( zone.getProperty( S_NETWORK_END ) ); 
  i_photon_end = 
    boost::lexical_cast<unsigned>( zone.getProperty( S_PHOTON_END ) ); 

  add_generic_molecules_to_nuc(
    Libnucnet__Net__getNuc( 
      Libnucnet__Zone__getNet( zone.getNucnetZone() )
    ),
    1, 
    i_network_end
  );

  add_generic_reactions_to_net(
    Libnucnet__Zone__getNet( zone.getNucnetZone() ),
    1,
    i_photon_end
  ); 

  if( i_network_end > i_photon_end )
    add_forward_generic_reactions_to_net(
      Libnucnet__Zone__getNet( zone.getNucnetZone() ),
      i_photon_end,
      i_network_end
    ); 

  return 1;

}

