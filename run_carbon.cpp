////////////////////////////////////////////////////////////////////////////////
// This file was originally written by Bradley S. Meyer.
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
//! \brief Example code for running a single zone network calculation.
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
#include "user/evolve.h"
#include "carbon_rate_functions.h"
#include "carbon_hydro.h"

//##############################################################################
// Define some parameters.
//##############################################################################

#define D_DT0          1.e-20       // Initial time step
#define D_REG_T        0.15         // Time step change regulator for dt update
#define D_REG_Y        0.15         // Abundance change regulator for dt update 
#define D_Y_MIN_DT     1.e-10       // Smallest y for dt update
#define S_SOLVER       nnt::s_ARROW // Solver type: ARROW or GSL

//##############################################################################
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step, k = 0;
  double d_t, d_dt;
  Libnucnet *p_my_nucnet = NULL, *p_my_output;
  nnt::Zone zone;

  //============================================================================
  // Get the nucnet.
  //============================================================================

  p_my_nucnet = get_nucnet( argc, argv );

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  register_my_rate_functions(
    Libnucnet__Net__getReac( Libnucnet__getNet( p_my_nucnet ) )
  );

  //============================================================================
  // Get the zone.
  //============================================================================

  if( !set_zone( p_my_nucnet, zone, argv ) )
  {
    std::cerr << "Couldn't set zone."  << std::endl;
    return EXIT_FAILURE;
  }

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

  //============================================================================
  // Sort the nuclei if using the arrow solver.
  //============================================================================

  if( S_SOLVER == nnt::s_ARROW )
  {

    Libnucnet__Nuc__setSpeciesCompareFunction(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
      (Libnucnet__Species__compare_function) nnt::species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    zone.updateProperty( nnt::s_SOLVER, nnt::s_ARROW );

    zone.updateProperty( nnt::s_ARROW_WIDTH, "2" );

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

  while ( d_t < boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) ) )
  {

    d_t += d_dt;

  //============================================================================
  // Set dt and t.
  //============================================================================

    zone.updateProperty(
      nnt::s_DTIME,
      boost::lexical_cast<std::string>( d_dt )
    );

    zone.updateProperty(
      nnt::s_TIME,
      boost::lexical_cast<std::string>( d_t )
    );

  //============================================================================
  // Update temperature and density.  Update dt and time again in case
  // changed in update_zone_properties.
  // Also update Ya, number of atoms per nucleon.
  //============================================================================

    update_zone_properties( zone );

    d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

    d_dt = boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) );

  //============================================================================
  // Evolve abundances.
  //============================================================================

    user::evolve( zone );

  //============================================================================
  // Print out abundances.
  //============================================================================

    if(
       (
         i_step % 
           boost::lexical_cast<int>( zone.getProperty( nnt::s_STEPS ) ) == 0 ||
         d_t >= boost::lexical_cast<double>( zone.getProperty( nnt::s_TEND ) )
       )
    )
    {
      Libnucnet__relabelZone(
        p_my_nucnet,
        zone.getNucnetZone(),
        ( boost::lexical_cast<std::string>( ++k ) ).c_str(),
        NULL,
        NULL
      );
      zone.printAbundances();
      zone.updateProperty(
        nnt::s_YE,
        boost::lexical_cast<std::string>(
          Libnucnet__Zone__computeZMoment( zone.getNucnetZone(), 1 )
        )
      );
      nnt::write_xml( p_my_output, zone.getNucnetZone() );
    }

  //============================================================================
  // Update timestep.
  //============================================================================

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

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

