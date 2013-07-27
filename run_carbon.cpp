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

#include "carbon_evolve.h"
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

void
my_print_abundances( nnt::Zone & );

int
print_abundance( Libnucnet__Species *, Libnucnet__Zone * );

int
print_bin_abundances( nnt::Zone & ); 

double
compute_bin_sum( nnt::Zone & );

void
bin_update_timestep(
  nnt::Zone &, double &, double, double, double
);

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

  std::cout << "species number = " <<
    Libnucnet__Nuc__getNumberOfSpecies(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      )
    ) << std::endl;

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

    double d_rho = 
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) );

    if( d_rho < 0. )
      std::cout << d_rho << std::endl;

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
        ( boost::lexical_cast<std::string>( ++k ) ).c_str(),
        NULL,
        NULL
      );

      if( 
        zone.hasProperty( "run bin" ) &&
        zone.getProperty( "run bin" ) == "yes" &&
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

    std::cout << "1 - Xsum_network = " <<
      1. - 
      Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) << "  ";

    double d_bin_mass = compute_bin_sum( zone );

    std::cout << "1 - Xsum = " <<
      1. - 
      Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) -
      //compute_bin_sum( zone ) <<
      d_bin_mass <<
      std::endl << std::endl;

  } else {

    std::cout << "1 - Xsum = " <<
      1. - 
      Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 ) -
      compute_bin_sum( zone ) <<
      std::endl << std::endl;

  }

  Libnucnet__Nuc__setSpeciesCompareFunction(
    Libnucnet__Net__getNuc( 
      Libnucnet__Zone__getNet( zone.getNucnetZone() ) 
    ),
    (Libnucnet__Species__compare_function) nnt::species_sort_function
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
  
  printf( "%16lu%16lu%16.6e%16.6e\n",
    (unsigned long) Libnucnet__Species__getZ( p_species ),
    (unsigned long) Libnucnet__Species__getA( p_species ),
    d_abund,
    Libnucnet__Zone__getSpeciesAbundanceChange(
      p_zone, p_species
    )
  );

  return 1;

}

//############################################################################
// print_bin_abundances()
//############################################################################

int
print_bin_abundances(
  nnt::Zone & zone
)
{

  for(
    int i = 1;
    i <= boost::lexical_cast<int>( zone.getProperty( "bin number" ) );
    i++
  )
  {

    printf( "%16lu%16lu%16.6e%16.6e\n",
      (unsigned long)( 
        boost::lexical_cast<int>( zone.getProperty( "bin start" ) ) *
        pow(
          boost::lexical_cast<int>( zone.getProperty( "bin factor" ) ),
          i - 1
        ) + 1
      ),
      (unsigned long)(
        boost::lexical_cast<int>( zone.getProperty( "bin start" ) ) *
        pow(
          boost::lexical_cast<int>( zone.getProperty( "bin factor" ) ),
          i
        )
      ),
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance", 
          boost::lexical_cast<std::string>( i ) 
        )
      ),
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance change", 
          boost::lexical_cast<std::string>( i ) 
        )
      )
    );

  }

  return 1;

}

//############################################################################
// compute_bin_sum()
//############################################################################

double
compute_bin_sum(
  nnt::Zone & zone
)
{

  double d_result = 0.;

  for(
    size_t i = 1;
    i <= boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) );
    i++
  )
  {

    double d_first, d_last;

    d_first = 
      boost::lexical_cast<double>( zone.getProperty( "bin start" ) ) *
      pow(
        boost::lexical_cast<double>( zone.getProperty( "bin factor" ) ),
        (double) i - 1
      ) + 1.; 

    d_last = 
      boost::lexical_cast<double>( zone.getProperty( "bin start" ) ) *
      pow(
        boost::lexical_cast<double>( zone.getProperty( "bin factor" ) ),
        (double) i
      ); 

    d_result +=
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance", 
          boost::lexical_cast<std::string>( i ) 
        )
      ) *
      ( d_last + d_first ) * ( d_last - d_first + 1. ) / 2.;

  }

  return d_result; 

}

//############################################################################
// bin_update_timestep()
//############################################################################

void
bin_update_timestep(
  nnt::Zone & zone,
  double & d_dt,
  double d_regt,
  double d_regy,
  double d_ymin
)
{

  double d_dtnew = ( 1. + d_regt ) * d_dt;
  double d_check, d_y, d_dy, d_tiny = 1.e-300;

  int i_bin_number = 
    boost::lexical_cast<int>( zone.getProperty( "bin number" ) );

  for( int i = 1; i <= i_bin_number; i++ ) {

    d_y = 
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance",
          boost::lexical_cast<std::string>( i )
        ) 
      );

    d_dy = 
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance change",
          boost::lexical_cast<std::string>( i )
        ) 
      );

    if( d_y > d_ymin && !WnMatrix__value_is_zero( d_dy ) ) {

      d_check =
        d_dt *
        d_regy *
        fabs( d_y / ( d_dy + d_tiny ) );

      if( d_check < d_dtnew ) {

        d_dtnew = d_check;

      }

    }

  }

  d_dt = d_dtnew; 

}

