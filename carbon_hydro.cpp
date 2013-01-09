//////////////////////////////////////////////////////////////////////////////
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
//! \brief Example code for carbon dust hydrodynamics.
////////////////////////////////////////////////////////////////////////////////

#define I_LO  9
#define I_HI  20

#define S_INCLUED_BIG_MOLECULES "no"

#include "carbon_hydro.h"

//##############################################################################
// get_nucnet().
//##############################################################################

Libnucnet *
get_nucnet( int argc, char **argv )
{

  Libnucnet * p_nucnet;

  //============================================================================
  // Check input.
  //============================================================================

  if( argc == 2 && strcmp( argv[1], "--help" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] <<
      " runs a carbon oxygen network calculation for matter expanding" <<
      " with time." << std::endl << std::endl;
    std::cout << "For a usage statement, type " << std::endl << std::endl;
    std::cout << argv[0] << " --usage" << std::endl << std::endl;
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " ../../data_pub/co_net.xml " <<
      "../../data_pub/co_zone.xml my_output.xml \"[z <= 20]\"" <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 4 || argc > 6 || strcmp( argv[1], "--usage" ) == 0 )
  {
    fprintf(
      stderr,
      "\nUsage: %s net_file zone_file out_file xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_file = input network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_file = input single zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  out_file = output data xml filename\n\n"
    );
    fprintf(
      stderr,
      "  xpath_nuc = nuclear xpath expression (optional--required if xpath_reac specified)\n\n"
    );
    fprintf(
      stderr, "  xpath_reac = reaction xpath expression (optional)\n\n"
    );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Validate input net file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__Net__is_valid_input_xml( argv[1] ) ) {
      fprintf( stderr, "Not valid libnucnet net input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Validate input zone file.
  //============================================================================

  if( strcmp( VALIDATE, "yes" ) == 0 )
  {
    if( !Libnucnet__is_valid_zone_data_xml( argv[2] ) ) {
      fprintf( stderr, "Not valid libnucnet zone data input!\n" );
      exit( EXIT_FAILURE );
    }
  }

  //============================================================================
  // Read and store input.
  //============================================================================

  p_nucnet = Libnucnet__new();

  if( argc == 4 )
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      NULL,
      NULL
    );

  }
  else if( argc == 5 )
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      argv[4],
      NULL
    );
  }
  else
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      argv[4],
      argv[5]
    );
  }

  if( strcmp( S_INCLUED_BIG_MOLECULES, "yes" ) == 0 ) {

    unsigned int i_lo = I_LO;
    unsigned int i_hi = I_HI;

    add_carbon_molecules( 
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) ),
      i_lo,
      i_hi
    );

    add_carbon_reactions(
      Libnucnet__getNet( p_nucnet ),
      i_lo,
      i_hi
    );

  }
 
  Libnucnet__assignZoneDataFromXml( p_nucnet, argv[2], NULL );

  return p_nucnet;

}

//##############################################################################
// add_carbon_molecules().
//##############################################################################

void
add_carbon_molecules(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  Libnucnet__Species * p_species;
  gsl_vector * p_t9, * p_log10_parf;

  for( unsigned int i = i_lo; i <= i_hi; i++ ) {

    p_t9 = gsl_vector_alloc( 1 );
    p_log10_parf = gsl_vector_alloc( 1 );

    gsl_vector_set( p_t9, 0, 1. );
    gsl_vector_set( p_log10_parf, 0, 0. );

    p_species = 
      Libnucnet__Species__new(
        i,
        i,
        "example",
        1,
        "C",
        0.,
        0.,
        p_t9,
        p_log10_parf
      );

    Libnucnet__Nuc__addSpecies( p_nuc, p_species );

    p_species = 
      Libnucnet__Species__new(
        i,
        i,
        "example",
        1,
        "R",
        0.,
        0.,
        p_t9,
        p_log10_parf
      );

    Libnucnet__Nuc__addSpecies( p_nuc, p_species );

    gsl_vector_free( p_t9 );
    gsl_vector_free( p_log10_parf );

  } 

}

//##############################################################################
// add_carbon_reactions().
//##############################################################################

void
add_carbon_reactions(
  Libnucnet__Net * p_net,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  Libnucnet__Reaction * p_reaction;

  for( unsigned int i = i_lo; i <= i_hi; i++ ) {

    p_reaction = Libnucnet__Reaction__new();

    Libnucnet__Reaction__setUserRateFunctionKey(
      p_reaction,
      "carbon condensation rate"
    );

    Libnucnet__Reaction__updateSource( p_reaction, "Meyer, notes" );

    //--------------------------------------------------------------------------
    // Add reaction: C_{i-1}^R + C -> C_i^R + gamma.
    //--------------------------------------------------------------------------

    Libnucnet__Reaction__addReactant( 
      p_reaction, 
      Libnucnet__Species__getName(
        Libnucnet__Nuc__getSpeciesByZA(
          Libnucnet__Net__getNuc( p_net ), 
          i - 1, 
          i - 1,
          "R"
        )
      )   
    );

    Libnucnet__Reaction__addReactant( p_reaction, "h1" );

    Libnucnet__Reaction__addProduct(
      p_reaction, 
      Libnucnet__Species__getName(
        Libnucnet__Nuc__getSpeciesByZA(
          Libnucnet__Net__getNuc( p_net ), 
          i, 
          i,
          "R"
        )
      )   
    );

    Libnucnet__Reaction__addProduct( p_reaction, "gamma" );

    Libnucnet__Reaction__updateUserRateFunctionProperty(
       p_reaction,
       "CarbonNumber",
       NULL,
       NULL,
       boost::lexical_cast<std::string>( i - 1 ).c_str()
    );

    Libnucnet__Reac__addReaction( 
      Libnucnet__Net__getReac( p_net ), 
      p_reaction 
    ); 

  }

}

//##############################################################################
// initialize_zone().
//##############################################################################

void
initialize_zone( nnt::Zone& zone, char ** argv )
{

  //============================================================================
  // Set output file.
  //============================================================================

  Libnucnet__Zone__updateProperty(
    zone.getNucnetZone(),
    S_OUTPUT_XML_FILE,
    NULL,
    NULL,
    argv[3]
  );

  //============================================================================
  // Set the initial evolution network.
  //============================================================================

  zone.updateProperty(
    nnt::s_RHO,
    zone.getProperty( nnt::s_RHO_0 )
  );

  zone.updateProperty(
    nnt::s_T9,
    zone.getProperty( nnt::s_T9_0 )
  );

  zone.updateProperty(
    nnt::s_TIME,
    zone.getProperty( nnt::s_TIME )
  );

}

//##############################################################################
// update_zone_properties().
//##############################################################################

void
update_zone_properties( nnt::Zone& zone )
{

  double d_t, d_tau;

  d_t = boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) );

  try
  {
    d_tau = boost::lexical_cast<double>( zone.getProperty( nnt::s_TAU ) );
  }
  catch( const boost::bad_lexical_cast& e )
  {
    if( zone.getProperty( nnt::s_TAU ) == "inf" )
    {
      return;
    }
    else
    {
      std::cerr << "Invalid tau." << std::endl;
    }
  }

  zone.updateProperty(
    nnt::s_RHO,
    boost::lexical_cast<std::string>(
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_RHO_0 ) 
      ) * 
      pow( d_tau / d_t, 3. )
    )
  );

  zone.updateProperty(
    nnt::s_T9,
    boost::lexical_cast<std::string>(
      boost::lexical_cast<double>(
        zone.getProperty( nnt::s_T9_0 ) 
      ) * 
      pow( 
        d_tau / d_t, 
        boost::lexical_cast<double>( zone.getProperty( S_POWER ) )
      )
    )
  );

  if(
    boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) < 1.e-10
  )
    zone.updateProperty(
      nnt::s_T9,
      boost::lexical_cast<std::string>( 1.e-10 )
    );

  if(
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) < 1.e-30
  )
    zone.updateProperty(
      nnt::s_RHO,
      boost::lexical_cast<std::string>( 1.e-30 )
    );

}

//##############################################################################
// set_zone().
//##############################################################################

int
set_zone( Libnucnet * p_nucnet, nnt::Zone& zone, char ** argv )
{

  if( !argv )
  {
    std::cerr << "Invalid input." << std::endl;
    return 0;
  }

  if( !Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" ) )
  {
    std::cerr << "Zone not found!" << std::endl;
    return 0;
  }

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" )
  );

  return 1;

}

