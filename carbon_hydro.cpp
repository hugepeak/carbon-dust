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
//! \brief Example code for carbon dust hydrodynamics.
////////////////////////////////////////////////////////////////////////////////

#include "carbon_hydro.h"

//##############################################################################
// Code for default expansion: rho(t) = rho_0 * (t/t_initial)^(-3).
//##############################################################################

#ifdef HYDRO_DEFAULT

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
    std::cout << argv[0] << " data/my_net.xml " <<
      "data/zone.xml out " <<
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

  Libnucnet__assignZoneDataFromXml( p_nucnet, argv[2], NULL );

  //============================================================================
  // Create decade molecules and reactions if desired. 
  //============================================================================

  nnt::Zone zone;

  if( !Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" ) )
  {
    std::cerr << "Zone not found!" << std::endl;
    return 0;
  }

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" )
  );

  if( 
    zone.hasProperty( S_RUN_DECADE ) &&
    zone.getProperty( S_RUN_DECADE ) == "yes"
  )
  {

    unsigned int i_exp_start = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_EXP_START ) );
    unsigned int i_exp = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_EXP ) );
    unsigned int i_base = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_BASE ) );
    unsigned int i_photon_end = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_PHOTON_END ) );

    Libnucnet__free( p_nucnet );

    p_nucnet = Libnucnet__new();

    my_user::add_molecules_to_nuc(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) ),
      1, 
      i_exp_start
    );

    my_user::add_reactions_to_net(
      Libnucnet__getNet( p_nucnet ),
      1,
      i_photon_end
    ); 
  
    my_user::add_forward_reactions_to_net(
      Libnucnet__getNet( p_nucnet ),
      i_photon_end,
      i_exp_start
    ); 

    for( 
      unsigned int i = 1; 
      i <= i_exp; 
      i++
    )
    {
  
      unsigned int i_index =
        i_exp_start *
        (unsigned int) pow( 
          i_base,
          i
        ); 
        
      my_user::add_molecules_to_nuc(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) ),
        i_index,
        i_index
      );
  
    }

    Libnucnet__assignZoneDataFromXml( p_nucnet, argv[2], NULL );

    return p_nucnet;

  }
    
  return p_nucnet;

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

#endif // HYDRO_DEFAULT

//##############################################################################
// Code for expansion from ascii trajectory file. 
//##############################################################################

#ifdef HYDRO_TRAJ

//##############################################################################
// globals.
//##############################################################################

std::vector<double> v_time, v_t9, v_log10_rho;

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
      std::endl;
    std::cout <<
      "according to data from a trajectory file." << std::endl << std::endl;
    std::cout << "For a usage statement, type " << std::endl << std::endl;
    std::cout << argv[0] << " --usage" << std::endl << std::endl;
    std::cout << "For an example usage, type " << std::endl << std::endl;
    std::cout << argv[0] << " --example" << std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc == 2 && strcmp( argv[1], "--example" ) == 0 )
  {
    std::cout << std::endl;
    std::cout << argv[0] << " data/my_net.xml data/zone_traj.xml " << 
      "data/traj.txt out " <<
      std::endl << std::endl;
    exit( EXIT_FAILURE );
  }

  if( argc < 5 || argc > 7 || strcmp( argv[1], "--usage" ) == 0 )
  {
    fprintf(
      stderr,
      "\nUsage: %s net_file zone_file traj_file out_file xpath_nuc xpath_reac\n\n",
      argv[0]
    );
    fprintf(
      stderr, "  net_file = input network data xml filename\n\n"
    );
    fprintf(
      stderr, "  zone_file = input single zone data xml filename\n\n"
    );
    fprintf(
      stderr, "  traj_file = trajectory text file\n\n"
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
  // Create nucnet. 
  //============================================================================

  p_nucnet = Libnucnet__new();

  if( argc == 5 )
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      NULL,
      NULL
    );
  }
  else if( argc == 6 )
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      argv[5],
      NULL
    );
  }
  else
  {
    Libnucnet__Net__updateFromXml(
      Libnucnet__getNet( p_nucnet ),
      argv[1],
      argv[5],
      argv[6]
    );
  }


  //============================================================================
  // Get zone data.
  //============================================================================

  Libnucnet__assignZoneDataFromXml( p_nucnet, argv[2], NULL );

  //============================================================================
  // Get trajectory data.
  //============================================================================

  get_trajectory_data( argv[3] );

  //============================================================================
  // Create decade molecules and reactions if desired. 
  //============================================================================

  nnt::Zone zone;

  if( !Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" ) )
  {
    std::cerr << "Zone not found!" << std::endl;
    return 0;
  }

  zone.setNucnetZone(
    Libnucnet__getZoneByLabels( p_nucnet, "0", "0", "0" )
  );

  if( 
    zone.hasProperty( S_RUN_DECADE ) &&
    zone.getProperty( S_RUN_DECADE ) == "yes"
  )
  {

    unsigned int i_exp_start = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_EXP_START ) );
    unsigned int i_exp = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_EXP ) );
    unsigned int i_base = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_BASE ) );
    unsigned int i_photon_end = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_PHOTON_END ) );

    Libnucnet__free( p_nucnet );

    p_nucnet = Libnucnet__new();

    my_user::add_molecules_to_nuc(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) ),
      1, 
      i_exp_start
    );

    my_user::add_reactions_to_net(
      Libnucnet__getNet( p_nucnet ),
      1,
      i_photon_end
    ); 
  
    my_user::add_forward_reactions_to_net(
      Libnucnet__getNet( p_nucnet ),
      i_photon_end,
      i_exp_start
    ); 

    for( 
      unsigned int i = 1; 
      i <= i_exp; 
      i++
    )
    {
  
      unsigned int i_index =
        i_exp_start *
        (unsigned int) pow( 
          i_base,
          i
        ); 
        
      my_user::add_molecules_to_nuc(
        Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) ),
        i_index,
        i_index
      );
  
    }

    Libnucnet__assignZoneDataFromXml( p_nucnet, argv[2], NULL );

    return p_nucnet;

  }
    
  //============================================================================
  // Done.
  //============================================================================

  return p_nucnet;

}

//##############################################################################
// get_trajectory_data().
//##############################################################################

void
get_trajectory_data( char *s_file )
{

  std::ifstream my_file;
  int i_x0;
  double d_x1, d_x2, d_x3;

  //============================================================================
  // Open file thermodynamics file.
  //============================================================================

  my_file.open( s_file );

  if( my_file.bad() )
  {
    fprintf( stderr, "Couldn't open file.\n" );
    exit( EXIT_FAILURE );
  }

  //============================================================================
  // Assign vectors.
  //============================================================================

  double d_factor_t = 1.e0;
  double d_factor_rho = 1.e0;

  while( my_file >> i_x0 >> d_x1 >> d_x2 >> d_x3 )
  {
    v_time.push_back( d_x1 * d_factor_t );
    v_t9.push_back( d_x2 );
    v_log10_rho.push_back( log10( d_x3 * d_factor_rho ) );
  }

  my_file.close();

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
    argv[4]
  );

  user::update_t9_rho_in_zone_by_interpolation(
    zone,
    "spline",
    v_time,
    v_t9,
    v_log10_rho
  );

  return;

}

//##############################################################################
// update_zone_properties().
//##############################################################################

void
update_zone_properties(
  nnt::Zone& zone
)
{

  user::update_t9_rho_in_zone_by_interpolation(
    zone,
    "spline",
    v_time,
    v_t9,
    v_log10_rho
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

#endif // HYDRO_TRAJ

