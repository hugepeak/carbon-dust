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
      stderr, "  out_file = output data hdf5 filename\n\n"
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
    zone.getProperty<std::string>( nnt::s_RHO_0 )
  );

  zone.updateProperty(
    nnt::s_T9,
    zone.getProperty<std::string>( nnt::s_T9_0 )
  );

  update_He_abundances( zone );

}

//##############################################################################
// update_zone_properties().
//##############################################################################

void
update_zone_properties( nnt::Zone& zone )
{

  double d_t, d_tau;

  d_t = zone.getProperty<double>( nnt::s_TIME );

  try
  {
    d_tau = zone.getProperty<double>( nnt::s_TAU );
  }
  catch( const boost::bad_lexical_cast& e )
  {
    if( zone.getProperty<std::string>( nnt::s_TAU ) == "inf" )
    {
      return;
    }
    else
    {
      std::cerr << "Invalid tau." << std::endl;
    }
  }

  //============================================================================
  // update density and temperature.
  //============================================================================

  zone.updateProperty(
    nnt::s_RHO,
    zone.getProperty<double>( nnt::s_RHO_0 ) 
      * 
      pow( d_tau / d_t, 3. )
  );

  zone.updateProperty(
    nnt::s_T9,
    zone.getProperty<double>( nnt::s_T9_0 ) 
      * 
      pow( d_tau / d_t, 1. )
  );

  //============================================================================
  // update density and temperature if they are too small.
  //============================================================================

  if(
    zone.getProperty<double>( nnt::s_T9 ) < 1.e-10
  )
    zone.updateProperty(
      nnt::s_T9,
      boost::lexical_cast<std::string>( 1.e-10 )
    );

  if(
    zone.getProperty<double>( nnt::s_RHO ) < 1.e-30
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

void 
update_He_abundances( nnt::Zone& zone )
{

  //============================================================================
  // Fix He and He+ abundances.
  // Atomic number A for He is 4 in Libnucnet.
  // When Libnucnet calculates the abundances, it uses mass fraction divided
  // by A. While here the mass fraction is really atom fraction, that is
  // the atom number in this kind of molecule divided by the total atoms number.
  // The code is designed to only have C and/or O in the molecules.
  // Z is C number and N is O number.
  // A is the sum of C and O number in a molecule.
  // So the abundance here is the atom fraction divided by the atom in the
  // molecule, which is the number of molecule per atom.
  // For He and He+, there is only one atom in the molecule, so we need to
  // multiply 4 back to its abundance.
  //============================================================================

  const char * species_list[2] = {"he4g","he4+"}; 

  BOOST_FOREACH( const char * species, species_list )
  {
 
    if( 
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
        species
      )
    ) {
   
      Libnucnet__Zone__updateSpeciesAbundance(
        zone.getNucnetZone(),
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          species
        ),
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc(
              Libnucnet__Zone__getNet( zone.getNucnetZone() )
            ),
            species
          )
        ) * 4.
      );

    }

  }

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


