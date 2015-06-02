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
//! \brief Example code for running a network calculation of carbon and oxygen.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <Libnucnet.h>

#include "user/hdf5_routines.h"

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
#define S_WIDTH        "1"          // Solver type: ARROW or GSL

#define S_TOTAL_ATOM_NUMBER    "total atom number"
#define S_XSEC                 "sigma v"
#define S_NUMBER_GRAINS        "number of grains"

#define I_MAX          8

//##############################################################################
// The grain.
//##############################################################################

struct Grain
{

  double Nprev, Nnext;
  double FormationTime;

  Grain( double _t, double _N ) :
    Nprev( _N ), Nnext( _N ), FormationTime( _t ) {}

};

//##############################################################################
// grain_function().
//##############################################################################

void
grain_function(
  nnt::Zone& zone,
  std::vector<Grain>& grains,
  WnMatrix * p_matrix,
  gsl_vector * p_rhs
)
{

  Libnucnet__Species * p_atom =
    Libnucnet__Nuc__getSpeciesByZA(
      Libnucnet__Net__getNuc(
        Libnucnet__Zone__getNet( zone.getNucnetZone() )
      ),
      1,
      1,
      "g"
    );

  size_t i_index = Libnucnet__Species__getIndex( p_atom );

  for( size_t i = 0; i < grains.size(); i++ )
  {

    double delta_t =
      GSL_MIN(
        zone.getProperty<double>( nnt::s_DTIME ),
        zone.getProperty<double>( nnt::s_TIME )
        -
        grains[i].FormationTime
      );

    grains[i].Nnext =
      pow(
        pow( grains[i].Nprev, 1. / 3. )
        +
        3. * zone.getProperty<double>( nnt::s_RHO ) * GSL_CONST_NUM_AVOGADRO *
        Libnucnet__Zone__getSpeciesAbundance(
          zone.getNucnetZone(),
          p_atom
        )
        *
        zone.getProperty<double>( S_XSEC )
        *
        delta_t,
        3.
      );

    gsl_vector_set(
      p_rhs,
      i_index,
      gsl_vector_get( p_rhs, i_index )
      -
      ( grains[i].Nnext - grains[i].Nprev ) /
      (
        zone.getProperty<double>( S_TOTAL_ATOM_NUMBER ) *
        zone.getProperty<double>( nnt::s_DTIME )
      )
    );

  }

}

//##############################################################################
// Xsum_eff().
//##############################################################################

double
Xsum_eff( nnt::Zone& zone, std::vector<Grain>& grains )
{

  double d_xsum_eff = 0;
  for( size_t i = 0; i < grains.size(); i++ )
  {
    d_xsum_eff += grains[i].Nnext /
      zone.getProperty<double>( S_TOTAL_ATOM_NUMBER );
  }
  d_xsum_eff +=  Libnucnet__Zone__computeAMoment( zone.getNucnetZone(), 1 );

  return d_xsum_eff;

}

//##############################################################################
// Xsum_check().
//##############################################################################

bool
Xsum_check( nnt::Zone& zone, std::vector<Grain>& grains )
{

  if( fabs( 1 - Xsum_eff( zone, grains ) ) < 1.e-12 )
    return true;
  else
    return false;

}

//############################################################################
// compute_number_of_new_grains().
//############################################################################

size_t
compute_number_of_new_grains(
  nnt::Zone& zone,
  Libnucnet__Species * p_last
)
{

  return
    static_cast<size_t>(
      Libnucnet__Zone__getSpeciesAbundance( zone.getNucnetZone(), p_last ) *
      zone.getProperty<double>( S_TOTAL_ATOM_NUMBER )
    );

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
// main().
//##############################################################################

int main( int argc, char * argv[] ) {

  int i_step;
  double d_t, d_dt;
  Libnucnet *p_my_nucnet = NULL;
  nnt::Zone zone;
  std::vector<Grain> grains;

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
    d_dt = zone.getProperty<double>( nnt::s_DTIME );
  else
  {
    d_dt = D_DT0;
    zone.updateProperty(
      nnt::s_DTIME,
      d_dt
    );
  }

  if( zone.hasProperty( nnt::s_TIME ) )
    d_t = zone.getProperty<double>( nnt::s_TIME );
  else
  {
    d_t = 0;
    zone.updateProperty(
      nnt::s_TIME,
      d_t
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
      (Libnucnet__Species__compare_function) my_species_sort_function
    );

    Libnucnet__Nuc__sortSpecies(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ) 
    );

    zone.updateProperty( nnt::s_SOLVER, nnt::s_ARROW );

    zone.updateProperty( nnt::s_ARROW_WIDTH, S_WIDTH );

  }

  //============================================================================
  // Create the output hdf5 file.  Do this after nuclide sort.
  //============================================================================

  user::hdf5::create_output(
    argv[3],
    p_my_nucnet
  );

  //============================================================================
  // Toggle off the reverse rates for reactions from detailed balance. 
  //============================================================================

  Libnucnet__Zone__toggleReverseRateDetailedBalance(
    zone.getNucnetZone(), "off"
  );

  //============================================================================
  // Functions.
  //============================================================================

  zone.updateFunction(
    nnt::s_MATRIX_MODIFICATION_FUNCTION,
    static_cast<boost::function<void( WnMatrix *, gsl_vector * )> >(
      boost::bind(
        grain_function,
        boost::ref( zone ),
        boost::ref( grains ),
        _1,
        _2
      )
    )
  );

  zone.updateFunction(
    nnt::s_SAFE_EVOLVE_CHECK_FUNCTION,
    static_cast<boost::function<bool()> >(
      boost::bind(
        Xsum_check, 
        boost::ref( zone ),
        boost::ref( grains )
      )
    )
  );

  //============================================================================
  // Evolve network while t < final t.
  //============================================================================

  i_step = 0;

  while ( 
    d_t < zone.getProperty<double>( nnt::s_TEND )
  )
  {

    d_t += d_dt;

    //==========================================================================
    // Set dt and t.
    //==========================================================================

    zone.updateProperty(
      nnt::s_DTIME,
      d_dt
    );

    zone.updateProperty(
      nnt::s_TIME,
      d_t
    );

    //==========================================================================
    // Update temperature and density.
    //==========================================================================

    update_zone_properties( zone );

    double d_rho = zone.getProperty<double>( nnt::s_RHO );

    if( d_rho < 0. ) {
      std::cout << "rho is " << d_rho << std::endl;
      exit(EXIT_FAILURE);
    }

    //==========================================================================
    // Evolve abundances.
    //==========================================================================

    user::safe_evolve( zone, d_dt, d_dt );

    //==========================================================================
    // Print out abundances.
    //==========================================================================

    if(
       (
         i_step % zone.getProperty<int>( nnt::s_STEPS ) == 0 ||
         d_t >= zone.getProperty<double>( nnt::s_TEND )
       )
    )
    {
      user::hdf5::append_zones(
        argv[3],
        p_my_nucnet
      );
      nnt::print_zone_abundances( zone );
      std::cout << "1 - Xsum_eff = " <<
        1 - Xsum_eff( zone, grains ) << std::endl;
    }

    //==========================================================================
    // Update old grains.
    //==========================================================================

    for( size_t i = 0; i < grains.size(); i++ )
    {
      grains[i].Nprev = grains[i].Nnext;
    }

    //==========================================================================
    // Add new grains.
    //==========================================================================

    Libnucnet__Species * p_last =
      Libnucnet__Nuc__getSpeciesByZA(
	Libnucnet__Net__getNuc( Libnucnet__getNet( p_my_nucnet ) ),
	I_MAX,
	I_MAX,
	"r"
      );

    size_t i_new = compute_number_of_new_grains( zone, p_last );

    for( size_t i = 0; i < i_new; i++ )
    {
      grains.push_back( Grain( d_t, I_MAX ) );
    }


    Libnucnet__Zone__updateSpeciesAbundance(
      zone.getNucnetZone(),
      p_last,
      Libnucnet__Zone__getSpeciesAbundance(
	zone.getNucnetZone(),
	p_last
      ) -
      static_cast<double>( i_new ) /
	zone.getProperty<double>( S_TOTAL_ATOM_NUMBER )
    );

    zone.updateProperty(
      S_NUMBER_GRAINS,
      grains.size()
    );
      
    //==========================================================================
    // Update timestep.
    //==========================================================================

    d_dt = zone.getProperty<double>( nnt::s_DTIME );

    Libnucnet__Zone__updateTimeStep(
      zone.getNucnetZone(),
      &d_dt,
      D_REG_T,
      D_REG_Y,
      D_Y_MIN_DT
    );

    if( zone.getProperty<double>( nnt::s_T9 ) > 10. )
      nnt::normalize_zone_abundances( zone );

    if(
      d_t + d_dt > zone.getProperty<double>( nnt::s_TEND )
    )
    {
      d_dt = zone.getProperty<double>( nnt::s_TEND ) - d_t;
    }

    i_step++;

  }  

  //============================================================================
  // Print out grains.
  //============================================================================

  for( size_t i = 0; i < grains.size(); i++ )
  {
    std::cout << i << "  " << grains[i].FormationTime <<
      "  " << grains[i].Nnext << std::endl;
  }

  //============================================================================
  // Clean up and exit.
  //============================================================================

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}

