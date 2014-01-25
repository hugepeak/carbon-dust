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
//! \brief Code for carbon rate functions. 
////////////////////////////////////////////////////////////////////////////////

#include "carbon_rate_functions.h"

//##############################################################################
// register_my_rate_functions()
//############################################################################//

void
register_my_rate_functions(
  Libnucnet__Reac * p_reac
)
{

  //============================================================================
  // Register arrhenius rates, Cherchneff et al. (2009).
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac, 
    S_ARRHENIUS_RATE,
    (Libnucnet__Reaction__userRateFunction) arrhenius_rate_function
  );

  //============================================================================
  // Register arrhenius inverse rates, Clayton PNews3 (2012). 
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac, 
    S_ARRHENIUS_INVERSE_RATE,
    (Libnucnet__Reaction__userRateFunction) arrhenius_inverse_rate_function
  );

  //============================================================================
  // Register carbon condensation rate function and the inverse (2013). 
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac, 
    S_CARBON_CONDENSATION_RATE,
    (Libnucnet__Reaction__userRateFunction) 
      carbon_condensation_rate_function
  );

  Libnucnet__Reac__registerUserRateFunction(
    p_reac, 
    S_CARBON_CONDENSATION_INVERSE_RATE,
    (Libnucnet__Reaction__userRateFunction) 
      carbon_condensation_inverse_rate_function
  );

  //============================================================================
  // Register compton electron rate function(2012). 
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac, 
    S_COMPTON_ELECTRON_RATE,
    (Libnucnet__Reaction__userRateFunction) compton_electron_rate_function
  );

  //============================================================================
  // Register single rate function(2012). 
  //============================================================================

  Libnucnet__Reac__registerUserRateFunction(
    p_reac, 
    S_SINGLE_RATE,
    (Libnucnet__Reaction__userRateFunction) single_rate_function
  );

}

//##############################################################################
// arrhenius_rate_function()
//############################################################################//

double
arrhenius_rate_function(
  Libnucnet__Reaction *p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_rate = 1.;
  const char * s_value;
  nnt::Zone zone = *(nnt::Zone *) p_data;

  //============================================================================
  // Check input.
  //============================================================================

  if( !p_reaction )
  {
    fprintf( stderr, "No reaction supplied.\n" );
    exit( EXIT_FAILURE );
  }

  if( d_t9 <= 0. )
  {
    fprintf( stderr, "Invalid temperature.\n" );
    exit( EXIT_FAILURE );
  }

  s_value =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      "A",
      NULL,
      NULL
    );

  d_rate = atof( s_value );

  s_value =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      "nu",
      NULL,
      NULL
    );

  d_rate *= 
    pow( 
      d_t9 * GSL_CONST_NUM_GIGA / 300., 
      atof( s_value ) 
    ); 

  s_value =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      "Ea",
      NULL,
      NULL
    );

  d_rate *= 
    exp( -atof( s_value ) / ( d_t9 * GSL_CONST_NUM_GIGA ) );

  //============================================================================
  // multiply by N_A for bimolecular reactants. 
  //============================================================================

  int i_count = 0;
  int * p_count = &i_count;

  Libnucnet__Reaction__iterateNuclideReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) count_reactants,
    p_count
  );

  if( i_count == 2 ) {
    d_rate *= 
      (
        GSL_CONST_NUM_AVOGADRO *
        carbon_compute_Ya( zone )
      );
  }

  if( i_count > 2 ) {
    std::cout << "Invalid reactant number for now." << std::endl;
    exit( EXIT_FAILURE );
  }

  return d_rate;

} 

//##############################################################################
// count_reactants()
//############################################################################//

int
count_reactants(
  Libnucnet__Reaction__Element * p_element,
  void * p_data
)
{

  int * p_count = (int *) p_data;

  if( p_element )
    (*p_count)++;

  return 1;

}

//##############################################################################
// arrhenius_inverse_rate_function()
//##############################################################################

double
arrhenius_inverse_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_tau;
  const char * s_value;
  nnt::Zone zone = *(nnt::Zone *) p_data;

  //============================================================================
  // Check input.
  //============================================================================

  if( !p_reaction )
  {
    fprintf( stderr, "No reaction supplied.\n" );
    exit( EXIT_FAILURE );
  }

  if( d_t9 <= 0. )
  {
    fprintf( stderr, "Invalid temperature.\n" );
    exit( EXIT_FAILURE );
  }

  s_value =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      "BondEnergy",
      NULL,
      NULL
    );

  d_tau = 
    1. / 
    arrhenius_rate_function( p_reaction, d_t9, p_data );

  //============================================================================
  // Calculate reduced mass. 
  //============================================================================

  zone.updateProperty( S_INVERSE_REDUCED_MASS, "0." );

  Libnucnet__Reaction__iterateNuclideProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) 
      compute_reduced_mass,
    &zone
  );

  double d_reduced_mass = 
    1. /
    boost::lexical_cast<double>( 
      zone.getProperty( S_INVERSE_REDUCED_MASS) 
    ) * 
    GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS;

  d_tau *=
    pow(
      gsl_pow_2( GSL_CONST_CGSM_PLANCKS_CONSTANT_H ) /
      (
        2. * M_PI * d_reduced_mass *
        GSL_CONST_CGSM_BOLTZMANN *
        d_t9 * GSL_CONST_NUM_GIGA
      ),
      1.5
    ) *
    exp(
      atof( s_value ) * GSL_CONST_CGSM_ELECTRON_VOLT /
      (
        GSL_CONST_CGSM_BOLTZMANN *
        d_t9 * GSL_CONST_NUM_GIGA
      )
    );

  return 1. / d_tau;

}

//##############################################################################
// compute_reduced_mass()
//##############################################################################

int
compute_reduced_mass(
  Libnucnet__Reaction__Element *p_reactant,
  void *p_data
)
{

  nnt::Zone zone = *(nnt::Zone *) p_data;

  double d_inverse_reduced_mass =
    boost::lexical_cast<double>( zone.getProperty( S_INVERSE_REDUCED_MASS ) );

  d_inverse_reduced_mass +=
    1. /
    boost::lexical_cast<double>(
      Libnucnet__Species__getA(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName( p_reactant )
        )
      ) * 16
      -
      Libnucnet__Species__getZ(
        Libnucnet__Nuc__getSpeciesByName(
          Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
          Libnucnet__Reaction__Element__getName( p_reactant )
        )
      ) * 4 
    );

  zone.updateProperty( 
    S_INVERSE_REDUCED_MASS, 
    boost::lexical_cast<std::string>( d_inverse_reduced_mass )
  );

  return 1;

}
  
//##############################################################################
// compute_bond_energy()
//##############################################################################

int
compute_bond_energy(
  Libnucnet__Reaction__Element *p_reactant,
  void *p_data
)
{

  nnt::Zone zone = *(nnt::Zone *) p_data;

  double d_bond_energy =
    boost::lexical_cast<double>( zone.getProperty( "bond energy" ) );

  d_bond_energy *=
    ( 
      1. +
      boost::lexical_cast<double>(
        Libnucnet__Species__getA(
          Libnucnet__Nuc__getSpeciesByName(
            Libnucnet__Net__getNuc(
              Libnucnet__Zone__getNet( zone.getNucnetZone() )
            ),
            Libnucnet__Reaction__Element__getName( p_reactant )
          )
        )
      ) * 1.e-1 
    );

  zone.updateProperty( 
    "bond energy", 
    boost::lexical_cast<std::string>( d_bond_energy )
  );

  return 1;

}
  
//##############################################################################
// carbon_condensation_rate_function()
//##############################################################################

double
carbon_condensation_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_rate, d_Rc = 1.7e-8;
  nnt::Zone zone = *(nnt::Zone *) p_data;

  //============================================================================
  // Check input.
  //==========================================================================//

  if( !p_reaction )
  {
    fprintf( stderr, "No reaction supplied.\n" );
    exit( EXIT_FAILURE );
  }

  if( d_t9 <= 0. )
  {
    fprintf( stderr, "Invalid temperature.\n" );
    exit( EXIT_FAILURE );
  }

  zone.updateProperty( S_CARBON_NUMBER, "1" );

  Libnucnet__Reaction__iterateNuclideReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) carbon_number_callback,
    p_data
  );

  Libnucnet__Reaction__iterateNuclideProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) carbon_number_callback,
    p_data
  );

  double d_carbon_number =
      boost::lexical_cast<double>( zone.getProperty( S_CARBON_NUMBER ) );
 
  d_rate = 
    pow( 
      d_carbon_number,
      2./3. 
    ) *
    M_PI *
    gsl_pow_2( d_Rc ) *
    sqrt(
      3. *
      GSL_CONST_CGSM_BOLTZMANN *
      d_t9 * GSL_CONST_NUM_GIGA /
      (
        12. / GSL_CONST_NUM_AVOGADRO
      )
    );

  d_rate *=
    (
      GSL_CONST_NUM_AVOGADRO *
      carbon_compute_Ya( zone )
    );

  return d_rate;
  
}

//##############################################################################
// carbon_number_callback()
//############################################################################//

int
carbon_number_callback(
  Libnucnet__Reaction__Element * p_element,
  void * p_data
)
{

  nnt::Zone zone = *(nnt::Zone *) p_data;

  unsigned int i_Z = 
    Libnucnet__Species__getZ(
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc(
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
        Libnucnet__Reaction__Element__getName( p_element )
      )
    );

  if( 
    i_Z >
    boost::lexical_cast<unsigned int>( 
      zone.getProperty( S_CARBON_NUMBER )
    )
  )
    zone.updateProperty( 
      S_CARBON_NUMBER, 
      boost::lexical_cast<std::string>( i_Z - 1 )
    );

  return 1;
           
}

//##############################################################################
// carbon_condensation_inverse_rate_function()
//##############################################################################

double
carbon_condensation_inverse_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_tau;
  const char * s_value;
  nnt::Zone zone = *(nnt::Zone *) p_data;

  //============================================================================
  // Check input.
  //============================================================================

  if( !p_reaction )
  {
    fprintf( stderr, "No reaction supplied.\n" );
    exit( EXIT_FAILURE );
  }

  if( d_t9 <= 0. )
  {
    fprintf( stderr, "Invalid temperature.\n" );
    exit( EXIT_FAILURE );
  }

  s_value =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      "BondEnergy",
      NULL,
      NULL
    );

  zone.updateProperty( "bond energy", s_value );

  Libnucnet__Reaction__iterateNuclideReactants(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) 
      compute_bond_energy,
    &zone
  );

  double d_bond_energy = 
    boost::lexical_cast<double>( zone.getProperty( "bond energy" ) );

  double d_condense_rate =
    carbon_condensation_rate_function( p_reaction, d_t9, p_data );

  d_tau = 
    1. / 
    d_condense_rate;

  //============================================================================
  // Calculate reduced mass. 
  //============================================================================

  zone.updateProperty( S_INVERSE_REDUCED_MASS, "0." );

  Libnucnet__Reaction__iterateNuclideProducts(
    p_reaction,
    (Libnucnet__Reaction__Element__iterateFunction) 
      compute_reduced_mass,
    &zone
  );

  double d_reduced_mass = 
    1. /
    boost::lexical_cast<double>( 
      zone.getProperty( S_INVERSE_REDUCED_MASS) 
    ) * 
    GSL_CONST_CGSM_UNIFIED_ATOMIC_MASS;

  d_tau *=
    pow(
      gsl_pow_2( GSL_CONST_CGSM_PLANCKS_CONSTANT_H ) /
      (
        2. * M_PI * d_reduced_mass *
        GSL_CONST_CGSM_BOLTZMANN *
        d_t9 * GSL_CONST_NUM_GIGA
      ),
      1.5
    ) *
    exp(
      d_bond_energy * GSL_CONST_CGSM_ELECTRON_VOLT /
      (
        GSL_CONST_CGSM_BOLTZMANN *
        d_t9 * GSL_CONST_NUM_GIGA
      )
    );

  return 1. / d_tau;

}

//##############################################################################
// compton_electron_rate_function()
//############################################################################//

double
compton_electron_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_tau;
  const char * s_value;
  nnt::Zone zone = *(nnt::Zone *) p_data;

  //============================================================================
  // Check input.
  //============================================================================

  if( !p_reaction )
  {
    fprintf( stderr, "No reaction supplied.\n" );
    exit( EXIT_FAILURE );
  }

  if( d_t9 <= 0. )
  {
    fprintf( stderr, "Invalid temperature.\n" );
    exit( EXIT_FAILURE );
  }

  s_value =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      "tau",
      NULL,
      NULL
    );

  d_tau = 
    atof( s_value ) *
    exp(
      (
        boost::lexical_cast<double>( zone.getProperty( nnt::s_TIME ) ) -
        1.e6
      ) /
      9.59e6 
    );

  return 1. / d_tau;
  
}

//##############################################################################
// single_rate_function()
//############################################################################//

double
single_rate_function(
  Libnucnet__Reaction * p_reaction,
  double d_t9,
  void *p_data
)
{

  double d_rate;
  const char * s_value;

  //============================================================================
  // Check input.
  //============================================================================
  
  if( !p_reaction )
  {
    fprintf( stderr, "No reaction supplied.\n" );
    exit( EXIT_FAILURE );
  }

  if( d_t9 <= 0. )
  {
    fprintf( stderr, "Invalid temperature.\n" );
    exit( EXIT_FAILURE );
  }

  if( p_data )
  {
    fprintf( stderr, "No extra data needed in single rate.\n" );
    exit( EXIT_FAILURE );
  }

  s_value =
    Libnucnet__Reaction__getUserRateFunctionProperty(
      p_reaction,
      "rate",
      NULL,
      NULL
    );

  d_rate = atof( s_value );

  return d_rate;
  
}

//##############################################################################
// update_my_rate_functions_data().
//##############################################################################

void
update_my_rate_functions_data(
  nnt::Zone &zone
)
{

  Libnucnet__Zone__updateDataForUserRateFunction(
    zone.getNucnetZone(),
    S_ARRHENIUS_RATE,
    &zone
  );

  Libnucnet__Zone__updateDataForUserRateFunction(
    zone.getNucnetZone(),
    S_ARRHENIUS_INVERSE_RATE,
    &zone
  );

  Libnucnet__Zone__updateDataForUserRateFunction(
    zone.getNucnetZone(),
    S_CARBON_CONDENSATION_RATE,
    &zone
  );

  Libnucnet__Zone__updateDataForUserRateFunction(
    zone.getNucnetZone(),
    S_CARBON_CONDENSATION_INVERSE_RATE,
    &zone
  );

  Libnucnet__Zone__updateDataForUserRateFunction(
    zone.getNucnetZone(),
    S_COMPTON_ELECTRON_RATE,
    &zone
  );

}
  
//##############################################################################
// carbon_compute_Ya().
//##############################################################################

double
carbon_compute_Ya(
  nnt::Zone &zone
)
{

  if( zone.hasProperty( S_YA ) )
    return boost::lexical_cast<double>( zone.getProperty( S_YA ) );

  zone.updateProperty( S_NUCLEON_NUMBER_PER_ATOM, "0" );

  Libnucnet__Nuc__iterateSpecies(
    Libnucnet__Net__getNuc( Libnucnet__Zone__getNet( zone.getNucnetZone() ) ),
    (Libnucnet__Species__iterateFunction) carbon_compute_Ya_callback,
    &zone 
  );

  double d_ya =
    1. / 
    boost::lexical_cast<double>( 
      zone.getProperty( S_NUCLEON_NUMBER_PER_ATOM ) 
    ); 

  zone.updateProperty( S_YA, boost::lexical_cast<std::string>( d_ya ) );

  return d_ya;

}
  
//##############################################################################
// carbon_compute_Ya_callback().
//##############################################################################

int
carbon_compute_Ya_callback(
  Libnucnet__Species * p_species,
  void * p_data 
)
{

  nnt::Zone zone = *(nnt::Zone *) p_data;

  double d_species_nucleon_number =
    (
      boost::lexical_cast<double>(
        Libnucnet__Species__getZ( p_species )
      ) * 12.
      +
      boost::lexical_cast<double>(
        Libnucnet__Species__getA( p_species ) -
        Libnucnet__Species__getZ( p_species )
      ) * 16.
    );
 
  double d_nucleon_number = 
    boost::lexical_cast<double>( 
      zone.getProperty( S_NUCLEON_NUMBER_PER_ATOM ) 
    );

  d_nucleon_number +=
    Libnucnet__Zone__getSpeciesAbundance(
      zone.getNucnetZone(),
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc( 
          Libnucnet__Zone__getNet( zone.getNucnetZone() )
        ),
        Libnucnet__Species__getName( p_species )
      )
    ) *
      d_species_nucleon_number;

  zone.updateProperty(
    S_NUCLEON_NUMBER_PER_ATOM,
    boost::lexical_cast<std::string>( d_nucleon_number )
  );
      
  return 1;

}
