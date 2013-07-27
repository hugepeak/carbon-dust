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
//! \brief Carbon molecule utility code. 
////////////////////////////////////////////////////////////////////////////////

#include "my_network_utilities.h"

namespace my_user
{

//##############################################################################
// add_default_molecules_to_nuc().
//##############################################################################

void
add_default_molecules_to_nuc( 
  Libnucnet__Nuc * p_nuc
)
{

  add_molecules_to_nuc( p_nuc, 1, (unsigned int) 1e1 );

}

//##############################################################################
// add_molecules_to_nuc().
//##############################################################################

void
add_molecules_to_nuc(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  gsl_vector * p_t9, * p_log10_partf;
  std::string s_state;

  p_t9 = gsl_vector_alloc( 1 );
  p_log10_partf = gsl_vector_alloc( 1 );

  gsl_vector_set( p_t9, 0, 1. );
  gsl_vector_set( p_log10_partf, 0, 0. );

  //============================================================================
  // Add Carbon and neutron. 
  //============================================================================

  s_state = "";

  for( unsigned int i = i_lo; i <= i_hi; i++ )
    add_molecule_to_nuc( p_nuc, i, i, s_state, p_t9, p_log10_partf );

  //add_molecule_to_nuc( p_nuc, 0, 1, s_state, p_t9, p_log10_partf );

  gsl_vector_free( p_t9 );
  gsl_vector_free( p_log10_partf );

}

//##############################################################################
// add_molecule_to_nuc().
//##############################################################################

void
add_molecule_to_nuc(
  Libnucnet__Nuc * p_nuc,
  unsigned int i_Z,
  unsigned int i_A,
  std::string s_state,
  gsl_vector * p_t9,
  gsl_vector * p_log10_parf
)
{
 
  Libnucnet__Species * p_species;

  std::string s_source = "example";
  int i_state = 0;
  double d_mass_excess = 0.;
  double d_spin = 0.;

  i_state = 0;
     
  p_species = 
    Libnucnet__Species__new(
      i_Z,
      i_A,
      s_source.c_str(),
      i_state,
      s_state.c_str(),
      d_mass_excess,
      d_spin,
      p_t9,
      p_log10_parf
    );

  Libnucnet__Nuc__addSpecies( p_nuc, p_species );

}

//##############################################################################
// add_default_reactions_to_net().
//##############################################################################

void
add_default_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  add_reactions_to_net(
    p_net,
    1,
    10
  );

  add_forward_reactions_to_net( p_net, 10, 1e1 );

/*
  add_reactions_to_net(
    p_net,
    1,
    (unsigned int)
      Libnucnet__Nuc__getNumberOfSpecies(
        Libnucnet__Net__getNuc( p_net )
      )
  );
*/

}

//##############################################################################
// add_reactions_to_net().
//##############################################################################

void
add_reactions_to_net(
  Libnucnet__Net * p_net,
  unsigned int i_lo,
  unsigned int i_hi
)
{ 
    
  add_forward_reactions_to_net( p_net, i_lo, i_hi );
  add_reverse_reactions_to_net( p_net, i_lo, i_hi );

}

//##############################################################################
// add_forward_reactions_to_net().
//##############################################################################

void
add_forward_reactions_to_net(
  Libnucnet__Net * p_net,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  ReactionData reaction_data;

  for( unsigned int i = i_lo; i < i_hi; i++ )
  {

    reaction_data.sSource = "Meyer notes 2013";
    reaction_data.sUserRateFunctionKey = "carbon condensation rate";

    reaction_data.vReactants.push_back( std::make_pair( i, i ) );
    reaction_data.vReactants.push_back( std::make_pair( 1, 1 ) );
    reaction_data.vProducts.push_back( std::make_pair( i + 1, i + 1 ) );

    add_reaction_to_net( p_net, reaction_data );

    reaction_data.clear();

  }

}
    
//##############################################################################
// add_reverse_reactions_to_net().
//##############################################################################

void
add_reverse_reactions_to_net(
  Libnucnet__Net * p_net,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  ReactionData reaction_data;

  for( unsigned int i = i_lo; i < i_hi; i++ )
  {

    reaction_data.sSource = "Meyer notes 2013";
    reaction_data.sUserRateFunctionKey = "carbon condensation inverse rate";
    reaction_data.vUserRateData.push_back(
      std::make_pair( "BondEnergy", "5.0" )
    );

    reaction_data.vReactants.push_back( std::make_pair( i + 1, i + 1 ) );
    reaction_data.vProducts.push_back( std::make_pair( i, i ) );
    reaction_data.vProducts.push_back( std::make_pair( 1, 1 ) );

    add_reaction_to_net( p_net, reaction_data );

    reaction_data.clear();
    
  }

}

//##############################################################################
// add_reaction_to_net().
//##############################################################################

void
add_reaction_to_net(
  Libnucnet__Net * p_net,
  ReactionData &reaction_data
)
{

  Libnucnet__Reaction * p_reaction;

  p_reaction = Libnucnet__Reaction__new();

  std::pair<unsigned int,unsigned int> reactant;

  BOOST_FOREACH( reactant, reaction_data.vReactants ) {

    if( 
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_net ),
        reactant.first,
        reactant.second,
        NULL 
      )
    ) 
      Libnucnet__Reaction__addReactant( 
        p_reaction, 
        Libnucnet__Species__getName(
          Libnucnet__Nuc__getSpeciesByZA(
            Libnucnet__Net__getNuc( p_net ),
            reactant.first,
            reactant.second,
            NULL 
          )
        )
      );
    else
    {
      std::cerr << "No " << reactant.first
                << reactant.second << " in network!" << std::endl;
      exit( EXIT_FAILURE );
    }
  
  }

  std::pair<unsigned int,unsigned int> product;

  BOOST_FOREACH( product, reaction_data.vProducts ) {

    if( 
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_net ),
        product.first,
        product.second,
        NULL 
      )
    ) 
      Libnucnet__Reaction__addProduct( 
        p_reaction, 
        Libnucnet__Species__getName(
          Libnucnet__Nuc__getSpeciesByZA(
            Libnucnet__Net__getNuc( p_net ),
            product.first,
            product.second,
            NULL 
          )
        )
      );
    else
    {
      std::cerr << "No " << product.first 
                << product.second << " in network!" << std::endl;
      exit( EXIT_FAILURE );
    }
  
  }

  Libnucnet__Reaction__updateSource( 
    p_reaction, 
    reaction_data.sSource.c_str() 
  );

  Libnucnet__Reaction__setUserRateFunctionKey(
    p_reaction,
    reaction_data.sUserRateFunctionKey.c_str()
  );

  std::pair<std::string,std::string> my_data;

  BOOST_FOREACH( my_data, reaction_data.vUserRateData ) {
 
    Libnucnet__Reaction__updateUserRateFunctionProperty(
       p_reaction,
       my_data.first.c_str(),
       NULL,
       NULL,
       my_data.second.c_str() 
    );

  }

  Libnucnet__Reac__addReaction( 
    Libnucnet__Net__getReac( p_net ), 
    p_reaction 
  ); 

}

//##############################################################################
// ReactionData::~ReactionData(). 
//##############################################################################

ReactionData::~ReactionData()
{

  this->clear();

}

//##############################################################################
// ReactionData::clear(). 
//##############################################################################

void
ReactionData::clear()
{

  this->vReactants.clear();
  this->vProducts.clear();
  this->sSource = "";
  this->sUserRateFunctionKey = "";
  this->vUserRateData.clear();

}
 
//##############################################################################
// update_bin_net(). 
//##############################################################################

int
update_bin_net(
  Libnucnet * p_nucnet
)
{

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
    zone.hasProperty( S_RUN_BIN ) &&
    zone.getProperty( S_RUN_BIN ) == "yes"
  )
  {

    unsigned int i_bin_start = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_BIN_START ) );
    unsigned int i_photon_end = 
      boost::lexical_cast<unsigned int>( zone.getProperty( S_PHOTON_END ) );

    Libnucnet__free( p_nucnet );

    p_nucnet = Libnucnet__new();

    add_molecules_to_nuc(
      Libnucnet__Net__getNuc( Libnucnet__getNet( p_nucnet ) ),
      1, 
      i_bin_start
    );

    add_reactions_to_net(
      Libnucnet__getNet( p_nucnet ),
      1,
      i_photon_end
    ); 
  
    add_forward_reactions_to_net(
      Libnucnet__getNet( p_nucnet ),
      i_photon_end,
      i_bin_start
    ); 

    return 1;

  }

  return 0;

}

} // my_user
