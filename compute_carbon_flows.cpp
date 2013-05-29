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
//! \brief Example code to compute the flows for a given zone (chosen by
//!    first label) in a network xml file.
////////////////////////////////////////////////////////////////////////////////

#include <Libnucnet.h>
#include "nnt/auxiliary.h"
#include "nnt/iter.h"
#include "nnt/string_defs.h"
#include "nnt/two_d_weak_rates.h"
#include "nnt/weak_detailed_balance.h"

#include "user/my_user_rate_functions.h"
#include "user/my_flow_utilities.h"
#include "carbon_rate_functions.h"

int main( int argc, char * argv[] ) {

  Libnucnet *p_my_nucnet;
  Libnucnet__NetView * p_net_view;
  std::pair<double,double> flows;

  /*============================================================================
  // Check input.
  //==========================================================================*/

   if ( argc < 3 || argc > 4 ) {
      fprintf(
        stderr, "\nUsage: %s xml_filename zone_xpath reac_xpath\n\n", argv[0]
      );
      fprintf(
        stderr, "  xml_filename = input xml filename\n\n"
      );
      fprintf(
        stderr, "  zone_xpath = XPath to select zones for flows\n\n"
      );
      fprintf(
        stderr, "  reac_xpath = reaction xpath (optional)\n\n"
      );
      return EXIT_FAILURE;
   }

  /*============================================================================
  // Read file and exit if not present.
  //==========================================================================*/

  p_my_nucnet = Libnucnet__new_from_xml( argv[1], NULL, NULL, argv[2] );

  if( !p_my_nucnet ) {
    fprintf( stderr, "Input data not read!\n" );
    return EXIT_FAILURE;
  }

  //============================================================================
  // Register rate functions.
  //============================================================================

  user::register_my_rate_functions(
    Libnucnet__Net__getReac(
      Libnucnet__getNet( p_my_nucnet )
    )
  );

  register_my_rate_functions(
    Libnucnet__Net__getReac(
      Libnucnet__getNet( p_my_nucnet )
    )
  );
    
  //============================================================================
  // Get the zones.
  //============================================================================

  Libnucnet__setZoneCompareFunction(
    p_my_nucnet,
    (Libnucnet__Zone__compare_function) nnt::zone_compare_by_first_label
  );

  nnt::zone_list_t zone_list = nnt::make_zone_list( p_my_nucnet );

  //============================================================================
  // Get valid net view.
  //============================================================================

  if( argc == 3 )
    p_net_view =
      Libnucnet__NetView__new( Libnucnet__getNet( p_my_nucnet ), "", "" );
  else
    p_net_view =
      Libnucnet__NetView__new( Libnucnet__getNet( p_my_nucnet ), "", argv[3] );

  //============================================================================
  // Set compare function.
  //============================================================================

  Libnucnet__Reac__setReactionCompareFunction(
    Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) ),
    (Libnucnet__Reaction__compare_function)
       nnt::compare_reactions_by_string
  );

  //============================================================================
  // Get the reactions.
  //============================================================================

  nnt::reaction_list_t reaction_list =
    nnt::make_reaction_list(
      Libnucnet__Net__getReac( Libnucnet__NetView__getNet( p_net_view ) )
    );

  Libnucnet__NetView__free( p_net_view );

  //============================================================================
  // Iterate the zones.
  //============================================================================

  BOOST_FOREACH( nnt::Zone zone, zone_list )
  {
   
    user::update_my_rate_functions_data( zone );
    update_my_rate_functions_data( zone );

    //==========================================================================
    // Print conditions.
    //==========================================================================

    if( zone.hasProperty( nnt::s_TIME ) )
      std::cout <<
        "time(s) = " << zone.getProperty( nnt::s_TIME ) << " " <<
        "t9 = " << zone.getProperty( nnt::s_T9 ) << " " <<
        "rho(g/cc) = " << zone.getProperty( nnt::s_RHO ) << " " <<
        std::endl;
    else
      std::cout <<
        "t9 = " << zone.getProperty( nnt::s_T9 ) << " " <<
        "rho(g/cc) = " << zone.getProperty( nnt::s_RHO ) << " " <<
        std::endl;

    fprintf( stdout, "\n\t\t\tReaction\t\t\t  Forward     Reverse     Net\n" );
    printf(
      "=======================================================   =========   =========   =========\n"
    );

    //==========================================================================
    // Print flows.
    //==========================================================================
  
    BOOST_FOREACH( nnt::Reaction reaction, reaction_list )
    {

      flows =
        user::compute_flows_for_reaction(
          zone,
          reaction.getNucnetReaction()
      );

      if( fabs( flows.first - flows.second ) > 1.e-50 )
      {

        fprintf(
      	  stdout,
	  "%-55s%12.8e%12.3e%12.3e\n",
          Libnucnet__Reaction__getString( reaction.getNucnetReaction() ),
          flows.first,
          flows.second,
          flows.first - flows.second
        );

      }

    }

    fprintf( stdout, "\n" );

  }

  //===========================================================================
  // Clean up and exit.
  //==========================================================================*/

  Libnucnet__free( p_my_nucnet );

  return EXIT_SUCCESS;

}
