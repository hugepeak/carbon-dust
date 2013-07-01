////////////////////////////////////////////////////////////////////////////////
// Copyright (c) 2011-2012 Clemson University.
//
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
// You should have received a copy of the GNU General Public License
// along with this software; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307
// USA
//
//////////////////////////////////////////////////////////////////////////////*/

#include "carbon_flow_utilities.h"

//############################################################################
// compute_flows_for_reaction().
//##########################################################################//

std::pair<double,double>
compute_flows_for_reaction(
  nnt::Zone &zone,
  Libnucnet__Reaction *p_reaction
)
{

  double d_forward_rate, d_reverse_rate = 0.;
  double d_f, d_r = 0.;
  size_t i_elements;

  d_forward_rate =
    Libnucnet__Reaction__computeRate(
      p_reaction,
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ),
      Libnucnet__Zone__getDataForUserRateFunction(
        zone.getNucnetZone(),
        Libnucnet__Reaction__getRateFunctionKey( p_reaction )
      )
    );

  //==========================================================================
  // Modify rates for reaction.
  //==========================================================================

  user::modify_rates_for_reaction(
    zone,
    p_reaction,
    d_forward_rate,
    d_reverse_rate
  );

  //==========================================================================
  // Iterate the reactants to get the forward flow.
  //==========================================================================

  nnt::reaction_element_list_t reactant_list =
    nnt::make_reaction_nuclide_reactant_list( p_reaction );

  i_elements = reactant_list.size();

  d_f =
    d_forward_rate *
    pow(
      boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ),
      (double) i_elements - 1.
    ) /
      Libnucnet__Reaction__getDuplicateReactantFactor( p_reaction );

  BOOST_FOREACH( nnt::ReactionElement reactant, reactant_list )
  {

    d_f *=
      Libnucnet__Zone__getSpeciesAbundance(
	zone.getNucnetZone(),
	Libnucnet__Nuc__getSpeciesByName(
	  Libnucnet__Net__getNuc(
            Libnucnet__Zone__getNet( zone.getNucnetZone() )
          ),
	  Libnucnet__Reaction__Element__getName(
            reactant.getNucnetReactionElement()
          )
        )
      );

  }

  //==========================================================================
  // Return the flows.
  //==========================================================================
  
  return std::make_pair( d_f, d_r ); 

}

