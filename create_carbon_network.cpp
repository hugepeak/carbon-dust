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
//! \brief Example code to create carbon network. 
////////////////////////////////////////////////////////////////////////////////

#include "carbon_molecule_utilities.h"
#include "carbon_reaction_utilities.h"

int main()
{

  Libnucnet__Net * p_net;

  p_net = Libnucnet__Net__new();

  add_default_molecules_to_nuc( 
    Libnucnet__Net__getNuc( p_net )
  );


  add_default_carbon_carbon_reactions_to_net( p_net );

  add_default_carbon_oxygen_reactions_to_net( p_net );

  add_default_ion_molecule_reactions_to_net( p_net );

  //add_default_electronic_recombination_reactions_to_net( p_net );

  //add_default_o2_reactions_to_net( p_net );

  add_default_co_reactions_to_net( p_net );

  add_default_single_rate_reactions_to_net( p_net );

  //add_default_carbon_condensation_reactions_to_net( p_net );

  Libnucnet__Net__writeToXmlFile( 
    p_net,
    "data/my_net.xml" 
  );

  Libnucnet__Net__free( p_net );

  return EXIT_SUCCESS;

}
