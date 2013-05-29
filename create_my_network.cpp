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

#include "my_molecule_utilities.h"
#include "my_reaction_utilities.h"

int main()
{

  unsigned int i_max = 10000;

  Libnucnet__Net * p_net;

  p_net = Libnucnet__Net__new();

  my_user::add_molecules_to_nuc( 
    Libnucnet__Net__getNuc( p_net ),
    1,
    i_max
  );

  my_user::add_reactions_to_net( p_net, 1, 10 );
  my_user::add_forward_reactions_to_net( p_net, 10, i_max );

  Libnucnet__Net__writeToXmlFile( 
    p_net,
    "data/my_net.xml" 
  );

  Libnucnet__Net__free( p_net );

  return EXIT_SUCCESS;

}
