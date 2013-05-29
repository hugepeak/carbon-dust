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
//! \brief A header file for the carbon utilities file.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
#include <utility>
#include <vector>
#include <Libnucnet.h>
#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>

namespace my_user
{

struct ReactionData {

  ~ReactionData();
  void clear();
  std::vector<std::pair<unsigned int,unsigned int> > vReactants;
  std::vector<std::pair<unsigned int,unsigned int> > vProducts;
  std::string sSource;
  std::string sUserRateFunctionKey;
  std::vector<std::pair<std::string,std::string> > vUserRateData;

};

//##############################################################################
// Prototypes.
//##############################################################################

void
add_default_reactions_to_net( Libnucnet__Net * );

void
add_reactions_to_net( 
  Libnucnet__Net *, unsigned int, unsigned int
);

void
add_forward_reactions_to_net( 
  Libnucnet__Net *, unsigned int, unsigned int
);

void
add_reverse_reactions_to_net( 
  Libnucnet__Net *, unsigned int, unsigned int
);

void
add_reaction_to_net(
  Libnucnet__Net *,
  ReactionData &
);

}
