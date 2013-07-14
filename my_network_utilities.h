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
//! \brief A header file for my network utilities file.
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
#include "nnt/auxiliary.h"

namespace my_user
{

#define S_PHOTON_END               "photon end"
#define S_RUN_BIN                  "run bin"
#define S_BIN_START                "bin start"
#define S_BIN_SIZE                 "bin size"
#define S_BIN_NUMBER               "bin number"

//##############################################################################
// Prototypes.
//##############################################################################

void
add_default_molecules_to_nuc(
  Libnucnet__Nuc *
);

void
add_molecules_to_nuc(
  Libnucnet__Nuc *, unsigned int, unsigned int
);

void
add_molecule_to_nuc(
  Libnucnet__Nuc *, 
  unsigned int, 
  unsigned int, 
  std::string, 
  gsl_vector *,
  gsl_vector *
);

struct ReactionData {

  ~ReactionData();
  void clear();
  std::vector<std::pair<unsigned int,unsigned int> > vReactants;
  std::vector<std::pair<unsigned int,unsigned int> > vProducts;
  std::string sSource;
  std::string sUserRateFunctionKey;
  std::vector<std::pair<std::string,std::string> > vUserRateData;

};

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

int
update_bin_net(
  Libnucnet *
);

} // my_user
