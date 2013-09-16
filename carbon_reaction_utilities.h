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

struct ReactionData {

  ~ReactionData();
  void clear();
  std::vector<std::string> vReactants;
  std::vector<std::string> vProducts;
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
add_default_carbon_carbon_reactions_to_net( Libnucnet__Net * );

void
add_default_carbon_oxygen_reactions_to_net( Libnucnet__Net * );

void
add_default_co_reactions_to_net( Libnucnet__Net * );

void
add_default_oxygen_reactions_to_net( Libnucnet__Net * );

void
add_default_ion_molecule_reactions_to_net( Libnucnet__Net * );

void
add_default_electronic_recombination_reactions_to_net( Libnucnet__Net * );

void
add_arrhenius_rate_to_net( 
  Libnucnet__Net *,
  std::string,
  std::string,
  double,
  double,
  double
);

void
add_arrhenius_rate_to_net( 
  Libnucnet__Net *,
  std::string,
  std::string,
  std::string,
  std::string,
  double,
  double,
  double
);

void
add_arrhenius_rate_to_net( 
  Libnucnet__Net *,
  std::string,
  std::string,
  std::string,
  std::string,
  std::string,
  double,
  double,
  double
);

void
add_arrhenius_inverse_rate_to_net( 
  Libnucnet__Net *,
  std::string,
  std::string,
  std::string,
  std::string,
  double,
  double,
  double,
  double
);

void
add_compton_electron_rate_to_net( 
  Libnucnet__Net *,
  std::string,
  std::string,
  std::string,
  std::string,
  std::string,
  double
);

void
add_default_isomer_reactions_to_net(
  Libnucnet__Net *
);

void
add_isomer_reaction_to_net( 
  Libnucnet__Net *,
  std::string,
  std::string,
  double
);

void
add_reaction_to_net(
  Libnucnet__Net *,
  ReactionData &
);

void
add_default_carbon_condensation_reactions_to_net(
  Libnucnet__Net *
);

void
add_carbon_condensation_reaction_to_net(
  Libnucnet__Net *,
  unsigned int
);

void
add_default_generic_reactions_to_net( Libnucnet__Net * );

void
add_generic_reactions_to_net( 
  Libnucnet__Net *, unsigned int, unsigned int
);

void
add_forward_generic_reactions_to_net( 
  Libnucnet__Net *, unsigned int, unsigned int
);

void
add_reverse_generic_reactions_to_net( 
  Libnucnet__Net *, unsigned int, unsigned int
);

