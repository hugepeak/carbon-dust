//////////////////////////////////////////////////////////////////////////////
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
//! \brief Utility code. 
////////////////////////////////////////////////////////////////////////////////

#include "carbon_reaction_utilities.h"

//##############################################################################
// add_default_reactions_to_net().
//##############################################################################

void
add_default_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  add_carbon_carbon_reactions_to_net( p_net );

  add_carbon_oxygen_reactions_to_net( p_net );

  //============================================================================
  // Add O + O -> O2 + gamma (n + n -> nn + gamma) and inverse. RA3.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "n", 
    "nn", "gamma",
    1.e-19, 0., 0. 
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "nn", "gamma",
    "n", "n", 
    1.e-19, 0., 0., 5.15
  );

  //============================================================================
  // Add C + O2 -> CO + O (h1 + nn -> h2 + n) and inverse. NN56 and NN37.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "nn", 
    "h2", "n",
    2.46e-12, 1.5, -613. 
  );

  add_arrhenius_rate_to_net(
    p_net,
    "h2", "n",
    "h1", "nn", 
    1.e-16, 0., 0. 
  );

  //============================================================================
  // Add C + CO -> C2 + O (h1 + h2 -> he2 + n). NN57.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "h2", 
    "he2", "n",
    5.43e-7, -1.5, 57200.0
  );

  //============================================================================
  // Add CO + gamma -> C + O (h2 + gamma -> n + h1). Inverse of RA4.
  //============================================================================
 
  add_arrhenius_inverse_rate_to_net(
    p_net,
    "h2", "gamma",
    "n", "h1", 
    1.58e-17, 0.3, 1297.4, 11.1
  );

  //============================================================================
  // Add CO + e -> C + O + e (h2 + electron -> h1 + n + electron). Clayton 2012.
  //============================================================================
 
  add_compton_electron_rate_to_net(
    p_net,
    "h2", "electron", 
    "h1", "n", "electron",
    1.e5 
  );

  //============================================================================
  // Add O2 + e -> O + O + e (nn + electron -> n + n + electron). Clayton 2013.
  //============================================================================
 
  add_compton_electron_rate_to_net(
    p_net,
    "nn", "electron", 
    "n", "n", "electron",
    5.e4 
  );

  //============================================================================
  // Add C8c -> C8r (o8C -> o8R). Clayton 2012.
  //============================================================================
 
  add_isomer_rate_to_net(
    p_net,
    "o8C",
    "o8R",
    10.
  );
    
}

//##############################################################################
// add_carbon_carbon_reactions_to_net(). 
//##############################################################################

void
add_carbon_carbon_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  //============================================================================
  // Add C + C -> C2 + gamma (h1 + h1 -> he2 + gamma), and the inverse. C1.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "h1", 
    "he2", "gamma",
    4.36e-18, 0.3, 161.3
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "he2", "gamma",
    "h1", "h1", 
    4.36e-18, 0.3, 161.3, 6.30
  );

  //============================================================================
  // Add C + C2 -> C3 + gamma (h1 + he2 -> li3 + gamma), and the inverse. C2.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "he2", 
    "li3", "gamma",
    1.e-17, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "li3", "gamma",
    "h1", "he2", 
    1.e-17, 0., 0., 7.21
  );

  //============================================================================
  // Add C + C3 -> C4 + gamma (h1 + li3 -> be4 + gamma), and the inverse. C3.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "li3", 
    "be4", "gamma",
    1.e-10, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "be4", "gamma",
    "h1", "li3", 
    1.e-10, 0., 0., 5.06
  );

  //============================================================================
  // Add C + C4 -> C5 + gamma (h1 + be4 -> b5 + gamma), and the inverse. C4.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "be4", 
    "b5", "gamma",
    1.e-13, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "b5", "gamma",
    "h1", "be4", 
    1.e-13, 0., 0., 7.06
  );

  //============================================================================
  // Add C + C5 -> C6 + gamma (h1 + b5 -> c6C + gamma), and the inverse. C6.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "b5", 
    "c6C", "gamma",
    1.e-10, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "c6C", "gamma",
    "h1", "b5", 
    1.e-10, 0., 0., 4.96
  );

  //============================================================================
  // Add C + C6 -> C7 + gamma (h1 + c6C -> n7C + gamma), and the inverse. C7.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "c6C", 
    "n7C", "gamma",
    1.e-13, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "n7C", "gamma",
    "h1", "c6C", 
    1.e-13, 0., 0., 7.22
  );

  //============================================================================
  // Add C + C7 -> C8 + gamma (h1 + n7C -> o8C + gamma), and the inverse. C10.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1", "n7C", 
    "o8C", "gamma",
    1.e-10, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "o8C", "gamma",
    "h1", "n7C", 
    1.e-10, 0., 0., 5.68
  );

}

//##############################################################################
// add_carbon_oxygen_reactions_to_net(). 
//##############################################################################

void
add_carbon_oxygen_reactions_to_net(
  Libnucnet__Net * p_net 
)
{

  ReactionData reaction_data;

  //============================================================================
  // Add O + C -> CO + gamma (n + h1 -> h2 + gamma). RA4.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "h1", 
    "h2", "gamma",
    1.58e-17, 0.3, 1297.4
  );

  //============================================================================
  // Add O + C2 -> CO + C (n + he2 -> h2 + h1). C38.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "he2", 
    "h2", "h1",
    5.99e-10, 0., 0.
  );

  //============================================================================
  // Add O + C3 -> CO + C2 (n + li3 -> h2 + he2). C39.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "li3", 
    "h2", "he2",
    1.e-11, 0.3, 1130.
  );

  //============================================================================
  // Add O + C4 -> CO + C3 (n + be4 -> h2 + li3). C40.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "be4", 
    "h2", "li3",
    3.e-10, 0., 1420.
  );

  //============================================================================
  // Add O + C5 -> CO + C4 (n + b5 -> h2 + be4). C41.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "b5", 
    "h2", "be4",
    1.e-11, -0.3, 1130.
  );

  //============================================================================
  // Add O + C6 -> CO + C5 (n + c6C -> h2 + b5). C42.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "c6C", 
    "h2", "b5",
    3.e-10, 0., 1420.
  );

  //============================================================================
  // Add O + C7 -> CO + C6 (n + n7C -> h2 + c6C). C43.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "n7C", 
    "h2", "c6C",
    1.e-11, -0.3, 1130.
  );

  //============================================================================
  // Add O + C8 -> CO + C7 (n + o8C -> h2 + n7C). C44.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "o8C", 
    "h2", "n7C",
    3.e-10, 0., 1420.
  );

}

//##############################################################################
// add_compton_electron_rate_to_net(). 
//##############################################################################

void
add_compton_electron_rate_to_net(
  Libnucnet__Net * p_net,
  std::string s_first_reactant,
  std::string s_second_reactant,
  std::string s_first_product,
  std::string s_second_product,
  std::string s_third_product,
  double d_tau 
)
{

  ReactionData reaction_data;

  reaction_data.sSource = "Clayton PNews3 (2012)";
  reaction_data.sUserRateFunctionKey = "compton electron rate";

  reaction_data.vReactants.push_back( s_first_reactant );
  reaction_data.vReactants.push_back( s_second_reactant );
  reaction_data.vProducts.push_back( s_first_product );
  reaction_data.vProducts.push_back( s_second_product );
  reaction_data.vProducts.push_back( s_third_product );

  reaction_data.vUserRateData.push_back( 
    std::make_pair( "tau", boost::lexical_cast<std::string>( d_tau) ) 
  ); 

  add_reaction_to_net( p_net, reaction_data );

}
  
//##############################################################################
// add_isomer_rate_to_net().
//##############################################################################

void
add_isomer_rate_to_net(
  Libnucnet__Net * p_net,
  std::string s_reactant,
  std::string s_product,
  double d_tau
)
{

  ReactionData reaction_data;

  reaction_data.sSource = "Clayton suggestion (2012)";
  reaction_data.sUserRateFunctionKey = "isomer rate";

  reaction_data.vReactants.push_back( s_reactant );
  reaction_data.vProducts.push_back( s_product );

  reaction_data.vUserRateData.push_back( 
    std::make_pair( "tau", boost::lexical_cast<std::string>( d_tau) ) 
  ); 

  add_reaction_to_net( p_net, reaction_data );

}

//##############################################################################
// add_arrhenius_rate_to_net().
//##############################################################################

void
add_arrhenius_rate_to_net(
  Libnucnet__Net * p_net,
  std::string s_first_reactant,
  std::string s_second_reactant,
  std::string s_first_product,
  std::string s_second_product,
  double d_A,
  double d_nu,
  double d_Ea
)
{

  ReactionData reaction_data;

  reaction_data.sSource = "Cherchneff et al. (2009)";
  reaction_data.sUserRateFunctionKey = "arrhenius rate";

  reaction_data.vReactants.push_back( s_first_reactant );
  reaction_data.vReactants.push_back( s_second_reactant );
  reaction_data.vProducts.push_back( s_first_product );
  reaction_data.vProducts.push_back( s_second_product );

  reaction_data.vUserRateData.push_back( 
    std::make_pair( "A", boost::lexical_cast<std::string>( d_A ) ) 
  ); 
  reaction_data.vUserRateData.push_back( 
    std::make_pair( "nu", boost::lexical_cast<std::string>( d_nu ) ) 
  ); 
  reaction_data.vUserRateData.push_back( 
    std::make_pair( "Ea", boost::lexical_cast<std::string>( d_Ea ) ) 
  ); 

  add_reaction_to_net( p_net, reaction_data );

}

//##############################################################################
// add_arrhenius_inverse_rate_to_net().
//##############################################################################

void
add_arrhenius_inverse_rate_to_net(
  Libnucnet__Net * p_net,
  std::string s_first_reactant,
  std::string s_second_reactant,
  std::string s_first_product,
  std::string s_second_product,
  double d_A,
  double d_nu,
  double d_Ea,
  double d_BondEnergy
)
{

  ReactionData reaction_data;

  reaction_data.sSource = "Clayton et al. (2001) and Cherchneff et al. (2009)";
  reaction_data.sUserRateFunctionKey = "arrhenius inverse rate";

  reaction_data.vReactants.push_back( s_first_reactant );
  reaction_data.vReactants.push_back( s_second_reactant );
  reaction_data.vProducts.push_back( s_first_product );
  reaction_data.vProducts.push_back( s_second_product );

  reaction_data.vUserRateData.push_back( 
    std::make_pair( "A", boost::lexical_cast<std::string>( d_A ) ) 
  ); 
  reaction_data.vUserRateData.push_back( 
    std::make_pair( "nu", boost::lexical_cast<std::string>( d_nu ) ) 
  ); 
  reaction_data.vUserRateData.push_back( 
    std::make_pair( "Ea", boost::lexical_cast<std::string>( d_Ea ) ) 
  ); 
  reaction_data.vUserRateData.push_back( 
    std::make_pair( 
      "BondEnergy", 
      boost::lexical_cast<std::string>( d_BondEnergy ) 
    ) 
  ); 

  add_reaction_to_net( p_net, reaction_data );

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

  BOOST_FOREACH( std::string s_element, reaction_data.vReactants ) {

    if( 
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc( p_net ),
        s_element.c_str()
      ) ||
      s_element == "gamma" || 
      s_element == "electron"
    ) 
      Libnucnet__Reaction__addReactant( 
        p_reaction, 
        s_element.c_str()
      );
    else
    {
      std::cerr << "No " << s_element << " in network!" << std::endl;
      exit( EXIT_FAILURE );
    }
  
  }

  BOOST_FOREACH( std::string s_element, reaction_data.vProducts ) {
 
    if( 
      Libnucnet__Nuc__getSpeciesByName(
        Libnucnet__Net__getNuc( p_net ),
        s_element.c_str()
      ) ||
      s_element == "gamma" || 
      s_element == "electron"
    ) 
      Libnucnet__Reaction__addProduct( 
        p_reaction, 
        s_element.c_str()
      );
    else
    {
      std::cerr << "No " << s_element << " in network!" << std::endl;
      exit( EXIT_FAILURE );
    }
  
  }

  Libnucnet__Reaction__setUserRateFunctionKey(
    p_reaction,
    reaction_data.sUserRateFunctionKey.c_str()
  );

  Libnucnet__Reaction__updateSource( 
    p_reaction, 
    reaction_data.sSource.c_str() 
  );

  std::pair<std::string,std::string> my_data;

  BOOST_FOREACH( my_data, reaction_data.vUserRateData) {
 
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

  this->vReactants.clear();
  this->vProducts.clear();
  this->sSource = "";
  this->sUserRateFunctionKey = "";
  this->vUserRateData.clear();

}
  
