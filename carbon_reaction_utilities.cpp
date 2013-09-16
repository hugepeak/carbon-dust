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

  add_default_carbon_carbon_reactions_to_net( p_net );

  add_default_carbon_oxygen_reactions_to_net( p_net );

  add_default_ion_molecule_reactions_to_net( p_net );

  add_default_electronic_recombination_reactions_to_net( p_net );

  add_default_oxygen_reactions_to_net( p_net );

  add_default_co_reactions_to_net( p_net );

  add_default_isomer_reactions_to_net( p_net );

  add_default_carbon_condensation_reactions_to_net( p_net );

}

//##############################################################################
// add_default_oxygen_reactions_to_net(). 
//##############################################################################

void
add_default_oxygen_reactions_to_net(
  Libnucnet__Net * p_net
)
{

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
  // Add O2 + e -> O + O + e (nn + electron -> n + n + electron). Clayton 2013.
  //============================================================================
 
  add_compton_electron_rate_to_net(
    p_net,
    "nn", "electron", 
    "n", "n", "electron",
    5.e4 
  );

}

//##############################################################################
// add_default_co_reactions_to_net(). 
//##############################################################################

void
add_default_co_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  //============================================================================
  // Add C + O2 -> CO + O (h1g + nn -> h2 + n) and inverse. NN56 and NN37.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "nn", 
    "h2", "n",
    2.46e-12, 1.5, -613. 
  );

  add_arrhenius_rate_to_net(
    p_net,
    "h2", "n",
    "h1g", "nn", 
    1.e-16, 0., 0. 
  );

  //============================================================================
  // Add C + CO -> C2 + O (h1g + h2 -> he2c + n). NN57.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "h2", 
    "he2c", "n",
    5.43e-7, -1.5, 57200.0
  );

  //============================================================================
  // Add CO + gamma -> C + O (h2 + gamma -> n + h1g). Inverse of RA4.
  //============================================================================
 
  add_arrhenius_inverse_rate_to_net(
    p_net,
    "h2", "gamma",
    "n", "h1g", 
    1.58e-17, 0.3, 1297.4, 11.1
  );

  //============================================================================
  // Add CO + e -> C + O + e (h2 + electron -> h1g + n + electron). 
  //   Clayton 2012.
  //============================================================================
 
  add_compton_electron_rate_to_net(
    p_net,
    "h2", "electron", 
    "h1g", "n", "electron",
    1.e5 
  );

}

//##############################################################################
// add_default_isomer_reactions_to_net(). 
//##############################################################################

void
add_default_isomer_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  //============================================================================
  // Add C8c -> C8r (o8c -> o8r). Clayton 2012.
  //============================================================================
 
  add_isomer_reaction_to_net(
    p_net,
    "o8c",
    "o8r",
    10.
  );

}
    
//##############################################################################
// add_isomer_reaction_to_net().
//##############################################################################

void
add_isomer_reaction_to_net(
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
// add_default_carbon_condensation_reactions_to_net(). 
//##############################################################################

void
add_default_carbon_condensation_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  //============================================================================
  // Add carbon condensation reactions up to 100.
  //============================================================================
 
  unsigned i_begin, i_end;

  i_begin = 8;
  i_end = 100;

  for( unsigned int i = i_begin; i < i_end; i++ )
    add_carbon_condensation_reaction_to_net( p_net, i );

}

//##############################################################################
// add_carbon_condensation_reaction_to_net(). 
//##############################################################################

void
add_carbon_condensation_reaction_to_net(
  Libnucnet__Net * p_net,
  unsigned int i_start
)
{

  ReactionData reaction_data;

  reaction_data.sSource = "Meyer 2013";
  reaction_data.sUserRateFunctionKey = "carbon condensation rate";

  reaction_data.vReactants.push_back( 
    boost::lexical_cast<std::string>(
      Libnucnet__Species__getName(
        Libnucnet__Nuc__getSpeciesByZA(
          Libnucnet__Net__getNuc( p_net ),
          i_start,
          i_start,
          "r"
        )
      )
    )
  );

  reaction_data.vReactants.push_back( "h1g" );

  reaction_data.vProducts.push_back(
    boost::lexical_cast<std::string>(
      Libnucnet__Species__getName(
        Libnucnet__Nuc__getSpeciesByZA(
          Libnucnet__Net__getNuc( p_net ),
          i_start + 1,
          i_start + 1,
          "r"
        )
      )
    )
  );

  add_reaction_to_net( p_net, reaction_data );

}
 
//##############################################################################
// add_default_carbon_carbon_reactions_to_net(). 
//##############################################################################

void
add_default_carbon_carbon_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  //============================================================================
  // Add C + C -> C2 + gamma (h1g + h1g -> he2c + gamma), and the inverse. C1.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "h1g", 
    "he2c", "gamma",
    4.36e-18, 0.3, 161.3
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "he2c", "gamma",
    "h1g", "h1g", 
    4.36e-18, 0.3, 161.3, 6.30
  );

  //============================================================================
  // Add C + C2 -> C3 + gamma (h1g + he2c -> li3c + gamma), and the inverse. C2.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "he2c", 
    "li3c", "gamma",
    1.e-17, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "li3c", "gamma",
    "h1g", "he2c", 
    1.e-17, 0., 0., 7.21
  );

  //============================================================================
  // Add C + C3 -> C4 + gamma (h1g + li3c -> be4c + gamma), and the inverse. C3.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "li3c", 
    "be4c", "gamma",
    1.e-10, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "be4c", "gamma",
    "h1g", "li3c", 
    1.e-10, 0., 0., 5.06
  );

  //============================================================================
  // Add C + C4 -> C5 + gamma (h1g + be4c -> b5c + gamma), and the inverse. C4.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "be4c", 
    "b5c", "gamma",
    1.e-13, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "b5c", "gamma",
    "h1g", "be4c", 
    1.e-13, 0., 0., 7.06
  );

  //============================================================================
  // Add C + C5 -> C6 + gamma (h1g + b5c -> c6c + gamma), and the inverse. C6.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "b5c", 
    "c6c", "gamma",
    1.e-10, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "c6c", "gamma",
    "h1g", "b5c", 
    1.e-10, 0., 0., 4.96
  );

  //============================================================================
  // Add C + C6 -> C7 + gamma (h1g + c6c -> n7c + gamma), and the inverse. C7.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "c6c", 
    "n7c", "gamma",
    1.e-13, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "n7c", "gamma",
    "h1g", "c6c", 
    1.e-13, 0., 0., 7.22
  );

  //============================================================================
  // Add C + C7 -> C8 + gamma (h1g + n7c -> o8c + gamma), and the inverse. C10.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1g", "n7c", 
    "o8c", "gamma",
    1.e-10, 0., 0.
  );

  add_arrhenius_inverse_rate_to_net(
    p_net,
    "o8c", "gamma",
    "h1g", "n7c", 
    1.e-10, 0., 0., 5.68
  );

}

//##############################################################################
// add_default_carbon_oxygen_reactions_to_net(). 
//##############################################################################

void
add_default_carbon_oxygen_reactions_to_net(
  Libnucnet__Net * p_net 
)
{

  //============================================================================
  // Add O + C -> CO + gamma (n + h1g -> h2 + gamma). RA4.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "h1g", 
    "h2", "gamma",
    1.58e-17, 0.3, 1297.4
  );

  //============================================================================
  // Add O + C2 -> CO + C (n + he2c -> h2 + h1g). C38.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "he2c", 
    "h2", "h1g",
    5.99e-10, 0., 0.
  );

  //============================================================================
  // Add O + C3 -> CO + C2 (n + li3c -> h2 + he2c). C39.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "li3c", 
    "h2", "he2c",
    1.e-11, 0.3, 1130.
  );

  //============================================================================
  // Add O + C4 -> CO + C3 (n + be4c -> h2 + li3c). C40.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "be4c", 
    "h2", "li3c",
    3.e-10, 0., 1420.
  );

  //============================================================================
  // Add O + C5 -> CO + C4 (n + b5c -> h2 + be4c). C41.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "b5c", 
    "h2", "be4c",
    1.e-11, -0.3, 1130.
  );

  //============================================================================
  // Add O + C6 -> CO + C5 (n + c6c -> h2 + b5c). C42.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "c6c", 
    "h2", "b5c",
    3.e-10, 0., 1420.
  );

  //============================================================================
  // Add O + C7 -> CO + C6 (n + n7c -> h2 + c6c). C43.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "n7c", 
    "h2", "c6c",
    1.e-11, -0.3, 1130.
  );

  //============================================================================
  // Add O + C8 -> CO + C7 (n + o8c -> h2 + n7c). C44.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "n", "o8c", 
    "h2", "n7c",
    3.e-10, 0., 1420.
  );

}

//##############################################################################
// add_default_ion_molecule_reactions_to_net(). 
//##############################################################################

void
add_default_ion_molecule_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  //============================================================================
  // Add He+ + CO -> C+ + O + He (he4+ + h2 -> h1+ + n + he4g). IM13.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "he4+", "h2", 
    "h1+", "n", "he4g",
    3.e-10, 0., 1420.
  );

}

//##############################################################################
// add_default_electronic_recombination_reactions_to_net(). 
//##############################################################################

void
add_default_electronic_recombination_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  //============================================================================
  // Add C+ + electron -> C (h1+ + electron -> h1g). ER5.
  //============================================================================
 
  add_arrhenius_rate_to_net(
    p_net,
    "h1+",
    "h1g",
    4.67e-12, -0.6, 0.
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
// add_arrhenius_rate_to_net().
//##############################################################################

void
add_arrhenius_rate_to_net(
  Libnucnet__Net * p_net,
  std::string s_first_reactant,
  std::string s_first_product,
  double d_A,
  double d_nu,
  double d_Ea
)
{

  ReactionData reaction_data;

  reaction_data.sSource = "Cherchneff et al. (2009)";
  reaction_data.sUserRateFunctionKey = "arrhenius rate";

  reaction_data.vReactants.push_back( s_first_reactant );
  reaction_data.vProducts.push_back( s_first_product );

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

void
add_arrhenius_rate_to_net(
  Libnucnet__Net * p_net,
  std::string s_first_reactant,
  std::string s_second_reactant,
  std::string s_first_product,
  std::string s_second_product,
  std::string s_third_product,
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
  reaction_data.vProducts.push_back( s_third_product );

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
// add_default_generic_reactions_to_net().
//##############################################################################

void
add_default_generic_reactions_to_net(
  Libnucnet__Net * p_net
)
{

  add_generic_reactions_to_net(
    p_net,
    1,
    (unsigned int)
      Libnucnet__Nuc__getNumberOfSpecies(
        Libnucnet__Net__getNuc( p_net )
      )
  );

}

//##############################################################################
// add_generic_reactions_to_net().
//##############################################################################

void
add_generic_reactions_to_net(
  Libnucnet__Net * p_net,
  unsigned int i_lo,
  unsigned int i_hi
)
{ 
    
  add_forward_generic_reactions_to_net( p_net, i_lo, i_hi );
  add_reverse_generic_reactions_to_net( p_net, i_lo, i_hi );

}

//##############################################################################
// add_forward_generic_reactions_to_net().
//##############################################################################

void
add_forward_generic_reactions_to_net(
  Libnucnet__Net * p_net,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  ReactionData reaction_data;

  std::string s_element;

  for( unsigned int i = i_lo; i < i_hi; i++ )
  {

    reaction_data.sSource = "Meyer notes 2013";
    reaction_data.sUserRateFunctionKey = "carbon condensation rate";

    reaction_data.vReactants.push_back( "h1" );

    if( 
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_net ),
        i,
        i,
        NULL 
      )
    ) 
      s_element =
        Libnucnet__Species__getName(
          Libnucnet__Nuc__getSpeciesByZA(
            Libnucnet__Net__getNuc( p_net ),
            i,
            i,
            NULL 
          )
        );
    else
    {
      std::cerr << "No number " << i << " atom in network!" << std::endl;
      exit( EXIT_FAILURE );
    }
  
    reaction_data.vReactants.push_back( s_element );

    if( 
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_net ),
        i + 1,
        i + 1,
        NULL 
      )
    ) 
      s_element =
        Libnucnet__Species__getName(
          Libnucnet__Nuc__getSpeciesByZA(
            Libnucnet__Net__getNuc( p_net ),
            i + 1,
            i + 1,
            NULL 
          )
        );
    else
    {
      std::cerr << "No number " << i << " atom in network!" << std::endl;
      exit( EXIT_FAILURE );
    }
  
    reaction_data.vProducts.push_back( s_element );

    add_reaction_to_net( p_net, reaction_data );

    reaction_data.clear();

  }

}
    
//##############################################################################
// add_reverse_generic_reactions_to_net().
//##############################################################################

void
add_reverse_generic_reactions_to_net(
  Libnucnet__Net * p_net,
  unsigned int i_lo,
  unsigned int i_hi
)
{

  ReactionData reaction_data;

  std::string s_element;

  for( unsigned int i = i_lo; i < i_hi; i++ )
  {

    reaction_data.sSource = "Meyer notes 2013";
    reaction_data.sUserRateFunctionKey = "carbon condensation inverse rate";
    reaction_data.vUserRateData.push_back(
      std::make_pair( "BondEnergy", "5.0" )
    );

    if( 
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_net ),
        i + 1,
        i + 1,
        NULL 
      )
    ) 
      s_element =
        Libnucnet__Species__getName(
          Libnucnet__Nuc__getSpeciesByZA(
            Libnucnet__Net__getNuc( p_net ),
            i + 1,
            i + 1,
            NULL 
          )
        );
    else
    {
      std::cerr << "No number " << i << " atom in network!" << std::endl;
      exit( EXIT_FAILURE );
    }
  
    reaction_data.vReactants.push_back( s_element );

    if( 
      Libnucnet__Nuc__getSpeciesByZA(
        Libnucnet__Net__getNuc( p_net ),
        i,
        i,
        NULL 
      )
    ) 
      s_element =
        Libnucnet__Species__getName(
          Libnucnet__Nuc__getSpeciesByZA(
            Libnucnet__Net__getNuc( p_net ),
            i,
            i,
            NULL 
          )
        );
    else
    {
      std::cerr << "No number " << i << " atom in network!" << std::endl;
      exit( EXIT_FAILURE );
    }
  
    reaction_data.vProducts.push_back( s_element );

    reaction_data.vProducts.push_back( "h1" );

    add_reaction_to_net( p_net, reaction_data );

    reaction_data.clear();
    
  }

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
 
