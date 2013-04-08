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
//! \brief A header file for the carbon molecule utilities file.
////////////////////////////////////////////////////////////////////////////////

//##############################################################################
// Includes.
//##############################################################################

#include <iostream>
//#include <string>
#include <Libnucnet.h>

//##############################################################################
// Prototypes.
//##############################################################################

void
add_default_molecules_to_nuc(
  Libnucnet__Nuc *
);

void
add_carbon_molecules_to_nuc(
  Libnucnet__Nuc *, unsigned int, unsigned int, unsigned int, unsigned int
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

