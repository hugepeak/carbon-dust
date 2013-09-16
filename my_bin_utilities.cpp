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
//! \brief Example code for bin utitilies.
////////////////////////////////////////////////////////////////////////////////

#include "my_bin_utilities.h"

//##############################################################################
// update_bin_rates().
//##############################################################################

void
update_bin_rates(
  nnt::Zone & zone,
  WnMatrix * p_matrix,
  gsl_vector * p_rhs
)
{

  if( zone.getProperty( nnt::s_SOLVER ) != nnt::s_ARROW ) {
    std::cerr << "Only work for arrow solver for now" << std::endl;
    exit( EXIT_FAILURE );
  }

  double d_begin, d_end; 
  double d_rate;

  //==========================================================================
  // Get parameters.
  //==========================================================================
    
  size_t i_bin_start =
    boost::lexical_cast<size_t>( zone.getProperty( "network end" ) ); 
  size_t i_bin_size = 
    boost::lexical_cast<size_t>( zone.getProperty( "bin factor" ) ); 
  size_t i_bin_number =
    boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) ); 

  if( i_bin_number == 0 )
    return;

  gsl_vector * p_abunds =
    Libnucnet__Zone__getAbundances(
      zone.getNucnetZone()
    );

  double d_na =
    carbon_compute_Ya( zone ) *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
    GSL_CONST_NUM_AVOGADRO;

  size_t i_size = WnMatrix__get_gsl_vector_size( p_rhs );

  //==========================================================================
  // bin_start + c -> b1 
  //==========================================================================
    
  d_begin = 
    (double) i_bin_start * pow( (double) i_bin_size, 0. );
  d_end = 
    (double) i_bin_start * pow( (double) i_bin_size, 1. );

  d_rate = 
    d_na *
    compute_effective_rate( 
      zone, d_begin, d_end
    );

  //--------------------------------------------------------------------------
  // Matrix elements.
  //--------------------------------------------------------------------------
  
  WnMatrix__assignElement(
    p_matrix,
    i_size - 1,
    i_size - 1,
    d_rate * gsl_vector_get( p_abunds, i_size - 1 )
  );

  WnMatrix__assignElement(
    p_matrix,
    i_size - 1,
    i_size,
    d_rate * gsl_vector_get( p_abunds, i_size - 2 )
  );

  WnMatrix__assignElement(
    p_matrix,
    i_size,
    i_size - 1,
    d_rate * gsl_vector_get( p_abunds, i_size - 1 ) * ( d_end - d_begin )
  );

  WnMatrix__assignElement(
    p_matrix,
    i_size,
    i_size,
    d_rate * gsl_vector_get( p_abunds, i_size - 2 ) * ( d_end - d_begin )
  );

  //--------------------------------------------------------------------------
  // Rhs vector. 
  //--------------------------------------------------------------------------
  
  gsl_vector_set( 
    p_rhs,
    i_size - 2,
    gsl_vector_get( p_rhs, i_size - 2 ) -
      d_rate *
      gsl_vector_get( p_abunds, i_size - 2 ) *
      gsl_vector_get( p_abunds, i_size - 1 )
  );

  gsl_vector_set( 
    p_rhs,
    i_size - 1,
    gsl_vector_get( p_rhs, i_size - 1 ) -
      d_rate *
      gsl_vector_get( p_abunds, i_size - 2 ) *
      gsl_vector_get( p_abunds, i_size - 1 ) * 
      ( d_end - d_begin )
  );

  //==========================================================================
  // Loop from b1 + c -> b2 to bN-1 + c -> bN.
  //==========================================================================
   
  for( size_t i = 1; i < i_bin_number; i++ )
  {

    d_begin = 
      (double) i_bin_start * pow( (double) i_bin_size, (double) i );
    d_end = 
      (double) i_bin_start * pow( (double) i_bin_size, (double) i + 1. );

    d_rate = 
      d_na *
      compute_effective_rate( 
        zone, d_begin, d_end
      );

    WnMatrix__assignElement(
      p_matrix,
      i_size,
      i_size,
      d_rate * 
        boost::lexical_cast<double>(
          zone.getProperty( 
            "bin abundance", boost::lexical_cast<std::string>( i ) 
          )
        ) *
        ( d_end - d_begin )
    );

    gsl_vector_set( 
      p_rhs,
      i_size - 1,
      gsl_vector_get( p_rhs, i_size - 1 ) -
        d_rate *
        boost::lexical_cast<double>(
          zone.getProperty( 
            "bin abundance", boost::lexical_cast<std::string>( i ) 
          )
        ) *
        gsl_vector_get( p_abunds, i_size - 1 ) *
        ( d_end - d_begin )
    );

  }

  //==========================================================================
  // Free. 
  //==========================================================================
  
  gsl_vector_free( p_abunds );

}

//##############################################################################
// compute_effective_rate().
//##############################################################################

double
compute_effective_rate(
  nnt::Zone & zone,
  double d_begin,
  double d_end 
)
{

  if( zone.hasProperty( "k1" ) ) {

    return
      boost::lexical_cast<double>( zone.getProperty( "k1" ) ) /
      (
        3. *
        ( pow( d_end, 1./3. ) - pow( d_begin, 1./3. ) )
      );

  } else {

    return
      compute_k1( zone, 12. ) /    
      (
        3. *
        ( pow( d_end, 1./3. ) - pow( d_begin, 1./3. ) )
      );

  }

}

//##############################################################################
// evolve_bin().
//##############################################################################

void
evolve_bin(
  nnt::Zone & zone,
  std::vector<double> & v_bin_old
)
{

  if( zone.getProperty( nnt::s_SOLVER ) != nnt::s_ARROW ) {
    std::cerr << "Only work for arrow solver for now" << std::endl;
    exit( EXIT_FAILURE );
  }

  double d_begin, d_end, d_rate_last;

  size_t i_bin_start =
    boost::lexical_cast<size_t>( zone.getProperty( "network end" ) ); 
  size_t i_bin_size = 
    boost::lexical_cast<size_t>( zone.getProperty( "bin factor" ) ); 
  size_t i_bin_number =
    boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) ); 

  if( i_bin_number == 0 )
    return;

  gsl_vector * p_abunds =
    Libnucnet__Zone__getAbundances( zone.getNucnetZone() );

  double d_na =
    carbon_compute_Ya( zone ) *
    boost::lexical_cast<double>( zone.getProperty( nnt::s_RHO ) ) *
    GSL_CONST_NUM_AVOGADRO;

  size_t i_size = WnMatrix__get_gsl_vector_size( p_abunds );

  //==========================================================================
  // "Allocate" new abundance vector.
  //==========================================================================
    
  std::vector<double> v_new( i_bin_number );

  //==========================================================================
  // Compute effective rates.
  //==========================================================================
    
  d_begin = (double) i_bin_start * pow( (double) i_bin_size, 0. );
  d_end = (double) i_bin_start * pow( (double) i_bin_size, 1. );

  d_rate_last =
    d_na *
    compute_effective_rate( 
      zone, d_begin, d_end
    );

  std::vector<double> v_rates;

  for( size_t i = 1; i < i_bin_number; i++ ) {

    d_begin = (double) i_bin_start * pow( (double) i_bin_size, (double) i );
    d_end = (double) i_bin_start * pow( (double) i_bin_size, (double) i + 1. );

    v_rates.push_back( 
      d_na *
      compute_effective_rate( 
        zone, d_begin, d_end
      )
    );

  }

  if( i_bin_number == 1 ) {
 
    v_new[0] =
      d_rate_last *
      gsl_vector_get( p_abunds, i_size - 2 ) *
      gsl_vector_get( p_abunds, i_size - 1 ) *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
      +
      v_bin_old[0];

  } else {
    
    v_new[0] =
      (
        d_rate_last *
        gsl_vector_get( p_abunds, i_size - 2 ) *
        gsl_vector_get( p_abunds, i_size - 1 ) *
        boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
        +
        v_bin_old[0]
      )
      /
      (
        1. +
        v_rates[0] *     
        gsl_vector_get( p_abunds, i_size - 1 ) *
        boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
      );
      
    for( size_t i = 1; i < i_bin_number - 1; i++ ) {
  
      v_new[i] = 
        (
          v_rates[i-1] *
          v_new[i-1] *
          gsl_vector_get( p_abunds, i_size - 1 ) *
          boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
          +
          v_bin_old[i]
        )
        /
        (
          1. +
          v_rates[i] *     
          gsl_vector_get( p_abunds, i_size - 1 ) *
          boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
        );
    
    }
  
    v_new[i_bin_number - 1] =
      v_rates[i_bin_number - 2] *
      v_new[i_bin_number - 2] *
      gsl_vector_get( p_abunds, i_size - 1 ) *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_DTIME ) )
      +
      v_bin_old[i_bin_number - 1];

   } // i_bin_number > 1

  //==========================================================================
  // Update bin abundances and abundance changes.
  //==========================================================================
  
  std::vector<double> v_work = get_bin_abundances( zone );

  for( size_t i = 1; i <= i_bin_number; i++ ) 
  {

    zone.updateProperty(
      "bin abundance", 
      boost::lexical_cast<std::string>( i ),
      boost::lexical_cast<std::string>( v_new[i-1] )
    );

    zone.updateProperty(
      "bin abundance change", 
      boost::lexical_cast<std::string>( i ),
      boost::lexical_cast<std::string>( v_new[i-1] - v_work[i-1] )
    );

  }

  //==========================================================================
  // Free. 
  //==========================================================================
  
  gsl_vector_free( p_abunds );

}

//##############################################################################
// check_bin_change().
//##############################################################################

double
check_bin_change(
  nnt::Zone & zone
)
{

  double d_check = 0., d_checkT;

  for(
    size_t i = 1;
    i <= boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) );
    i++
  )
  {

    if(
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance", boost::lexical_cast<std::string>( i ) 
        )
      ) > 1.e-10
    )
    {

      d_checkT = 
        fabs(
          boost::lexical_cast<double>( 
            zone.getProperty( 
              "bin abundance change", 
              boost::lexical_cast<std::string>( i ) 
            )
          ) /   
          boost::lexical_cast<double>( 
            zone.getProperty( 
              "bin abundance", 
              boost::lexical_cast<std::string>( i ) 
            )
          ) 
        ); 

      if( d_checkT > d_check )
        d_check = d_checkT;

    }

  }

  return d_check;

}

//##############################################################################
// initialize_bin().
//##############################################################################

void
initialize_bin(
  nnt::Zone & zone
)
{

  for(
    unsigned i = 1;
    i <= boost::lexical_cast<unsigned>( zone.getProperty( "bin number" ) );
    i++
  )
  {

    zone.updateProperty( 
      "bin abundance", 
      boost::lexical_cast<std::string>( i ), 
      "0." 
    );

    zone.updateProperty( 
      "bin abundance change", 
      boost::lexical_cast<std::string>( i ), 
      "0." 
    );

  }

  return;

}

//##############################################################################
// update_bin_properties().
//##############################################################################

void
update_bin_properties(
  nnt::Zone & zone
)
{

  compute_k1( zone, 12. );

  return;

}

//##############################################################################
// compute_k1().
//##############################################################################

double
compute_k1(
  nnt::Zone & zone,
  double d_A
)
{

  double d_Rc = 1.7e-8;

  double d_result =
    M_PI *
    gsl_pow_2( d_Rc ) *
    sqrt(
      3. *
      GSL_CONST_CGSM_BOLTZMANN *
      boost::lexical_cast<double>( zone.getProperty( nnt::s_T9 ) ) *
      GSL_CONST_NUM_GIGA /
      (
        d_A / GSL_CONST_NUM_AVOGADRO
      )
    );

  zone.updateProperty( "k1", boost::lexical_cast<std::string>( d_result ) );

  return d_result;

}

//############################################################################
// bin_update_timestep()
//############################################################################

void
bin_update_timestep(
  nnt::Zone & zone,
  double & d_dt,
  double d_regt,
  double d_regy,
  double d_ymin
)
{

// need to fill this later.
/*
  if ( d_y > p_extra_data->dYmin && !WnMatrix__value_is_zero( d_dy ) ) {

    d_check =
      p_extra_data->dDt *
      p_extra_data->dRegy *
      fabs( d_y / ( d_dy + d_tiny) );

    if ( d_check < p_extra_data->dDtnew ) {
      p_extra_data->dDtnew = d_check;
    }

  }
*/

}

//############################################################################
// compute_bin_sum()
//############################################################################

double
compute_bin_sum(
  nnt::Zone & zone
)
{

  double d_result = 0.;

  for(
    size_t i = 1;
    i <= boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) );
    i++
  )
  {

    d_result +=
      boost::lexical_cast<double>( zone.getProperty( "network end" ) ) *
      pow(
        boost::lexical_cast<double>( zone.getProperty( "bin factor" ) ),
        (double) i
      ) *
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance", boost::lexical_cast<std::string>( i ) 
        )
      );

  }

  return d_result; 

}

//############################################################################
// print_bin_abundances()
//############################################################################

int
print_bin_abundances(
  nnt::Zone & zone
)
{

  for(
    size_t i = 1;
    i <= boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) );
    i++
  )
  {

    printf( "%20lu%20lu%16.6e%16.6e\n",
      (unsigned long)( 
        boost::lexical_cast<double>( zone.getProperty( "network end" ) ) *
        pow(
          boost::lexical_cast<double>( zone.getProperty( "bin factor" ) ),
          (double) i - 1.
        ) + 1.
      ),
      (unsigned long)(
        boost::lexical_cast<double>( zone.getProperty( "network end" ) ) *
        pow(
          boost::lexical_cast<double>( zone.getProperty( "bin factor" ) ),
          (double) i
        )
      ),
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance", boost::lexical_cast<std::string>( i ) 
        )
      ),
      boost::lexical_cast<double>( 
        zone.getProperty( 
          "bin abundance", boost::lexical_cast<std::string>( i ) 
        )
      ) *
        boost::lexical_cast<double>( zone.getProperty( "network end" ) ) *
        pow(
          boost::lexical_cast<double>( zone.getProperty( "bin factor" ) ),
          (double) i
        )
    );

  }

  return 1;

}

//##############################################################################
// get_bin_abundances().
//##############################################################################

std::vector<double>
get_bin_abundances(
  nnt::Zone & zone
)
{

  std::vector<double> v_abunds;

  for( 
    size_t i = 1; 
    i <= boost::lexical_cast<size_t>( zone.getProperty( "bin number" ) ); 
    i++ 
  ) 
  {

    v_abunds.push_back( 
      boost::lexical_cast<double>(
        zone.getProperty( 
          "bin abundance", 
          boost::lexical_cast<std::string>( i ) 
        )
      )
    );

  }

  return v_abunds;

}

