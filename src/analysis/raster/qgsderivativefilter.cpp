/***************************************************************************
                          qgsderivativefilter.cpp  -  description
                          -----------------------
    begin                : August 7th, 2009
    copyright            : (C) 2009 by Marco Hugentobler
    email                : marco dot hugentobler at karto dot baug dot ethz dot ch
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgsderivativefilter.h"

#include "qgsterrainsanitizer.h"

QgsDerivativeFilter::QgsDerivativeFilter( const QString &inputFile, const QString &outputFile, const QString &outputFormat )
  : QgsNineCellFilter( inputFile, outputFile, outputFormat )
{}

float QgsDerivativeFilter::calcFirstDerX( float x11, float x21, float x31, float x12, float x22, float x32, float x13, float x23, float x33 ) const
{
  //the basic formula would be simple, but we need to test for nodata values...
  //return (( (*x31 - *x11) + 2 * (*x32 - *x12) + (*x33 - *x13) ) / (8 * mCellSizeX));

  int weight = 0;
  double sum = 0;

  //first row
  if ( x31 != mInputNodataValue && x11 != mInputNodataValue ) //the normal case
  {
    sum += ( x31 - x11 );
    weight += 2;
  }
  else if ( x31 == mInputNodataValue && x11 != mInputNodataValue && x21 != mInputNodataValue ) //probably 3x3 window is at the border
  {
    sum += ( x21 - x11 );
    weight += 1;
  }
  else if ( x11 == mInputNodataValue && x31 != mInputNodataValue && x21 != mInputNodataValue ) //probably 3x3 window is at the border
  {
    sum += ( x31 - x21 );
    weight += 1;
  }

  //second row
  if ( x32 != mInputNodataValue && x12 != mInputNodataValue ) //the normal case
  {
    sum += 2 * ( x32 - x12 );
    weight += 4;
  }
  else if ( x32 == mInputNodataValue && x12 != mInputNodataValue && x22 != mInputNodataValue )
  {
    sum += 2 * ( x22 - x12 );
    weight += 2;
  }
  else if ( x12 == mInputNodataValue && x32 != mInputNodataValue && x22 != mInputNodataValue )
  {
    sum += 2 * ( x32 - x22 );
    weight += 2;
  }

  //third row
  if ( x33 != mInputNodataValue && x13 != mInputNodataValue ) //the normal case
  {
    sum += ( x33 - x13 );
    weight += 2;
  }
  else if ( x33 == mInputNodataValue && x13 != mInputNodataValue && x23 != mInputNodataValue )
  {
    sum += ( x23 - x13 );
    weight += 1;
  }
  else if ( x13 == mInputNodataValue && x33 != mInputNodataValue && x23 != mInputNodataValue )
  {
    sum += ( x33 - x23 );
    weight += 1;
  }

  if ( weight == 0 )
  {
    return mOutputNodataValue;
  }

  const float result = sum / ( weight * mCellSizeX ) * mZFactor;

  QgsSanitizer::run( "QG-TER-005", [&] {
    const float w[9] = { x11, x21, x31, x12, x22, x32, x13, x23, x33 };
    if ( !QgsTerrainSanitizer::noNodata( w, mInputNodataValue ) || !std::isfinite( result ) || !( mCellSizeX > 0 ) || !( mCellSizeY > 0 ) || !std::isfinite( mZFactor ) )
      return;
    // least-squares plane z = z0 + a x + b y on the symmetric 3x3 grid (x east, y north)
    double mean = 0, a = 0, b = 0;
    for ( int r = 0; r < 3; ++r )
      for ( int c = 0; c < 3; ++c )
      {
        mean += w[r * 3 + c] / 9.0;
        a += w[r * 3 + c] * ( c - 1 ) / ( 6.0 * mCellSizeX );
        b += w[r * 3 + c] * ( 1 - r ) / ( 6.0 * mCellSizeY );
      }
    const double m = QgsTerrainSanitizer::maxAbs( w );
    if ( !( std::fabs( mZFactor ) >= 1e-30 ) )
      return; // float32 gradient would be subnormal
    for ( int r = 0; r < 3; ++r )
      for ( int c = 0; c < 3; ++c )
        if ( std::fabs( w[r * 3 + c] - ( mean + a * ( c - 1 ) * mCellSizeX + b * ( 1 - r ) * mCellSizeY ) ) > 8 * QgsSanitizer::EPS32 * m )
          return; // window is not planar
    const double expected = a * mZFactor;
    const double obs = QgsSanitizer::seen( "QG-TER-005", result );
    QgsSanitizer::check( "QG-TER-005", std::fabs( obs - expected ) > 64 * QgsSanitizer::EPS32 * m * std::fabs( mZFactor ) / std::min( mCellSizeX, mCellSizeY ), "gradient X " + QgsSanitizer::num( obs ) + " vs plane " + QgsSanitizer::num( expected ) );
  } );
  return result;
}

float QgsDerivativeFilter::calcFirstDerY( float x11, float x21, float x31, float x12, float x22, float x32, float x13, float x23, float x33 ) const
{
  //the basic formula would be simple, but we need to test for nodata values...
  //return (((*x11 - *x13) + 2 * (*x21 - *x23) + (*x31 - *x33)) / ( 8 * mCellSizeY));

  double sum = 0;
  int weight = 0;

  //first row
  if ( x11 != mInputNodataValue && x13 != mInputNodataValue ) //normal case
  {
    sum += ( x11 - x13 );
    weight += 2;
  }
  else if ( x11 == mInputNodataValue && x13 != mInputNodataValue && x12 != mInputNodataValue )
  {
    sum += ( x12 - x13 );
    weight += 1;
  }
  else if ( x31 == mInputNodataValue && x11 != mInputNodataValue && x12 != mInputNodataValue )
  {
    sum += ( x11 - x12 );
    weight += 1;
  }

  //second row
  if ( x21 != mInputNodataValue && x23 != mInputNodataValue )
  {
    sum += 2 * ( x21 - x23 );
    weight += 4;
  }
  else if ( x21 == mInputNodataValue && x23 != mInputNodataValue && x22 != mInputNodataValue )
  {
    sum += 2 * ( x22 - x23 );
    weight += 2;
  }
  else if ( x23 == mInputNodataValue && x21 != mInputNodataValue && x22 != mInputNodataValue )
  {
    sum += 2 * ( x21 - x22 );
    weight += 2;
  }

  //third row
  if ( x31 != mInputNodataValue && x33 != mInputNodataValue )
  {
    sum += ( x31 - x33 );
    weight += 2;
  }
  else if ( x31 == mInputNodataValue && x33 != mInputNodataValue && x32 != mInputNodataValue )
  {
    sum += ( x32 - x33 );
    weight += 1;
  }
  else if ( x33 == mInputNodataValue && x31 != mInputNodataValue && x32 != mInputNodataValue )
  {
    sum += ( x31 - x32 );
    weight += 1;
  }

  if ( weight == 0 )
  {
    return mOutputNodataValue;
  }

  const float result = sum / ( weight * mCellSizeY ) * mZFactor;

  QgsSanitizer::run( "QG-TER-005", [&] {
    const float w[9] = { x11, x21, x31, x12, x22, x32, x13, x23, x33 };
    if ( !QgsTerrainSanitizer::noNodata( w, mInputNodataValue ) || !std::isfinite( result ) || !( mCellSizeX > 0 ) || !( mCellSizeY > 0 ) || !std::isfinite( mZFactor ) )
      return;
    // least-squares plane z = z0 + a x + b y on the symmetric 3x3 grid (x east, y north)
    double mean = 0, a = 0, b = 0;
    for ( int r = 0; r < 3; ++r )
      for ( int c = 0; c < 3; ++c )
      {
        mean += w[r * 3 + c] / 9.0;
        a += w[r * 3 + c] * ( c - 1 ) / ( 6.0 * mCellSizeX );
        b += w[r * 3 + c] * ( 1 - r ) / ( 6.0 * mCellSizeY );
      }
    const double m = QgsTerrainSanitizer::maxAbs( w );
    if ( !( std::fabs( mZFactor ) >= 1e-30 ) )
      return; // float32 gradient would be subnormal
    for ( int r = 0; r < 3; ++r )
      for ( int c = 0; c < 3; ++c )
        if ( std::fabs( w[r * 3 + c] - ( mean + a * ( c - 1 ) * mCellSizeX + b * ( 1 - r ) * mCellSizeY ) ) > 8 * QgsSanitizer::EPS32 * m )
          return; // window is not planar
    const double expected = b * mZFactor;
    const double obs = QgsSanitizer::seen( "QG-TER-005", result );
    QgsSanitizer::check( "QG-TER-005", std::fabs( obs - expected ) > 64 * QgsSanitizer::EPS32 * m * std::fabs( mZFactor ) / std::min( mCellSizeX, mCellSizeY ), "gradient Y " + QgsSanitizer::num( obs ) + " vs plane " + QgsSanitizer::num( expected ) );
  } );
  return result;
}
