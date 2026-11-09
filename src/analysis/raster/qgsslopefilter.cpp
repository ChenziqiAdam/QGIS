/***************************************************************************
                          qgsslopefilter.h  -  description
                          --------------------------------
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

#include "qgsslopefilter.h"

#include "qgsterrainsanitizer.h"

#include <cmath>

QgsSlopeFilter::QgsSlopeFilter( const QString &inputFile, const QString &outputFile, const QString &outputFormat )
  : QgsDerivativeFilter( inputFile, outputFile, outputFormat )
{}

float QgsSlopeFilter::processNineCellWindow( float *x11, float *x21, float *x31, float *x12, float *x22, float *x32, float *x13, float *x23, float *x33 )
{
  const float derX = calcFirstDerX( *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 );
  const float derY = calcFirstDerY( *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 );

  if ( derX == mOutputNodataValue || derY == mOutputNodataValue )
  {
    return mOutputNodataValue;
  }

  const float result = std::atan( std::sqrt( derX * derX + derY * derY ) ) * 180.0 / M_PI;

  QgsSanitizer::run( "QG-TER-001", [&] {
    const float w[9] = { *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 };
    if ( !QgsTerrainSanitizer::noNodata( w, mInputNodataValue ) || !std::isfinite( result ) || mCellSizeX != mCellSizeY || !( mCellSizeX > 0 ) || !std::isfinite( mZFactor ) )
      return;
    float rot[9]; // rotated by 90 degrees: rot[r][c] = w[2 - c][r]
    for ( int r = 0; r < 3; ++r )
      for ( int c = 0; c < 3; ++c )
        rot[r * 3 + c] = w[( 2 - c ) * 3 + r];
    const float r2 = QGS_TERRAIN_EVAL( rot );
    const double a = QgsSanitizer::seen( "QG-TER-001", result );
    const double tol = 1e-4 + 64 * QgsTerrainSanitizer::RAD2DEG * QgsSanitizer::EPS32 * QgsTerrainSanitizer::maxAbs( w ) * std::fabs( mZFactor ) / mCellSizeX;
    QgsSanitizer::check( "QG-TER-001", std::fabs( a - r2 ) > tol, "slope " + QgsSanitizer::num( a ) + " vs rotated window " + QgsSanitizer::num( r2 ) );
  } );
  QgsSanitizer::run( "QG-TER-003", [&] {
    const float w[9] = { *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 };
    QgsTerrainSanitizer::shiftInvariance( w, mInputNodataValue, mOutputNodataValue, mCellSizeX, mCellSizeY, mZFactor, result, QgsTerrainSanitizer::Slope, QGS_TERRAIN_EVAL );
  } );
  return result;
}
