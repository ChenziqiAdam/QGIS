/***************************************************************************
                          qgsaspectfilter.cpp  -  description
                          -----------------------------------
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

#include "qgsaspectfilter.h"

#include "qgsterrainsanitizer.h"

#include <cmath>

QgsAspectFilter::QgsAspectFilter( const QString &inputFile, const QString &outputFile, const QString &outputFormat )
  : QgsDerivativeFilter( inputFile, outputFile, outputFormat )
{}

float QgsAspectFilter::processNineCellWindow( float *x11, float *x21, float *x31, float *x12, float *x22, float *x32, float *x13, float *x23, float *x33 )
{
  const float derX = calcFirstDerX( *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 );
  const float derY = calcFirstDerY( *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 );

  if ( derX == mOutputNodataValue || derY == mOutputNodataValue || ( derX == 0.0 && derY == 0.0 ) )
  {
    return mOutputNodataValue;
  }
  else
  {
    const float result = 180.0 + std::atan2( derX, derY ) * 180.0 / M_PI;

    QgsSanitizer::run( "QG-TER-002", [&] {
      const float w[9] = { *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 };
      if ( !QgsTerrainSanitizer::noNodata( w, mInputNodataValue ) || !std::isfinite( result ) || !( mCellSizeX > 0 ) || !( mCellSizeY > 0 ) || !std::isfinite( mZFactor ) )
        return;
      double gx, gy;
      QgsTerrainSanitizer::hornGradient( w, mCellSizeX, mCellSizeY, gx, gy );
      const double g = std::hypot( gx, gy ) * std::fabs( mZFactor );
      const double errG = 64 * QgsSanitizer::EPS32 * QgsTerrainSanitizer::maxAbs( w ) * std::fabs( mZFactor ) / std::min( mCellSizeX, mCellSizeY );
      if ( !( g > 1e3 * errG ) )
        return; // aspect undefined or dominated by rounding
      float mirrored[9]; // east-west mirror: m[r][c] = w[r][2 - c]
      for ( int r = 0; r < 3; ++r )
        for ( int c = 0; c < 3; ++c )
          mirrored[r * 3 + c] = w[r * 3 + 2 - c];
      const float r2 = QGS_TERRAIN_EVAL( mirrored );
      if ( static_cast<double>( r2 ) == mOutputNodataValue || !std::isfinite( r2 ) )
        return;
      const double a = QgsSanitizer::seen( "QG-TER-002", result );
      QgsSanitizer::check( "QG-TER-002", QgsTerrainSanitizer::circularDiffDeg( r2, 360.0 - a ) > 1e-3 + QgsTerrainSanitizer::RAD2DEG * errG / g, "aspect " + QgsSanitizer::num( a ) + " but mirrored window gives " + QgsSanitizer::num( r2 ) );
    } );
    QgsSanitizer::run( "QG-TER-003", [&] {
      const float w[9] = { *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 };
      QgsTerrainSanitizer::shiftInvariance( w, mInputNodataValue, mOutputNodataValue, mCellSizeX, mCellSizeY, mZFactor, result, QgsTerrainSanitizer::Aspect, QGS_TERRAIN_EVAL );
    } );
    return result;
  }
}
