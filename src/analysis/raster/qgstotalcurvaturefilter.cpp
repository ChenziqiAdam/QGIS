/***************************************************************************
                          qgstotalcurvaturefilter.h  -  description
                             -------------------
    begin                : August 21th, 2009
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

#include "qgstotalcurvaturefilter.h"

#include "qgsterrainsanitizer.h"

QgsTotalCurvatureFilter::QgsTotalCurvatureFilter( const QString &inputFile, const QString &outputFile, const QString &outputFormat )
  : QgsNineCellFilter( inputFile, outputFile, outputFormat )
{}

float QgsTotalCurvatureFilter::processNineCellWindow( float *x11, float *x21, float *x31, float *x12, float *x22, float *x32, float *x13, float *x23, float *x33 )
{
  //return nodata if one value is the nodata value
  if ( *x11 == mInputNodataValue
       || *x21 == mInputNodataValue
       || *x31 == mInputNodataValue
       || *x12 == mInputNodataValue
       || *x22 == mInputNodataValue
       || *x32 == mInputNodataValue
       || *x13 == mInputNodataValue
       || *x23 == mInputNodataValue
       || *x33 == mInputNodataValue )
  {
    return mOutputNodataValue;
  }

  const double cellSizeAvg = ( mCellSizeX + mCellSizeY ) / 2.0;
  const double dxx = ( *x32 - 2 * *x22 + *x12 ) / ( mCellSizeX * mCellSizeX );
  const double dxy = ( -*x11 + *x31 + *x13 - *x33 ) / ( 4 * cellSizeAvg * cellSizeAvg );
  const double dyy = ( *x21 - 2 * *x22 + *x23 ) / ( mCellSizeY * mCellSizeY );

  const float result = dxx * dxx + 2 * dxy * dxy + dyy * dyy;

  QgsSanitizer::run( "QG-TER-003", [&] {
    const float w[9] = { *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 };
    QgsTerrainSanitizer::shiftInvariance( w, mInputNodataValue, mOutputNodataValue, mCellSizeX, mCellSizeY, mZFactor, result, QgsTerrainSanitizer::Curvature, QGS_TERRAIN_EVAL );
  } );
  QgsSanitizer::run( "QG-TER-006", [&] {
    const float w[9] = { *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 };
    if ( !QgsTerrainSanitizer::noNodata( w, mInputNodataValue ) || !std::isfinite( result ) || static_cast<double>( result ) == mOutputNodataValue )
      return;
    float doubled[9];
    for ( int i = 0; i < 9; ++i )
    {
      doubled[i] = 2.0f * w[i];
      if ( !std::isfinite( doubled[i] ) || static_cast<double>( doubled[i] ) == mInputNodataValue )
        return;
    }
    const float r2 = QGS_TERRAIN_EVAL( doubled );
    if ( !std::isfinite( r2 ) || static_cast<double>( r2 ) == mOutputNodataValue )
      return;
    const double a = QgsSanitizer::seen( "QG-TER-006", result );
    QgsSanitizer::check( "QG-TER-006", std::fabs( r2 - 4.0 * a ) > 16 * QgsSanitizer::EPS32 * std::fabs( a ), "curvature " + QgsSanitizer::num( a ) + " but doubled elevations give " + QgsSanitizer::num( r2 ) );
  } );
  return result;
}
