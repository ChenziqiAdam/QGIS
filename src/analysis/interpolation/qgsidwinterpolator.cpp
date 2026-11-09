/***************************************************************************
                              qgsidwinterpolator.cpp
                              ----------------------
  begin                : Marco 10, 2008
  copyright            : (C) 2008 by Marco Hugentobler
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

#include "qgsidwinterpolator.h"

#include <algorithm>
#include <cmath>
#include <limits>

#include "qgis.h"
#include "qgsgeometryutils_base.h"
#include "qgssanitizer.h"

QgsIDWInterpolator::QgsIDWInterpolator( const QList<LayerData> &layerData )
  : QgsInterpolator( layerData )
{}

int QgsIDWInterpolator::interpolatePoint( double x, double y, double &result, QgsFeedback *feedback )
{
  if ( !mDataIsCached )
  {
    cacheBaseData( feedback );
  }

  double sumCounter = 0;
  double sumDenominator = 0;

  for ( const QgsInterpolatorVertexData &vertex : std::as_const( mCachedBaseData ) )
  {
    double distance = QgsGeometryUtilsBase::distance2D( vertex.x, vertex.y, x, y );
    if ( qgsDoubleNear( distance, 0.0 ) )
    {
      result = vertex.z;
      return 0;
    }
    double currentWeight = 1 / ( std::pow( distance, mDistanceCoefficient ) );
    sumCounter += ( currentWeight * vertex.z );
    sumDenominator += currentWeight;
  }

  if ( sumDenominator == 0.0 )
  {
    return 1;
  }

  result = sumCounter / sumDenominator;

  // shared precondition: coefficient in (0, 10], finite data, query away from every node; returns the independent
  // Shepard reference (long double, weights relative to the nearest vertex) and the data range
  const auto reference = [&]( long double &ref, double &zMin, double &zMax, double &zAbs ) {
    if ( !( mDistanceCoefficient > 0 && mDistanceCoefficient <= 10 ) || mCachedBaseData.isEmpty() || !std::isfinite( x ) || !std::isfinite( y ) )
      return false;
    double dMin = std::numeric_limits<double>::max();
    zMin = std::numeric_limits<double>::max();
    zMax = -zMin;
    zAbs = 0;
    for ( const QgsInterpolatorVertexData &v : std::as_const( mCachedBaseData ) )
    {
      const double d = std::hypot( v.x - x, v.y - y );
      if ( !std::isfinite( d ) || !std::isfinite( v.z ) || !( d > 0 ) )
        return false;
      dMin = std::min( dMin, d );
      zMin = std::min( zMin, v.z );
      zMax = std::max( zMax, v.z );
      zAbs = std::max( zAbs, std::fabs( v.z ) );
    }
    // floating-point range: the implementation sums z / d^p unscaled, which must not overflow double
    if ( !( std::fabs( zAbs ) * std::pow( dMin, -mDistanceCoefficient ) < 1e150 ) )
      return false;
    long double num = 0, den = 0;
    for ( const QgsInterpolatorVertexData &v : std::as_const( mCachedBaseData ) )
    {
      const long double w = std::pow( static_cast<long double>( dMin ) / std::hypot( v.x - x, v.y - y ), static_cast<long double>( mDistanceCoefficient ) );
      num += w * v.z;
      den += w;
    }
    ref = num / den;
    return true;
  };
  QgsSanitizer::run( "QG-INT-002", [&] {
    long double ref;
    double zMin, zMax, zAbs;
    if ( !reference( ref, zMin, zMax, zAbs ) )
      return;
    const double slack = 8 * mCachedBaseData.size() * QgsSanitizer::EPS * zAbs;
    const double v = QgsSanitizer::fault( "QG-INT-002" ) ? zMax + 1.0 + std::fabs( zMax ) : result;
    QgsSanitizer::check( "QG-INT-002", !( v >= zMin - slack && v <= zMax + slack ), "IDW value " + QgsSanitizer::num( v ) + " outside data range [" + QgsSanitizer::num( zMin ) + "," + QgsSanitizer::num( zMax ) + "]" );
  } );
  return 0;
}
