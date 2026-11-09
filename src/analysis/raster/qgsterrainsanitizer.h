/***************************************************************************
  qgsterrainsanitizer.h
  ---------------------
  Shared helpers for the terrain sanitizers (QG-TER-*). Research
  instrumentation: see qgssanitizer.h. Windows are row-major
  w[0..8] = x11 x21 x31 x12 x22 x32 x13 x23 x33 (rows top to bottom,
  columns left to right; the filters' x_ij are x_{column,row}).
 ***************************************************************************/
#ifndef QGSTERRAINSANITIZER_H
#define QGSTERRAINSANITIZER_H

#include <algorithm>
#include <cmath>
#include <string>

#include "qgssanitizer.h"

//! Evaluator of the enclosing filter's processNineCellWindow on a copy of a row-major window.
#define QGS_TERRAIN_EVAL \
  [this]( const float *w ) { return QgsTerrainSanitizer::callWindow( [this]( float *a, float *b, float *c, float *d, float *e, float *f, float *g, float *h, float *i ) { return this->processNineCellWindow( a, b, c, d, e, f, g, h, i ); }, w ); }

namespace QgsTerrainSanitizer
{
  constexpr double RAD2DEG = 57.29577951308232;

  //! Calls a filter's nine-argument evaluation on a copy of the window.
  template<class F> float callWindow( F &&f, const float *w )
  {
    float c[9];
    std::copy( w, w + 9, c );
    return f( &c[0], &c[1], &c[2], &c[3], &c[4], &c[5], &c[6], &c[7], &c[8] );
  }

  inline float maxAbs( const float *w )
  {
    float m = 0;
    for ( int i = 0; i < 9; ++i )
      m = std::max( m, std::fabs( w[i] ) );
    return m;
  }

  inline bool noNodata( const float *w, double nodata )
  {
    for ( int i = 0; i < 9; ++i )
      if ( !std::isfinite( w[i] ) || static_cast<double>( w[i] ) == nodata )
        return false;
    return true;
  }

  //! Horn gradient (east, north) of a full window in double precision, before the z-factor.
  inline void hornGradient( const float *w, double cellX, double cellY, double &gx, double &gy )
  {
    gx = ( ( double( w[2] ) + 2.0 * w[5] + w[8] ) - ( double( w[0] ) + 2.0 * w[3] + w[6] ) ) / ( 8.0 * cellX );
    gy = ( ( double( w[0] ) + 2.0 * w[1] + w[2] ) - ( double( w[6] ) + 2.0 * w[7] + w[8] ) ) / ( 8.0 * cellY );
  }

  //! circular difference of two angles in degrees, in [0, 180]
  inline double circularDiffDeg( double a, double b )
  {
    double d = std::fmod( a - b, 360.0 );
    if ( d < 0 )
      d += 360.0;
    return std::min( d, 360.0 - d );
  }

  //! kinds for the elevation-shift law
  enum Kind { Slope, Aspect, Hillshade, Ruggedness, Curvature };

  /**
   * QG-TER-003: derivatives depend on elevation differences, not absolute elevation.
   * Re-evaluates the window with every cell shifted by c = max|z|.
   */
  template<class F> void shiftInvariance( const float *w, double inNodata, double outNodata, double cellX, double cellY, double zfactor, float result, Kind kind, F &&eval )
  {
    if ( !noNodata( w, inNodata ) || !std::isfinite( result ) || static_cast<double>( result ) == outNodata || !std::isfinite( zfactor ) )
      return;
    const double m = maxAbs( w );
    const float c = static_cast<float>( std::max( m, 1.0 ) ); // |c| <= 10 max|z| (1 for an all-zero window)
    float shifted[9];
    for ( int i = 0; i < 9; ++i )
    {
      shifted[i] = w[i] + c;
      if ( static_cast<double>( shifted[i] ) == inNodata )
        return; // precondition: no shifted cell hits the nodata value
    }
    const double cell = std::min( cellX, cellY );
    if ( !( cell > 0 ) )
      return;
    const float r2 = eval( shifted );
    if ( !std::isfinite( r2 ) || static_cast<double>( r2 ) == outNodata )
      return;
    const double mag = m + c;
    const double a = QgsSanitizer::seen( "QG-TER-003", result ), b = r2;
    double diff = 0, tol = 0;
    switch ( kind )
    {
      case Slope:
        diff = std::fabs( a - b );
        tol = 1e-4 + 64 * RAD2DEG * QgsSanitizer::EPS32 * mag * std::fabs( zfactor ) / cell + 4 * QgsSanitizer::EPS32 * std::fabs( a );
        break;
      case Aspect:
      {
        double gx, gy;
        hornGradient( w, cellX, cellY, gx, gy );
        const double g = std::hypot( gx, gy ) * std::fabs( zfactor );
        const double errG = 64 * QgsSanitizer::EPS32 * mag * std::fabs( zfactor ) / cell;
        if ( !( g > 1e3 * errG ) )
          return; // direction undefined or dominated by rounding
        diff = circularDiffDeg( a, b );
        tol = 1e-3 + RAD2DEG * errG / g + 4 * QgsSanitizer::EPS32 * std::fabs( a );
        break;
      }
      case Hillshade:
        diff = std::fabs( a - b );
        tol = 255.0 * ( 64 * QgsSanitizer::EPS32 + 8 * QgsSanitizer::EPS32 * mag * std::fabs( zfactor ) / cell ) + 1e-3 + 4 * QgsSanitizer::EPS32 * std::fabs( a );
        break;
      case Ruggedness:
        diff = std::fabs( a - b );
        tol = 64 * QgsSanitizer::EPS32 * mag + 4 * QgsSanitizer::EPS32 * std::fabs( a );
        break;
      case Curvature:
      {
        const double e = 64 * QgsSanitizer::EPS32 * mag / ( cell * cell );
        diff = std::fabs( a - b );
        tol = 4 * e * ( std::sqrt( std::fabs( a ) ) + e ) + 4 * QgsSanitizer::EPS32 * std::fabs( a );
        break;
      }
    }
    QgsSanitizer::check( "QG-TER-003", diff > tol, "terrain value changed from " + QgsSanitizer::num( a ) + " to " + QgsSanitizer::num( b ) + " under a uniform elevation shift of " + QgsSanitizer::num( c ) );
  }
} // namespace QgsTerrainSanitizer

#endif
