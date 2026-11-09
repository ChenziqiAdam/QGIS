/***************************************************************************
                          qgshillshadefilter.h  -  description
                          --------------------------------
    begin                : September 26th, 2011
    copyright            : (C) 2011 by Marco Hugentobler
    email                : marco dot hugentobler at sourcepole dot ch
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgshillshadefilter.h"

#include "qgsterrainsanitizer.h"

#include <cmath>

QgsHillshadeFilter::QgsHillshadeFilter( const QString &inputFile, const QString &outputFile, const QString &outputFormat, double lightAzimuth, double lightAngle )
  : QgsDerivativeFilter( inputFile, outputFile, outputFormat )
  , mLightAzimuth( static_cast<float>( lightAzimuth ) )
  , mLightAngle( static_cast<float>( lightAngle ) )
  , mCosZenithRad( std::cos( static_cast<float>( lightAngle * M_PI ) / 180.0f ) )
  , mSinZenithRad( std::sin( static_cast<float>( lightAngle * M_PI ) / 180.0f ) )
  , mAzimuthRad( static_cast<float>( lightAzimuth * M_PI ) / 180.0f )
{}

float QgsHillshadeFilter::processNineCellWindow( float *x11, float *x21, float *x31, float *x12, float *x22, float *x32, float *x13, float *x23, float *x33 )
{
  const float derX = calcFirstDerX( *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 );
  const float derY = calcFirstDerY( *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 );

  if ( derX == mOutputNodataValue || derY == mOutputNodataValue )
  {
    return mOutputNodataValue;
  }

  const float slope_rad = std::atan( std::sqrt( derX * derX + derY * derY ) );
  float aspect_rad = 0;
  if ( derX == 0 && derY == 0 ) //aspect undefined, take a neutral value. Better solutions?
  {
    aspect_rad = mAzimuthRad / 2.0f;
  }
  else
  {
    aspect_rad = M_PI + std::atan2( derX, derY );
  }
  const float result = std::max( 0.0f, 255.0f * ( ( mCosZenithRad * std::cos( slope_rad ) ) + ( mSinZenithRad * std::sin( slope_rad ) * std::cos( mAzimuthRad - aspect_rad ) ) ) );

  QgsSanitizer::run( "QG-TER-003", [&] {
    const float w[9] = { *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 };
    QgsTerrainSanitizer::shiftInvariance( w, mInputNodataValue, mOutputNodataValue, mCellSizeX, mCellSizeY, mZFactor, result, QgsTerrainSanitizer::Hillshade, QGS_TERRAIN_EVAL );
  } );
  QgsSanitizer::run( "QG-TER-004", [&] {
    // Lambert's cosine law: illumination = n . s for the surface normal n (independent Horn gradient)
    // and the unit sun vector s (azimuth clockwise from north; the light angle is the vertical angle of the
    // sun above the horizon, as in QgsHillshadeRenderer and the GDAL convention, so zenith = 90 - angle)
    const float w[9] = { *x11, *x21, *x31, *x12, *x22, *x32, *x13, *x23, *x33 };
    if ( !QgsTerrainSanitizer::noNodata( w, mInputNodataValue ) || !std::isfinite( result ) || !( mCellSizeX > 0 ) || !( mCellSizeY > 0 ) || !std::isfinite( mZFactor ) || !( mLightAngle > 0 && mLightAngle <= 90 ) || !std::isfinite( mLightAzimuth ) )
      return;
    double gx, gy;
    QgsTerrainSanitizer::hornGradient( w, mCellSizeX, mCellSizeY, gx, gy );
    gx *= mZFactor;
    gy *= mZFactor;
    const double zen = ( 90.0 - static_cast<double>( mLightAngle ) ) * M_PI / 180.0, az = static_cast<double>( mLightAzimuth ) * M_PI / 180.0;
    const double ns = ( std::cos( zen ) - std::sin( zen ) * ( gx * std::sin( az ) + gy * std::cos( az ) ) ) / std::sqrt( 1.0 + gx * gx + gy * gy );
    const double expected = 255.0 * std::max( 0.0, ns );
    const double obs = QgsSanitizer::seen( "QG-TER-004", result );
    const double tol = 255.0 * ( 64 * QgsSanitizer::EPS32 + 8 * QgsSanitizer::EPS32 * QgsTerrainSanitizer::maxAbs( w ) * std::fabs( mZFactor ) / std::min( mCellSizeX, mCellSizeY ) ) + 1e-3;
    QgsSanitizer::check( "QG-TER-004", std::fabs( obs - expected ) > tol, "hillshade " + QgsSanitizer::num( obs ) + " vs Lambert " + QgsSanitizer::num( expected ) );
  } );
  return result;
}

void QgsHillshadeFilter::setLightAzimuth( float azimuth )
{
  mLightAzimuth = azimuth;
  mAzimuthRad = azimuth * static_cast<float>( M_PI ) / 180.0f;
}

void QgsHillshadeFilter::setLightAngle( float angle )
{
  mLightAngle = angle;
  mCosZenithRad = std::cos( angle * static_cast<float>( M_PI ) / 180.0f );
  mSinZenithRad = std::sin( angle * static_cast<float>( M_PI ) / 180.0f );
}

#ifdef HAVE_OPENCL

void QgsHillshadeFilter::addExtraRasterParams( std::vector<float> &params )
{
  params.push_back( mCosZenithRad ); // cos_zenith_rad 5
  params.push_back( mSinZenithRad ); // sin_zenith_rad 6
  params.push_back( mAzimuthRad );   // azimuth_rad 7
}

#endif
