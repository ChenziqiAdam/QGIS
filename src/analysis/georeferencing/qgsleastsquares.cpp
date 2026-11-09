/***************************************************************************
     qgsleastsquares.cpp
     --------------------------------------
    Date                 : Sun Sep 16 12:03:37 AKDT 2007
    Copyright            : (C) 2007 by Gary E. Sherman
    Email                : sherman at mrcc dot com
 ***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "qgsconfig.h"
#include "qgsleastsquares.h"

#include <cmath>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "qgsexception.h"
#include "qgssanitizer.h"

#include <QObject>

#ifdef HAVE_GSL
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#endif

static void linearImpl( const QVector<QgsPointXY> &sourceCoordinates, const QVector<QgsPointXY> &destinationCoordinates, QgsPointXY &origin, double &pixelXSize, double &pixelYSize )
{
  const int n = destinationCoordinates.size();
  if ( n < 2 )
  {
    throw std::domain_error( QObject::tr( "Fit to a linear transform requires at least 2 points." ).toLocal8Bit().constData() );
  }

  double sumPx( 0 ), sumPy( 0 ), sumPx2( 0 ), sumPy2( 0 ), sumPxMx( 0 ), sumPyMy( 0 ), sumMx( 0 ), sumMy( 0 );
  for ( int i = 0; i < n; ++i )
  {
    sumPx += sourceCoordinates.at( i ).x();
    sumPy += sourceCoordinates.at( i ).y();
    sumPx2 += std::pow( sourceCoordinates.at( i ).x(), 2 );
    sumPy2 += std::pow( sourceCoordinates.at( i ).y(), 2 );
    sumPxMx += sourceCoordinates.at( i ).x() * destinationCoordinates.at( i ).x();
    sumPyMy += sourceCoordinates.at( i ).y() * destinationCoordinates.at( i ).y();
    sumMx += destinationCoordinates.at( i ).x();
    sumMy += destinationCoordinates.at( i ).y();
  }

  const double deltaX = n * sumPx2 - std::pow( sumPx, 2 );
  const double deltaY = n * sumPy2 - std::pow( sumPy, 2 );

  const double aX = ( sumPx2 * sumMx - sumPx * sumPxMx ) / deltaX;
  const double aY = ( sumPy2 * sumMy - sumPy * sumPyMy ) / deltaY;
  const double bX = ( n * sumPxMx - sumPx * sumMx ) / deltaX;
  const double bY = ( n * sumPyMy - sumPy * sumMy ) / deltaY;

  origin.setX( aX );
  origin.setY( aY );

  pixelXSize = std::fabs( bX );
  pixelYSize = std::fabs( bY );
}


static void helmertImpl( const QVector<QgsPointXY> &sourceCoordinates, const QVector<QgsPointXY> &destinationCoordinates, QgsPointXY &origin, double &pixelSize, double &rotation )
{
#ifndef HAVE_GSL
  ( void ) sourceCoordinates;
  ( void ) destinationCoordinates;
  ( void ) origin;
  ( void ) pixelSize;
  ( void ) rotation;
  throw QgsNotSupportedException( QObject::tr( "Calculating a helmert transformation requires a QGIS build based GSL" ) );
#else
  const int n = destinationCoordinates.size();
  if ( n < 2 )
  {
    throw std::domain_error( QObject::tr( "Fit to a Helmert transform requires at least 2 points." ).toLocal8Bit().constData() );
  }

  double A = 0;
  double B = 0;
  double C = 0;
  double D = 0;
  double E = 0;
  double F = 0;
  double G = 0;
  double H = 0;
  double I = 0;
  double J = 0;
  for ( int i = 0; i < n; ++i )
  {
    A += sourceCoordinates.at( i ).x();
    B += sourceCoordinates.at( i ).y();
    C += destinationCoordinates.at( i ).x();
    D += destinationCoordinates.at( i ).y();
    E += destinationCoordinates.at( i ).x() * sourceCoordinates.at( i ).x();
    F += destinationCoordinates.at( i ).y() * sourceCoordinates.at( i ).y();
    G += std::pow( sourceCoordinates.at( i ).x(), 2 );
    H += std::pow( sourceCoordinates.at( i ).y(), 2 );
    I += destinationCoordinates.at( i ).x() * sourceCoordinates.at( i ).y();
    J += sourceCoordinates.at( i ).x() * destinationCoordinates.at( i ).y();
  }

  /* The least squares fit for the parameters { a, b, x0, y0 } is the solution
     to the matrix equation Mx = b, where M and b is given below. I *think*
     that this is correct but I derived it myself late at night. Look at
     helmert.jpg if you suspect bugs. */

  double MData[] = { A, -B, ( double ) n, 0., B, A, 0., ( double ) n, G + H, 0., A, B, 0., G + H, -B, A };

  double bData[] = { C, D, E + F, J - I };

  // we want to solve the equation M*x = b, where x = [a b x0 y0]
  gsl_matrix_view M = gsl_matrix_view_array( MData, 4, 4 );
  const gsl_vector_view b = gsl_vector_view_array( bData, 4 );
  gsl_vector *x = gsl_vector_alloc( 4 );
  gsl_permutation *p = gsl_permutation_alloc( 4 );
  int s;
  gsl_linalg_LU_decomp( &M.matrix, p, &s );
  gsl_linalg_LU_solve( &M.matrix, p, &b.vector, x );
  gsl_permutation_free( p );

  origin.setX( gsl_vector_get( x, 2 ) );
  origin.setY( gsl_vector_get( x, 3 ) );
  pixelSize = std::sqrt( std::pow( gsl_vector_get( x, 0 ), 2 ) + std::pow( gsl_vector_get( x, 1 ), 2 ) );
  rotation = std::atan2( gsl_vector_get( x, 1 ), gsl_vector_get( x, 0 ) );

  gsl_vector_free( x );
#endif
}

namespace
{
  //! second-moment conditioning of a point set (uncentered over centered), spread = RMS distance to the centroid
  bool momentConditioning( const QVector<QgsPointXY> &pts, double &cond, double &spread, double &maxAbs )
  {
    const int n = pts.size();
    if ( n < 2 )
      return false;
    long double cx = 0, cy = 0, raw = 0;
    maxAbs = 0;
    for ( const QgsPointXY &p : pts )
    {
      if ( !std::isfinite( p.x() ) || !std::isfinite( p.y() ) )
        return false;
      cx += p.x();
      cy += p.y();
      raw += static_cast<long double>( p.x() ) * p.x() + static_cast<long double>( p.y() ) * p.y();
      maxAbs = std::max( { maxAbs, std::fabs( p.x() ), std::fabs( p.y() ) } );
    }
    cx /= n;
    cy /= n;
    long double centered = 0;
    for ( const QgsPointXY &p : pts )
      centered += ( p.x() - cx ) * ( p.x() - cx ) + ( p.y() - cy ) * ( p.y() - cy );
    if ( !( centered > 0 ) )
      return false;
    cond = static_cast<double>( raw / centered );
    spread = std::sqrt( static_cast<double>( centered / n ) );
    return std::isfinite( cond ) && spread > 0;
  }
} // namespace

void QgsLeastSquares::linear( const QVector<QgsPointXY> &sourceCoordinates, const QVector<QgsPointXY> &destinationCoordinates, QgsPointXY &origin, double &pixelXSize, double &pixelYSize )
{
  linearImpl( sourceCoordinates, destinationCoordinates, origin, pixelXSize, pixelYSize );

  QgsSanitizer::run( "QG-GCP-002", [&] {
    // per-axis fit X = a + b x: the returned parameters must be a stationary point of the squared residuals
    const int n = destinationCoordinates.size();
    double cond, spread, maxAbs;
    if ( n < 2 || sourceCoordinates.size() != n || !momentConditioning( sourceCoordinates, cond, spread, maxAbs ) || !std::isfinite( origin.x() ) || !std::isfinite( origin.y() ) || !std::isfinite( pixelXSize ) || !std::isfinite( pixelYSize ) )
      return;
    for ( int axis = 0; axis < 2; ++axis )
    {
      long double mx = 0, mX = 0, cov = 0;
      for ( int i = 0; i < n; ++i )
      {
        mx += axis == 0 ? sourceCoordinates[i].x() : sourceCoordinates[i].y();
        mX += axis == 0 ? destinationCoordinates[i].x() : destinationCoordinates[i].y();
      }
      mx /= n;
      mX /= n;
      for ( int i = 0; i < n; ++i )
        cov += ( ( axis == 0 ? sourceCoordinates[i].x() : sourceCoordinates[i].y() ) - mx ) * ( ( axis == 0 ? destinationCoordinates[i].x() : destinationCoordinates[i].y() ) - mX );
      const double a = axis == 0 ? origin.x() : origin.y();
      const double b = ( cov < 0 ? -1.0 : 1.0 ) * ( axis == 0 ? pixelXSize : pixelYSize ); // the API returns the magnitude of the slope
      long double ga = 0, gb = 0, terms = 0;
      for ( int i = 0; i < n; ++i )
      {
        const long double x = axis == 0 ? sourceCoordinates[i].x() : sourceCoordinates[i].y();
        const long double X = axis == 0 ? destinationCoordinates[i].x() : destinationCoordinates[i].y();
        const long double r = a + b * x - X;
        ga += r;
        gb += r * x;
        terms += ( 1 + std::fabs( x ) ) * ( std::fabs( a ) + std::fabs( b * x ) + std::fabs( X ) );
      }
      const double g = QgsSanitizer::seen( "QG-GCP-002", static_cast<double>( std::max( std::fabs( ga ), std::fabs( gb ) ) ) );
      QgsSanitizer::check( "QG-GCP-002", g > 1e3 * QgsSanitizer::EPS * cond * static_cast<double>( terms ), "linear fit gradient " + QgsSanitizer::num( g ) + " on axis " + std::to_string( axis ) );
    }
  } );
}

void QgsLeastSquares::helmert( const QVector<QgsPointXY> &sourceCoordinates, const QVector<QgsPointXY> &destinationCoordinates, QgsPointXY &origin, double &pixelSize, double &rotation )
{
  helmertImpl( sourceCoordinates, destinationCoordinates, origin, pixelSize, rotation );

#ifdef HAVE_GSL
  QgsSanitizer::run( "QG-GCP-002", [&] {
    const int n = destinationCoordinates.size();
    double cond, spread, maxSrc;
    if ( n < 2 || sourceCoordinates.size() != n || !momentConditioning( sourceCoordinates, cond, spread, maxSrc ) || !std::isfinite( origin.x() ) || !std::isfinite( origin.y() ) || !std::isfinite( pixelSize ) || !std::isfinite( rotation ) )
      return;
    const long double a = pixelSize * std::cos( rotation ), b = pixelSize * std::sin( rotation );
    long double g[4] = { 0, 0, 0, 0 }, terms = 0;
    for ( int i = 0; i < n; ++i )
    {
      const long double x = sourceCoordinates[i].x(), y = sourceCoordinates[i].y(), X = destinationCoordinates[i].x(), Y = destinationCoordinates[i].y();
      const long double rx = a * x - b * y + origin.x() - X, ry = b * x + a * y + origin.y() - Y;
      g[0] += rx * x + ry * y;
      g[1] += -rx * y + ry * x;
      g[2] += rx;
      g[3] += ry;
      terms += ( 1 + std::fabs( x ) + std::fabs( y ) ) * ( std::fabs( X ) + std::fabs( Y ) + pixelSize * ( std::fabs( x ) + std::fabs( y ) ) + std::fabs( origin.x() ) + std::fabs( origin.y() ) );
    }
    double worst = 0;
    for ( const long double v : g )
      worst = std::max( worst, static_cast<double>( std::fabs( v ) ) );
    worst = QgsSanitizer::seen( "QG-GCP-002", worst );
    QgsSanitizer::check( "QG-GCP-002", worst > 1e3 * QgsSanitizer::EPS * cond * static_cast<double>( terms ), "Helmert fit gradient " + QgsSanitizer::num( worst ) );
  } );

  QgsSanitizer::run( "QG-GCP-003", [&] {
    // translating every destination point must translate the origin and leave scale and rotation unchanged
    const int n = destinationCoordinates.size();
    double cond, spread, maxSrc, condD, spreadD, maxDst;
    if ( n < 2 || sourceCoordinates.size() != n || !momentConditioning( sourceCoordinates, cond, spread, maxSrc ) || !momentConditioning( destinationCoordinates, condD, spreadD, maxDst ) || !( cond < 1e10 ) || !( pixelSize > 1e-12 ) || !std::isfinite( origin.x() ) || !std::isfinite( origin.y() ) )
      return;
    const double t = 10.0 * ( maxDst + std::hypot( origin.x(), origin.y() ) + 1.0 );
    QVector<QgsPointXY> moved;
    for ( const QgsPointXY &p : destinationCoordinates )
      moved << QgsPointXY( p.x() + t, p.y() - t );
    QgsPointXY origin2;
    double scale2 = 0, rotation2 = 0;
    helmertImpl( sourceCoordinates, moved, origin2, scale2, rotation2 );
    const double dSlope = 1e3 * QgsSanitizer::EPS * cond * ( maxDst + t ) / spread;
    double dRot = std::fabs( std::fmod( rotation2 - rotation, 2 * M_PI ) );
    dRot = std::min( dRot, 2 * M_PI - dRot );
    const double s = QgsSanitizer::seen( "QG-GCP-003", pixelSize );
    const double dOrigin = std::hypot( origin2.x() - origin.x() - t, origin2.y() - origin.y() + t );
    const bool bad = std::fabs( scale2 - s ) > dSlope || dRot > dSlope / pixelSize || dOrigin > 1e3 * QgsSanitizer::EPS * cond * ( t + std::hypot( origin.x(), origin.y() ) + maxDst ) + dSlope * maxSrc;
    QgsSanitizer::check( "QG-GCP-003", bad, "Helmert fit changed under destination translation: scale " + QgsSanitizer::num( s ) + " -> " + QgsSanitizer::num( scale2 ) + ", rotation " + QgsSanitizer::num( rotation ) + " -> " + QgsSanitizer::num( rotation2 ) );
  } );
#endif
}

#if 0
void QgsLeastSquares::affine( QVector<QgsPointXY> mapCoords,
                              QVector<QgsPointXY> pixelCoords )
{
  int n = mapCoords.size();
  if ( n < 4 )
  {
    throw std::domain_error( QObject::tr( "Fit to an affine transform requires at least 4 points." ).toLocal8Bit().constData() );
  }

  double A = 0, B = 0, C = 0, D = 0, E = 0, F = 0,
         G = 0, H = 0, I = 0, J = 0, K = 0;
  for ( int i = 0; i < n; ++i )
  {
    A += pixelCoords[i].x();
    B += pixelCoords[i].y();
    C += mapCoords[i].x();
    D += mapCoords[i].y();
    E += std::pow( pixelCoords[i].x(), 2 );
    F += std::pow( pixelCoords[i].y(), 2 );
    G += pixelCoords[i].x() * pixelCoords[i].y();
    H += pixelCoords[i].x() * mapCoords[i].x();
    I += pixelCoords[i].y() * mapCoords[i].y();
    J += pixelCoords[i].x() * mapCoords[i].y();
    K += mapCoords[i].x() * pixelCoords[i].y();
  }

  /* The least squares fit for the parameters { a, b, c, d, x0, y0 } is the
     solution to the matrix equation Mx = b, where M and b is given below.
     I *think* that this is correct but I derived it myself late at night.
     Look at affine.jpg if you suspect bugs. */

  double MData[] = { A,    B,    0,    0, ( double ) n,    0,
                     0,    0,    A,    B,    0, ( double ) n,
                     E,    G,    0,    0,    A,    0,
                     G,    F,    0,    0,    B,    0,
                     0,    0,    E,    G,    0,    A,
                     0,    0,    G,    F,    0,    B
                   };

  double bData[] = { C,    D,    H,    K,    J,    I };

  // we want to solve the equation M*x = b, where x = [a b c d x0 y0]
  gsl_matrix_view M = gsl_matrix_view_array( MData, 6, 6 );
  gsl_vector_view b = gsl_vector_view_array( bData, 6 );
  gsl_vector *x = gsl_vector_alloc( 6 );
  gsl_permutation *p = gsl_permutation_alloc( 6 );
  int s;
  gsl_linalg_LU_decomp( &M.matrix, p, &s );
  gsl_linalg_LU_solve( &M.matrix, p, &b.vector, x );
  gsl_permutation_free( p );

}
#endif

/**
 * Scales the given coordinates so that the center of gravity is at the origin and the mean distance to the origin is sqrt(2).
 *
 * Also returns 3x3 homogeneous matrices which can be used to normalize and de-normalize coordinates.
 */
void normalizeCoordinates( const QVector<QgsPointXY> &coords, QVector<QgsPointXY> &normalizedCoords, double normalizeMatrix[9], double denormalizeMatrix[9] )
{
  // Calculate center of gravity
  double cogX = 0.0, cogY = 0.0;
  for ( int i = 0; i < coords.size(); i++ )
  {
    cogX += coords[i].x();
    cogY += coords[i].y();
  }
  cogX *= 1.0 / coords.size();
  cogY *= 1.0 / coords.size();

  // Calculate mean distance to origin
  double meanDist = 0.0;
  for ( int i = 0; i < coords.size(); i++ )
  {
    const double X = ( coords[i].x() - cogX );
    const double Y = ( coords[i].y() - cogY );
    meanDist += std::sqrt( X * X + Y * Y );
  }
  meanDist *= 1.0 / coords.size();

  const double OOD = meanDist * M_SQRT1_2;
  const double D = 1.0 / OOD;
  normalizedCoords.resize( coords.size() );
  for ( int i = 0; i < coords.size(); i++ )
  {
    normalizedCoords[i] = QgsPointXY( ( coords[i].x() - cogX ) * D, ( coords[i].y() - cogY ) * D );
  }

  normalizeMatrix[0] = D;
  normalizeMatrix[1] = 0.0;
  normalizeMatrix[2] = -cogX * D;
  normalizeMatrix[3] = 0.0;
  normalizeMatrix[4] = D;
  normalizeMatrix[5] = -cogY * D;
  normalizeMatrix[6] = 0.0;
  normalizeMatrix[7] = 0.0;
  normalizeMatrix[8] = 1.0;

  denormalizeMatrix[0] = OOD;
  denormalizeMatrix[1] = 0.0;
  denormalizeMatrix[2] = cogX;
  denormalizeMatrix[3] = 0.0;
  denormalizeMatrix[4] = OOD;
  denormalizeMatrix[5] = cogY;
  denormalizeMatrix[6] = 0.0;
  denormalizeMatrix[7] = 0.0;
  denormalizeMatrix[8] = 1.0;
}

// Fits a homography to the given corresponding points, and
// return it in H (row-major format).
void QgsLeastSquares::projective( const QVector<QgsPointXY> &sourceCoordinates, const QVector<QgsPointXY> &destinationCoordinates, double H[9] )
{
#ifndef HAVE_GSL
  ( void ) sourceCoordinates;
  ( void ) destinationCoordinates;
  ( void ) H;
  throw QgsNotSupportedException( QObject::tr( "Calculating a projective transformation requires a QGIS build based GSL" ) );
#else
  Q_ASSERT( sourceCoordinates.size() == destinationCoordinates.size() );

  if ( destinationCoordinates.size() < 4 )
  {
    throw std::domain_error( QObject::tr( "Fitting a projective transform requires at least 4 corresponding points." ).toLocal8Bit().constData() );
  }

  QVector<QgsPointXY> sourceCoordinatesNormalized;
  QVector<QgsPointXY> destinationCoordinatesNormalized;

  double normSource[9], denormSource[9];
  double normDest[9], denormDest[9];
  normalizeCoordinates( sourceCoordinates, sourceCoordinatesNormalized, normSource, denormSource );
  normalizeCoordinates( destinationCoordinates, destinationCoordinatesNormalized, normDest, denormDest );

  // GSL does not support a full SVD, so we artificially add a linear dependent row
  // to the matrix in case the system is underconstrained.
  const uint m = std::max( 9u, ( uint ) destinationCoordinatesNormalized.size() * 2u );
  const uint n = 9;
  gsl_matrix *S = gsl_matrix_alloc( m, n );

  for ( int i = 0; i < destinationCoordinatesNormalized.size(); i++ )
  {
    gsl_matrix_set( S, i * 2, 0, sourceCoordinatesNormalized[i].x() );
    gsl_matrix_set( S, i * 2, 1, sourceCoordinatesNormalized[i].y() );
    gsl_matrix_set( S, i * 2, 2, 1.0 );

    gsl_matrix_set( S, i * 2, 3, 0.0 );
    gsl_matrix_set( S, i * 2, 4, 0.0 );
    gsl_matrix_set( S, i * 2, 5, 0.0 );

    gsl_matrix_set( S, i * 2, 6, -destinationCoordinatesNormalized[i].x() * sourceCoordinatesNormalized[i].x() );
    gsl_matrix_set( S, i * 2, 7, -destinationCoordinatesNormalized[i].x() * sourceCoordinatesNormalized[i].y() );
    gsl_matrix_set( S, i * 2, 8, -destinationCoordinatesNormalized[i].x() * 1.0 );

    gsl_matrix_set( S, i * 2 + 1, 0, 0.0 );
    gsl_matrix_set( S, i * 2 + 1, 1, 0.0 );
    gsl_matrix_set( S, i * 2 + 1, 2, 0.0 );

    gsl_matrix_set( S, i * 2 + 1, 3, sourceCoordinatesNormalized[i].x() );
    gsl_matrix_set( S, i * 2 + 1, 4, sourceCoordinatesNormalized[i].y() );
    gsl_matrix_set( S, i * 2 + 1, 5, 1.0 );

    gsl_matrix_set( S, i * 2 + 1, 6, -destinationCoordinatesNormalized[i].y() * sourceCoordinatesNormalized[i].x() );
    gsl_matrix_set( S, i * 2 + 1, 7, -destinationCoordinatesNormalized[i].y() * sourceCoordinatesNormalized[i].y() );
    gsl_matrix_set( S, i * 2 + 1, 8, -destinationCoordinatesNormalized[i].y() * 1.0 );
  }

  if ( destinationCoordinatesNormalized.size() == 4 )
  {
    // The GSL SVD routine only supports matrices with rows >= columns (m >= n)
    // Unfortunately, we can't use the SVD of the transpose (i.e. S^T = (U D V^T)^T = V D U^T)
    // to work around this, because the solution lies in the right nullspace of S, and
    // gsl only supports a thin SVD of S^T, which does not return these vectors.

    // HACK: duplicate last row to get a 9x9 equation system
    for ( int j = 0; j < 9; j++ )
    {
      gsl_matrix_set( S, 8, j, gsl_matrix_get( S, 7, j ) );
    }
  }

  // Solve Sh = 0 in the total least squares sense, i.e.
  // with Sh = min and |h|=1. The solution "h" is given by the
  // right singular eigenvector of S corresponding, to the smallest
  // singular value (via SVD).
  gsl_matrix *V = gsl_matrix_alloc( n, n );
  gsl_vector *singular_values = gsl_vector_alloc( n );
  gsl_vector *work = gsl_vector_alloc( n );

  // V = n x n
  // U = m x n (thin SVD)  U D V^T
  gsl_linalg_SV_decomp( S, V, singular_values, work );

  // Columns of V store the right singular vectors of S
  for ( unsigned int i = 0; i < n; i++ )
  {
    H[i] = gsl_matrix_get( V, i, n - 1 );
  }

  gsl_matrix *prodMatrix = gsl_matrix_alloc( 3, 3 );

  gsl_matrix_view Hmatrix = gsl_matrix_view_array( H, 3, 3 );
  const gsl_matrix_view normSourceMatrix = gsl_matrix_view_array( normSource, 3, 3 );
  const gsl_matrix_view denormDestMatrix = gsl_matrix_view_array( denormDest, 3, 3 );

  // Change coordinate frame of image and pre-image from normalized to destination and source coordinates.
  // H' = denormalizeMapCoords*H*normalizePixelCoords
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, &Hmatrix.matrix, &normSourceMatrix.matrix, 0.0, prodMatrix );
  gsl_blas_dgemm( CblasNoTrans, CblasNoTrans, 1.0, &denormDestMatrix.matrix, prodMatrix, 0.0, &Hmatrix.matrix );

  gsl_matrix_free( prodMatrix );
  gsl_matrix_free( S );
  gsl_matrix_free( V );
  gsl_vector_free( singular_values );
  gsl_vector_free( work );
#endif
}
