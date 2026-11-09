/***************************************************************************
                              qgstininterpolator.cpp
                              ----------------------
  begin                : March 10, 2008
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

#include "qgstininterpolator.h"
#include "qgssanitizer.h"

#include "CloughTocherInterpolator.h"
#include "LinTriangleInterpolator.h"
#include "NormVecDecorator.h"
#include "qgscurve.h"
#include "qgscurvepolygon.h"
#include "qgsdualedgetriangulation.h"
#include "qgsfeature.h"
#include "qgsfeatureiterator.h"
#include "qgsfeedback.h"
#include "qgsgeometry.h"
#include "qgsmulticurve.h"
#include "qgsmultisurface.h"
#include "qgspoint.h"
#include "qgsvariantutils.h"
#include "qgsvectorlayer.h"

QgsTinInterpolator::QgsTinInterpolator( const QList<LayerData> &inputData, QgsTinInterpolator::TinInterpolation interpolation, QgsFeedback *feedback )
  : QgsInterpolator( inputData )
  , mFeedback( feedback )
  , mInterpolation( interpolation )
{}

QgsTinInterpolator::~QgsTinInterpolator()
{
  delete mTriangulation;
  delete mTriangleInterpolator;
}

int QgsTinInterpolator::interpolatePoint( double x, double y, double &result, QgsFeedback * )
{
  if ( !mIsInitialized )
  {
    initialize();
  }

  if ( !mTriangleInterpolator )
  {
    return 1;
  }

  QgsPoint r( 0, 0, 0 );
  if ( !mTriangleInterpolator->calcPoint( x, y, r ) )
  {
    return 2;
  }
  result = r.z();

  QgsSanitizer::run( "QG-INT-003", [&] {
    // piecewise-linear interpolation must reproduce a plane exactly
    if ( mInterpolation != QgsTinInterpolator::TinInterpolation::Linear || !mTriangulation || !std::isfinite( result ) )
      return;
    const int n = mTriangulation->pointsCount();
    if ( n < 3 || n > 5000 )
      return;
    long double mx = 0, my = 0, mz = 0;
    for ( int i = 0; i < n; ++i )
    {
      const QgsPoint *p = mTriangulation->point( i );
      if ( !p || !std::isfinite( p->x() ) || !std::isfinite( p->y() ) || !std::isfinite( p->z() ) )
        return;
      mx += p->x();
      my += p->y();
      mz += p->z();
    }
    mx /= n;
    my /= n;
    mz /= n;
    long double sxx = 0, sxy = 0, syy = 0, sxz = 0, syz = 0;
    double zAbs = 0, coordAbs = 0;
    double xmin = 1e300, xmax = -1e300, ymin = 1e300, ymax = -1e300;
    for ( int i = 0; i < n; ++i )
    {
      const QgsPoint *p = mTriangulation->point( i );
      const long double dx = p->x() - mx, dy = p->y() - my, dz = p->z() - mz;
      sxx += dx * dx;
      sxy += dx * dy;
      syy += dy * dy;
      sxz += dx * dz;
      syz += dy * dz;
      zAbs = std::max( zAbs, std::fabs( p->z() ) );
      coordAbs = std::max( { coordAbs, std::fabs( p->x() ), std::fabs( p->y() ) } );
      xmin = std::min( xmin, p->x() );
      xmax = std::max( xmax, p->x() );
      ymin = std::min( ymin, p->y() );
      ymax = std::max( ymax, p->y() );
    }
    const long double det = sxx * syy - sxy * sxy;
    if ( !( det > 1e-12L * ( sxx + syy ) * ( sxx + syy ) ) )
      return; // collinear or degenerate sample
    const double a = static_cast<double>( ( sxz * syy - syz * sxy ) / det );
    const double b = static_cast<double>( ( syz * sxx - sxz * sxy ) / det );
    const double extent = std::max( xmax - xmin, ymax - ymin );
    const double scale = zAbs + ( std::fabs( a ) + std::fabs( b ) ) * ( coordAbs + extent );
    for ( int i = 0; i < n; ++i )
    {
      const QgsPoint *p = mTriangulation->point( i );
      if ( std::fabs( p->z() - static_cast<double>( mz + a * ( p->x() - mx ) + b * ( p->y() - my ) ) ) > 64 * QgsSanitizer::EPS * scale )
        return; // data are not planar
    }
    const double expected = static_cast<double>( mz + a * ( x - mx ) + b * ( y - my ) );
    const double v = QgsSanitizer::seen( "QG-INT-003", result );
    QgsSanitizer::check( "QG-INT-003", std::fabs( v - expected ) > 1e3 * QgsSanitizer::EPS * scale * ( 1.0 + coordAbs / std::max( extent, 1e-300 ) ), "TIN value " + QgsSanitizer::num( v ) + " vs plane " + QgsSanitizer::num( expected ) );
  } );
  return 0;
}

QgsFields QgsTinInterpolator::triangulationFields()
{
  return QgsTriangulation::triangulationFields();
}

void QgsTinInterpolator::setTriangulationSink( QgsFeatureSink *sink )
{
  mTriangulationSink = sink;
}

void QgsTinInterpolator::initialize()
{
  QgsDualEdgeTriangulation *dualEdgeTriangulation = new QgsDualEdgeTriangulation( 100000 );
  if ( mInterpolation == QgsTinInterpolator::TinInterpolation::CloughTocher )
  {
    NormVecDecorator *dec = new NormVecDecorator();
    dec->addTriangulation( dualEdgeTriangulation );
    mTriangulation = dec;
  }
  else
  {
    mTriangulation = dualEdgeTriangulation;
  }

  //get number of features if we use a progress bar
  long long nFeatures = 0;
  long long nProcessedFeatures = 0;
  if ( mFeedback )
  {
    for ( const LayerData &layer : std::as_const( mLayerData ) )
    {
      if ( layer.source )
      {
        nFeatures += layer.source->featureCount();
      }
    }
  }

  const QgsCoordinateReferenceSystem crs = !mLayerData.empty() ? mLayerData.at( 0 ).source->sourceCrs() : QgsCoordinateReferenceSystem();

  QgsFeature f;
  for ( const LayerData &layer : std::as_const( mLayerData ) )
  {
    if ( layer.source )
    {
      QgsAttributeList attList;
      switch ( layer.valueSource )
      {
        case QgsInterpolator::ValueSource::Attribute:
          attList.push_back( layer.interpolationAttribute );
          break;

        case QgsInterpolator::ValueSource::M:
        case QgsInterpolator::ValueSource::Z:
          break;
      }

      QgsFeatureIterator fit = layer.source->getFeatures( QgsFeatureRequest().setSubsetOfAttributes( attList ).setDestinationCrs( crs, layer.transformContext ) );

      while ( fit.nextFeature( f ) )
      {
        if ( mFeedback )
        {
          if ( mFeedback->isCanceled() )
          {
            break;
          }
          if ( nFeatures > 0 )
            mFeedback->setProgress( 100.0 * static_cast<double>( nProcessedFeatures ) / nFeatures );
        }
        insertData( f, layer.valueSource, layer.interpolationAttribute, layer.sourceType );
        ++nProcessedFeatures;
      }
    }
  }

  if ( mInterpolation == QgsTinInterpolator::TinInterpolation::CloughTocher )
  {
    NormVecDecorator *dec = dynamic_cast<NormVecDecorator *>( mTriangulation );
    if ( dec )
    {
      auto ctInterpolator = std::make_unique<CloughTocherInterpolator>();
      dec->estimateFirstDerivatives( mFeedback );
      ctInterpolator->setTriangulation( dec );
      mTriangleInterpolator = ctInterpolator.release();
      dec->setTriangleInterpolator( mTriangleInterpolator );
    }
  }
  else //linear
  {
    mTriangleInterpolator = new LinTriangleInterpolator( dualEdgeTriangulation );
  }
  mIsInitialized = true;

  //debug
  if ( mTriangulationSink )
  {
    dualEdgeTriangulation->saveTriangulation( mTriangulationSink, mFeedback );
  }
}

int QgsTinInterpolator::insertData( const QgsFeature &f, QgsInterpolator::ValueSource source, int attr, QgsInterpolator::SourceType type )
{
  QgsGeometry g = f.geometry();
  if ( g.isNull() || g.isEmpty() )
  {
    return 2;
  }

  //check attribute value
  double attributeValue = 0;
  bool attributeConversionOk = false;
  switch ( source )
  {
    case QgsInterpolator::ValueSource::Attribute:
    {
      QVariant attributeVariant = f.attribute( attr );
      if ( QgsVariantUtils::isNull( attributeVariant ) ) //attribute not found, something must be wrong (e.g. NULL value)
      {
        return 3;
      }
      attributeValue = attributeVariant.toDouble( &attributeConversionOk );
      if ( !attributeConversionOk || std::isnan( attributeValue ) ) //don't consider vertices with attributes like 'nan' for the interpolation
      {
        return 4;
      }
      break;
    }

    case QgsInterpolator::ValueSource::M:
      if ( !g.constGet()->isMeasure() )
        return 3;
      else
        break;

    case QgsInterpolator::ValueSource::Z:
      if ( !g.constGet()->is3D() )
        return 3;
      else
        break;
  }

  switch ( type )
  {
    case QgsInterpolator::SourceType::Points:
    {
      if ( addPointsFromGeometry( g, source, attributeValue ) != 0 )
        return -1;
      break;
    }

    case QgsInterpolator::SourceType::BreakLines:
    case QgsInterpolator::SourceType::StructureLines:
    {
      switch ( QgsWkbTypes::geometryType( g.wkbType() ) )
      {
        case Qgis::GeometryType::Point:
        {
          if ( addPointsFromGeometry( g, source, attributeValue ) != 0 )
            return -1;
          break;
        }

        case Qgis::GeometryType::Line:
        case Qgis::GeometryType::Polygon:
        {
          // need to extract all rings from input geometry
          std::vector<const QgsCurve *> curves;
          if ( QgsWkbTypes::geometryType( g.wkbType() ) == Qgis::GeometryType::Polygon )
          {
            std::vector<const QgsCurvePolygon *> polygons;
            if ( g.isMultipart() )
            {
              const QgsMultiSurface *ms = qgsgeometry_cast<const QgsMultiSurface *>( g.constGet() );
              for ( int i = 0; i < ms->numGeometries(); ++i )
              {
                polygons.emplace_back( qgsgeometry_cast<const QgsCurvePolygon *>( ms->geometryN( i ) ) );
              }
            }
            else
            {
              polygons.emplace_back( qgsgeometry_cast<const QgsCurvePolygon *>( g.constGet() ) );
            }

            for ( const QgsCurvePolygon *polygon : polygons )
            {
              if ( !polygon )
                continue;

              if ( polygon->exteriorRing() )
                curves.emplace_back( polygon->exteriorRing() );

              for ( int i = 0; i < polygon->numInteriorRings(); ++i )
              {
                curves.emplace_back( polygon->interiorRing( i ) );
              }
            }
          }
          else
          {
            if ( g.isMultipart() )
            {
              const QgsMultiCurve *mc = qgsgeometry_cast<const QgsMultiCurve *>( g.constGet() );
              for ( int i = 0; i < mc->numGeometries(); ++i )
              {
                curves.emplace_back( mc->curveN( i ) );
              }
            }
            else
            {
              curves.emplace_back( qgsgeometry_cast<const QgsCurve *>( g.constGet() ) );
            }
          }

          for ( const QgsCurve *curve : curves )
          {
            if ( !curve )
              continue;

            QgsPointSequence linePoints;
            curve->points( linePoints );
            for ( QgsPoint &point : linePoints )
            {
              switch ( source )
              {
                case QgsInterpolator::ValueSource::Attribute:
                  if ( point.is3D() )
                    point.setZ( attributeValue );
                  else
                    point.addZValue( attributeValue );
                  break;

                case QgsInterpolator::ValueSource::M:
                  if ( point.is3D() )
                    point.setZ( point.m() );
                  else
                    point.addZValue( point.m() );
                  break;

                case QgsInterpolator::ValueSource::Z:
                  break;
              }
            }
            mTriangulation->addLine( linePoints, type );
          }
          break;
        }
        case Qgis::GeometryType::Unknown:
        case Qgis::GeometryType::Null:
          break;
      }
      break;
    }
  }

  return 0;
}


int QgsTinInterpolator::addPointsFromGeometry( const QgsGeometry &g, QgsInterpolator::ValueSource source, double attributeValue )
{
  // loop through all vertices and add to triangulation
  for ( auto point = g.vertices_begin(); point != g.vertices_end(); ++point )
  {
    QgsPoint p = *point;
    double z = 0;
    switch ( source )
    {
      case QgsInterpolator::ValueSource::Attribute:
        z = attributeValue;
        break;

      case QgsInterpolator::ValueSource::Z:
        z = p.z();
        break;

      case QgsInterpolator::ValueSource::M:
        z = p.m();
        break;
    }
    if ( mTriangulation->addPoint( QgsPoint( p.x(), p.y(), z ) ) == -100 )
    {
      return -1;
    }
  }
  return 0;
}
