/***************************************************************************
  qgssanitizer.h
  ---------------
  Scientific sanitizer runtime (research instrumentation, header only).

  Inactive unless the environment variable QGIS_SANITIZE_LOG names a file.
  Each alarm appends one JSON line {"id":..,"detail":..}; a checker that
  throws appends to "<log>.err" and is otherwise swallowed. Checks observe
  state and never change results. QGIS_SANITIZE_FAULT=<ID> perturbs the
  value that checker observes (instrumentation self-test only).
 ***************************************************************************/
#ifndef QGSSANITIZER_H
#define QGSSANITIZER_H

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <string>

namespace QgsSanitizer
{
  inline const char *logPath()
  {
    static const char *p = std::getenv( "QGIS_SANITIZE_LOG" );
    return ( p && *p ) ? p : nullptr;
  }

  inline bool enabled() { return logPath() != nullptr; }

  inline bool &busy()
  {
    thread_local bool b = false;
    return b;
  }

  struct Guard
  {
    bool prev;
    Guard() : prev( busy() ) { busy() = true; }
    ~Guard() { busy() = prev; }
    Guard( const Guard & ) = delete;
    Guard &operator=( const Guard & ) = delete;
  };

  inline void append( const std::string &path, const std::string &line )
  {
    if ( FILE *f = std::fopen( path.c_str(), "a" ) )
    {
      std::fputs( line.c_str(), f );
      std::fputc( '\n', f );
      std::fclose( f );
    }
  }

  inline std::string esc( const std::string &s )
  {
    std::string o;
    for ( const char c : s )
    {
      if ( c == '"' || c == '\\' )
      {
        o += '\\';
        o += c;
      }
      else if ( static_cast<unsigned char>( c ) >= 32 )
        o += c;
    }
    return o;
  }

  inline void alarm( const char *id, const std::string &detail )
  {
    if ( const char *p = logPath() )
      append( p, std::string( "{\"id\":\"" ) + id + "\",\"detail\":\"" + esc( detail ) + "\"}" );
  }

  inline std::string num( double v )
  {
    char b[48];
    std::snprintf( b, sizeof b, "%.17g", v );
    return b;
  }

  //! True if QGIS_SANITIZE_FAULT selects this checker (self-test only).
  inline bool fault( const char *id )
  {
    static const char *p = std::getenv( "QGIS_SANITIZE_FAULT" );
    return p && std::string( p ) == id;
  }

  //! Value the checker observes: unchanged unless the fault is selected.
  inline double seen( const char *id, double v ) { return fault( id ) ? v * 1.1 + 1.0 : v; }

  //! Runs a checker body: disabled, re-entrant and exception safe.
  template<class F> inline void run( const char *id, F &&body )
  {
    if ( !enabled() || busy() )
      return;
    Guard g;
    try
    {
      body();
    }
    catch ( const std::exception &e )
    {
      append( std::string( logPath() ) + ".err", std::string( id ) + ": " + e.what() );
    }
    catch ( ... )
    {
      append( std::string( logPath() ) + ".err", std::string( id ) + ": unknown exception" );
    }
  }

  //! Records an alarm when `violated` and everything observed was finite.
  inline void check( const char *id, bool violated, const std::string &detail )
  {
    // optional coverage log: one line per evaluated predicate (QGIS_SANITIZE_REACH=<file>)
    static const char *reach = std::getenv( "QGIS_SANITIZE_REACH" );
    if ( reach && *reach )
      append( reach, id );
    if ( violated )
      alarm( id, detail );
  }

  constexpr double EPS = 2.220446049250313e-16;

  /**
   * QG-GEO-007: a planar projection of a point by `distance` along `azimuthDeg` (clockwise from north) must be
   * at that distance and azimuth from the origin. Shared by every implementation of the projection.
   */
  inline void checkProject( double ax, double ay, double rx, double ry, double distance, double azimuthDeg )
  {
    run( "QG-GEO-007", [&] {
      if ( !( distance > 0 ) || !std::isfinite( distance ) || !( std::fabs( azimuthDeg ) <= 720.0 ) || !( distance >= 1e-200 ) || !std::isfinite( azimuthDeg ) || !std::isfinite( ax ) || !std::isfinite( ay ) || !std::isfinite( rx ) || !std::isfinite( ry ) )
        return;
      const double ddx = rx - ax, ddy = ry - ay;
      double az = std::atan2( ddx, ddy ) * 180.0 / 3.14159265358979323846 - azimuthDeg;
      az -= 360.0 * std::round( az / 360.0 ); // circular difference in [-180, 180]
      const double scale = std::max( std::fabs( ax ), std::fabs( ay ) );
      const double dist = seen( "QG-GEO-007", std::hypot( ddx, ddy ) );
      const double tolAz = 64 * EPS * ( 1.0 + scale / distance ) * 180.0 / 3.14159265358979323846 + 1e-12 * std::fabs( azimuthDeg );
      const double tolD = 64 * EPS * ( scale + distance );
      check( "QG-GEO-007", std::fabs( az ) > tolAz || std::fabs( dist - distance ) > tolD, "azimuth error " + num( az ) + " deg, distance " + num( dist ) + " vs " + num( distance ) );
    } );
  }
  constexpr double EPS32 = 1.1920928955078125e-07;
} // namespace QgsSanitizer

#endif
