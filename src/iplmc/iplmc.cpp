#include "../iplmc.hpp"

#include <cmath>
#include <boost/foreach.hpp>
#include <boost/random/linear_congruential.hpp>
#include <boost/random/uniform_01.hpp>

#include "../diagnostic/tracer.hpp"

namespace {

  using namespace iplmcfd;

  void taylor_green_particle( 
    const mix_t& mix, real_t x, real_t y, Particle& p )
  {

    // static types
    static bool first_call(true);
    static cpr_t::phi phi_fuel( mix ), phi_air( mix ), phi_equil( mix );
    static mix_t::pressure pressure;

    // init static types
    if ( first_call )
    {
      // pressure is 1 atm
      pressure = mix.atm2dynescc( 1 );

      typedef mix_t::molefrac X_t;
      typedef std::vector< X_t > X_a;

      // equilibrium composition
      X_a Xequil( mix.nspec(), X_t(0.) );
      mix_t::temperature Tequil( 300 );                     // T = 300 K
      Xequil[mix.index("CH4")] = 0.5;     // Reaction:
      Xequil[mix.index("O2")] = 1.0;      // .5CH4 + O2 + 3.76 --> Products
      Xequil[mix.index("N2")] = 3.76;
      mix.normalize( Xequil ); 
      mix.equil( Xequil, pressure, Tequil, mix_t::PH ); 
      mix.convert( Xequil, phi_equil.Y() ); // mole fractions to mass fractions
      phi_equil.T() = Tequil;

      // fuel particle
      X_a Xfuel( mix.nspec(), X_t(0.) );
      mix_t::temperature Tfuel( 1500 );
      Xfuel[mix.index("CH4")] = 0.6;
      Xfuel[mix.index("O2")]  = 0.4 * 0.210084;
      Xfuel[mix.index("N2")]  = 0.4 * ( 1 - 0.210084 );
      mix.normalize( Xfuel );
      mix.equil( Xfuel, pressure, Tfuel, mix_t::PH ); 
      mix.normalize( Xfuel );
      mix.convert( Xfuel, phi_fuel.Y() );
      phi_fuel.T() = Tfuel;

      // air particle
      X_a Xair( mix.nspec(), X_t(0.) );
      Xair[mix.index("O2")]  = 0.210084; 
      Xair[mix.index("N2")]  = 1 - 0.210084;
      mix.normalize( Xair );
      mix.convert( Xair, phi_air.Y() );
      phi_air.T() = Tfuel;

      first_call = false;
    }

    /*
     * 2pi +----------+----------+
     *     |          |          |
     *     | fuel/air | air      |
     *     | @ 300 K  | @ 1000 K |
     *     |          |          |
     *     +----------+----------+
     *     |          |          |
     *     | air      | fuel/air |
     *     | @ 1000 K | @ 300 K  |
     *     |          |          |
     *   0 +----------+----------+
     *     0                    2pi
     */

    /*
    real_t pi = 4*std::atan( 1. );
    if ( ( x - pi )*( y - pi ) > 0 ) // air composition
    {
      p.f = 0.;
      p.phi = phi_air;
      p.mass =  mix.rho( pressure, p.phi.T(), p.phi.Y() );
    }
    else // fuel/air mixture
    {
      p.f = 1.;
      p.phi = phi_equil;
      p.mass =  mix.rho( pressure, p.phi.T(), p.phi.Y() );
    }
    */

    /*
     * 2pi    +------------------------+
     *        |  eq  |   fuel   |  eq  |
     *        |      |          |      |
     * 3pi/4  +------+----------+------+
     *        |      |          |      |
     *        | fuel |  equil   | fuel |
     *        |      |          |      |
     *        |      |          |      |
     * pi/2   +------+----------+------+
     *        |      |          |      |
     *        |  eq  |   fuel   |  eq  |
     * 0      +------+----------+------+
     *        0     pi/2      3pi/2    2pi
     */

    // it's ugly, but it's an initialization!
    real_t pi = 4*std::atan( 1. );
    if ( y < ( pi / 2 ) ) // bottom row
    {

      if ( x < ( pi / 2 ) || x > ( 3 * pi / 2 ) ) // left, right columns
      {
        p.f = 0.;
        p.phi = phi_air;
        p.mass =  mix.rho( pressure, p.phi.T(), p.phi.Y() );
      }
      else // middle column
      {
        p.f = 1.;
        p.phi = phi_fuel;
        p.mass =  mix.rho( pressure, p.phi.T(), p.phi.Y() );
      }

    }
    else if ( y > ( 3 * pi / 2 ) ) // top row
    {

      if ( x < ( pi / 2 ) || x > ( 3 * pi / 2 ) ) // left, right columns
      {
        p.f = 0.;
        p.phi = phi_air;
        p.mass =  mix.rho( pressure, p.phi.T(), p.phi.Y() );
      }
      else // middle column
      {
        p.f = 1.;
        p.phi = phi_fuel;
        p.mass =  mix.rho( pressure, p.phi.T(), p.phi.Y() );
      } 

    }
    else // middle row
    {

      if ( x < ( pi / 2 ) || x > ( 3 * pi / 2 ) ) // left, right columns
      {
        p.f = 1.;
        p.phi = phi_fuel;
        p.mass =  mix.rho( pressure, p.phi.T(), p.phi.Y() );
      }
      else // middle column
      {
        p.f = 0.;
        p.phi = phi_air;
        p.mass =  mix.rho( pressure, p.phi.T(), p.phi.Y() );
      } 

    }

  }

} // unnamed namespace

namespace iplmcfd {

  // +--------------+
  // | Iplmc::Iplmc |
  // +--------------+
  Iplmc::Iplmc( const simparam& p, const Domain& d, const Subdomain& s )
    : m_simparam( p ), 
    m_domain(d), 
    m_subdomain(s),
    m_dt( 0 ),
    m_time( 0 )
  {
    tracer::scope _( "Iplmc::Iplmc" );

    m_chem.reset( new mix_t( p.chem, p.chem_therm, p.chem_tran ) );
    m_cpr.reset( new cpr_t( *m_chem ) );
    
    // set some chemkinpp parameters
    m_cpr->solver( p.chem_solver );
    m_cpr->reltol( p.chem_reltol );
    m_cpr->abstol( p.chem_abstol );
    
    // set cloning and clustering limits
    m_min_ppc  = p.ppc[0];
    m_mean_ppc = p.ppc[1];
    m_max_ppc  = p.ppc[2];

    /* The Taylor-Green Analytic solution is from the solution of the 2D
     * incompressible Navier-Stokes equations.  As such, we get u and v and p 
     * from direct function calls.  p is a mechanical pressure and thus not 
     * really relevant to the MC solver at all.  u and v are dependent on nu.  
     * nu varies with temperature, but in the derivation of the analytic 
     * solution it must be constant.  So we calculate a reference nu here.
     */

    // thermodynamic reference pressure is just 1 atm
    m_pref = m_chem->atm2dynescc( 1. );

    //! \todo remove hard-coding of reference viscosity
    // let's use air at some average temp (equil + air) / 2
    // ( 2200 K + 300 K ) / 2 = 1250 K
    std::vector< mix_t::molefrac > Xair( m_chem->nspec(), mix_t::molefrac( 0. ) );
    Xair[m_chem->index("O2")]  = 0.210084; 
    Xair[m_chem->index("N2")]  = 1 - 0.210084;
    m_chem->normalize( Xair );
    mix_t::temperature Tair( 1250. );
    real_t mu = m_chem->mu( Tair, Xair );
    real_t rho = m_chem->rho( m_pref, Tair, Xair );

    m_nuref = mu / rho;

    m_tgv.reset( new TaylorGreenVortex( p.Re, 1, m_nuref ) );

    // allocate cell map
    BOOST_FOREACH( Subdomain::gidlid_mp::const_reference gidlid, 
      m_subdomain.homes() 
      )
    {
      Cell cell( *m_chem );
      cell.ctr.reset();
      m_cells.insert( cell_mp::value_type( gidlid.first, cell ) );
    }
      
    // allocate neighbor cell map
    BOOST_FOREACH( Subdomain::gidrank_mp::const_reference gidrank, 
      m_subdomain.neighs() )
    {
      Cell cell( *m_chem );
      cell.ctr.reset();
      m_neigh_cells.insert( cell_mp::value_type( gidrank.first, cell ) );
    }

    m_p2p.reset( 
      new Particle_p2p( m_subdomain, m_cells, m_neigh_cells, *m_chem ) );

  }

  // +--------------+
  // | Iplmc::omega |
  // +--------------+
  real_t Iplmc::omega( const Cell& cell )
  {
    // mu
    static std::vector< mix_t::molefrac > X( cell.scalars.phi.Y().size() );
    m_chem->convert( cell.scalars.phi.Y(), X );
    real_t mu = m_chem->mu( cell.scalars.phi.T(), X );

    // density
    real_t rho = m_pref / cell.scalars.RT;

    real_t numer = simparam::instance().Comega * mu;
    real_t denom = rho * simparam::instance().Sc * 
      m_domain.delta() * m_domain.delta();

    return numer / denom;
  }

  // +------------------+
  // | Iplmc::update_dt |
  // +------------------+
  real_t Iplmc::update_dt()
  {
    tracer::scope _( "Iplmc::update_dt" );

    static real_t Re_lambda = m_simparam.Re_lambda;
    static real_t stable = m_simparam.stable;
    static const real_t& dx = m_domain.dx()[XDIM];
    static const real_t& dy = m_domain.dx()[YDIM];
    static real_t inv_dx2 = 1 / ( dx*dx ) + 1 / ( dy*dy );

    m_ke = 0;
    m_dt = MAX_REAL;
    BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
    {
      Cell& cell = gidcell.second;
      real_v vel;

      m_tgv->uvw( m_time, m_domain.xyz( gidcell.first ), vel );

      real_t denom = 
        std::abs( vel[XDIM] ) / dx + 
        std::abs( vel[YDIM] ) / dy + 
        Re_lambda * m_nuref * inv_dx2;

      real_t dt = 1 / denom;
      m_dt = std::min( m_dt, dt );  

      // This is a good opportunity to log the integral of kinetic energy
      m_ke += .5 * ( vel[XDIM] * vel[XDIM] +
        vel[YDIM] * vel[YDIM] ) * m_domain.dV();
    }

    MPI_Allreduce( MPI_IN_PLACE, &m_dt, 1, mpi_real, MPI_MIN, MPI_COMM_WORLD );
    MPI_Allreduce( MPI_IN_PLACE, &m_ke, 1, mpi_real, MPI_SUM, MPI_COMM_WORLD );

    // The z-direction is bogus, so lets normalize the KE sum
    m_ke = m_ke / m_domain.shape()[ZDIM];

    m_dt = stable*m_dt;

    return m_dt;

  }

  // +-------------------+
  // | Iplmc::initialize |
  // +-------------------+
  void Iplmc::initialize()
  {
    tracer::scope _( "Iplmc::initialize" );

    // Create random generator
    typedef boost::rand48  rng_t;
    boost::uniform_01< rng_t, real_t > 
      rand( rng_t( rng_t::result_type( m_subdomain.my_rank() ) ) );

    Particle p( *m_chem );
    BOOST_FOREACH( cell_mp::reference gidcell, m_cells )
    {
      gid_t I = gidcell.first;
      Cell& cell = gidcell.second;
      
      // get cell mid-point coordinate (in non-dim real x coordinate)
      real_v xloc = m_domain.xyz( m_domain.gid2ind( I ) );
      
      taylor_green_particle( *m_chem, xloc[XDIM], xloc[YDIM], p );

      // Add particles
      for( int j=0; j!=m_mean_ppc; ++j)
      {
        cell.contained().push_back( p );
        cell.contained().back().X = rand(), rand(), rand();
        cell.contained().back().mass /= m_mean_ppc;
        ++(cell.ctr.camein);
      }

    }

    // scalar average here to output the initial condition
    scalar_average();

  }

  // +------------------------+
  // | Iplmc::register_fields |
  // +------------------------+
  void Iplmc::register_fields( Field& fields )
  {
    tracer::scope _( "Iplmc::regsiter_fields" );
    
    // allocate and initialize fields to be output
    m_colors.resize( m_subdomain.nhomes(), MAX_REAL );
    m_u.resize(      m_subdomain.nhomes(), MAX_REAL );
    m_v.resize(      m_subdomain.nhomes(), MAX_REAL );
    m_w.resize(      m_subdomain.nhomes(), MAX_REAL );
    m_f.resize(      m_subdomain.nhomes(), MAX_REAL );
    m_T.resize(      m_subdomain.nhomes(), MAX_REAL );
    m_O2.resize(     m_subdomain.nhomes(), MAX_REAL );
    m_CH4.resize(    m_subdomain.nhomes(), MAX_REAL );
    m_OH.resize(     m_subdomain.nhomes(), MAX_REAL );
    m_CO2.resize(    m_subdomain.nhomes(), MAX_REAL );
    m_NP.resize(     m_subdomain.nhomes(), MAX_REAL );
    m_cputime.resize(m_subdomain.nhomes(), MAX_REAL );

    // register with Fields class
    fields.add( m_colors, "colors" );
    fields.add( m_u,      "u"      );
    fields.add( m_v,      "v"      );
    fields.add( m_w,      "w"      );
    fields.add( m_f,      "f"      );
    fields.add( m_T,      "T"      );
    fields.add( m_O2,     "O2"     );
    fields.add( m_CH4,    "CH4"    );
    fields.add( m_OH,     "OH"     );
    fields.add( m_CO2,    "CO2"    );
    fields.add( m_NP,     "NP"     );
    fields.add( m_cputime, "cputime" );

  }

  // +----------------------+
  // | Iplmc::update_fields |
  // +----------------------+
  void Iplmc::update_fields()
  {
    tracer::scope _( "Iplmc::update_fields" );

    // update ensembled fields
    scalar_average();

    lid_t i = 0;
    BOOST_FOREACH( cell_mp::const_reference gidcell, m_cells )
    {
      const Cell& cell = gidcell.second;
      real_v vel;

      // get velocity
      m_tgv->uvw( m_time, m_domain.xyz( gidcell.first ), vel );

      m_colors[i]  = m_subdomain.my_rank();
      m_u[i]       = vel[XDIM];
      m_v[i]       = vel[YDIM];
      m_w[i]       = vel[ZDIM];
      m_f[i]       = cell.scalars.f;
      m_T[i]       = cell.scalars.phi.T();
      m_O2[i]      = cell.scalars.phi[ m_chem->index("O2") ];
      m_CH4[i]     = cell.scalars.phi[ m_chem->index("CH4") ];
      m_OH[i]      = cell.scalars.phi[ m_chem->index("OH") ];
      m_CO2[i]     = cell.scalars.phi[ m_chem->index("CO2") ];
      m_NP[i]      = cell.contained().size();
      m_cputime[i] = cell.ctr.cputime;

      i++;
    }

  }

  // +------------------+
  // | Iplmc::cell_wgts |
  // +------------------+
  void Iplmc::cell_wgts( std::vector< float >& wgts )
  {

    assert( wgts.size() == m_cells.size() );

    std::size_t i=0;
    BOOST_FOREACH( cell_mp::const_reference gidcell, m_cells )
      wgts[i++] = gidcell.second.ctr.cputime;
  }

} // namespace iplmcfd
