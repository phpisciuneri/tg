//#include "../iplmc.hpp"
//
//#include <limits>
//#include <cmath>
//#include <boost/foreach.hpp>
//#include <mytools/numeric/indsub.hpp>
//
//#include "../stats.hpp"
//#include "../cpublock.hpp"
//#include "../iplmcfd.hpp"
//#include "../qref.hpp"
//
//namespace iplmcfd {
//
//
////  +------------+
////  |  ENSEMBLE  |
////  +------------+
//
//void Iplmc::ensemble()
//{
//
//	const Chemistry::mix_t& mix = Chemistry::instance().mixture();
//
//   for ( size_t i=0; i<m_subdomain.nhomes(); ++i )
//   {
//
//		Intersolver::mc2fd_t& exports = m_intersolver.mc2fd[i];
//
//		// zero the exports before we ensemble average
//		m_stats[i] = 0;
//		exports = 0;
//
//		cpublock cpub( m_cells[i].ctr.cputime, m_cells[i].contained().size() );
//
//		BOOST_FOREACH( particle& p , m_cells[i].contained() )
//		{
//			real R   = mix.R(p.phi.Y());
//			if ( m_integ == SZ ) R = qref::instance()[qref::RT] / qref::instance()[qref::T];
//
//			real RT  = p.phi.T()*R;
//			real rho = // TODO fix this pressure being 1 atm arbitrarily
//					Chemistry::mix_t::atm2dynescc()/RT;
//
//			real hdot = 0.;
//			// overload with rho for the taylor-green case
//      exports[Intersolver::SRT]  += rho*p.mass;
//			exports[Intersolver::MASS] += p.mass;
//			exports[Intersolver::MU]   += // this is a nasty implementation here for now: MU is set to ensemble
//					p.phi.T()*p.mass;      // averaged temperature. Later inside compress::normalize_imports a
//			                               // simple formula will be used to set mu from this average temperature.
//
//			// compute other ensemble values of interest: means, RMS, etc.
//			stats::ensemble( m_stats[i], p, RT, rho, hdot );
//
//		}
//
//	}
//
//}
//
//} // namespace iplmcfd
