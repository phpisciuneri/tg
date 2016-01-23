#ifndef IPLMCFD_ACTIVITY_LOG_HPP_
#define IPLMCFD_ACTIVITY_LOG_HPP_

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include <fstream>
#include <sstream>
#include <string>
#include <boost/format.hpp>
#include <boost/array.hpp>

#include "../defs.hpp"
#include "../subdomain.hpp"
#include "../particle.hpp"

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
namespace iplmcfd {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



//! Singleton for general code-wide activity log


class activity_log : boost::noncopyable
{

	//===========================================================================
	// TYPEDEFS/ENUMS -----------------------------------------------------------
public: 
	enum {
		TOTAL, ERASED, CREATED, INCAME, OUTGONE, ROGUE, SHIFTED, MIGRATED,
		EMPTY_CELL, BULLIED_CELL, UNDERPOP_CELL, CROWDED_CELL, N_ACTIVITY
	};

	enum LogLevel { LOG_NOT, LOG_SUM, LOG_ALL, LOG_NLEVEL };
	//===========================================================================



	//===========================================================================
	// SINGLETON ----------------------------------------------------------------
private:
	activity_log();
	~activity_log(){end_iteration();}
	//===========================================================================


	//===========================================================================
	// INTERFACE ----------------------------------------------------------------
public:
	//___________________________________________________________________________
	//! Access the singleton object
	static activity_log& instance()
	{
		static activity_log log;
		return log;
	}


	//___________________________________________________________________________
	//! Attach upcoming activity this particle
	void attach_particle( const Particle& )
	{
		m_counter[TOTAL]++;
	}

	//___________________________________________________________________________
	//! Various particle activities
	//@{
	void shifted_particle () { ++m_counter[SHIFTED]; }
	void migrated_particle() { ++m_counter[MIGRATED];}
	void created_particle () { ++m_counter[CREATED]; }
	void erased_particle  () { ++m_counter[ERASED];  }
	void incame_particle  () { ++m_counter[INCAME];  }
	void outgone_particle () { ++m_counter[OUTGONE]; }
	void rogue_particle   (const int_v& shift) {
		++m_counter[ROGUE];
		if (m_loglevel < LOG_ALL ) return;
		m_logfile << (m_frogue % m_current_cell_id % shift[0] % shift[1] % shift[2]);
	}
	//@}


	//___________________________________________________________________________
	//! Attach upcoming activity this cell
	void attach_cell( gid_t cell_id ) {  m_current_cell_id = cell_id; }

	//___________________________________________________________________________
	//! Various cell activities
	//@{
	void empty_cell(size_t nAdd)
	{
		++m_counter[EMPTY_CELL];
		m_counter[CREATED] += nAdd;
		if (m_loglevel < LOG_ALL ) return;
		m_logfile << m_fempty % m_current_cell_id;
	}
	void bullied_cell(size_t nHeavy, size_t nCloned)
	{
		++m_counter[BULLIED_CELL];
		m_counter[CREATED] += nCloned;
		if (m_loglevel < LOG_ALL ) return;
		m_logfile << m_fbullied % m_current_cell_id %  nHeavy % nCloned;
	}
	void crowded_cell(size_t nMerged)
	{
		++m_counter[CROWDED_CELL];
		m_counter[ERASED] += nMerged;
		if (m_loglevel < LOG_ALL ) return;
		m_logfile << m_fcrowded % m_current_cell_id %  nMerged;
	}
	void underpop_cell(size_t nAdded)
	{
		++m_counter[UNDERPOP_CELL];
		m_counter[CREATED] += nAdded;
		if (m_loglevel < LOG_ALL ) return;
		m_logfile << m_funderpop % m_current_cell_id %  nAdded;
	}
	//@}



	//___________________________________________________________________________
	void begin_iteration(size_t iter);
	void end_iteration();


	//___________________________________________________________________________
	//! Set log detail level
	void set_log_level( unsigned level  )
	{
		if ( level >= LOG_NLEVEL ) throw(std::invalid_argument("Invalid carlo_log_level"));
		m_loglevel = (LogLevel) level;
	}

	//===========================================================================






	//===========================================================================
	// DATA MEMBERS -------------------------------------------------------------
private:
	std::ofstream m_logfile;
	std::ofstream m_sumfile;
	boost::array<size_t,N_ACTIVITY> m_counter;
	boost::format m_frogue;
	boost::format m_fempty;
	boost::format m_fbullied;
	boost::format m_fcrowded;
	boost::format m_funderpop;
	gid_t m_current_cell_id;
	LogLevel m_loglevel;
	bool m_flush;
	size_t m_iter;
	//===========================================================================
};







//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
} // namespace iplmcfd
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#endif // IPLMCFD_ACTIVITY_LOG_HPP_
