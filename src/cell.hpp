#ifndef IPLMCFD_CELL_HPP_
#define IPLMCFD_CELL_HPP_

#include <list>
#include <boost/array.hpp>
#include <mytools/mpi/opaque_objects.hpp>

#include "defs.hpp"
#include "domain.hpp"
#include "particle.hpp"
#include "scalar_averages.hpp"

namespace iplmcfd 
{

  //!
  //! \file cell.hpp
  //! \class Cell
  //! \brief
  //!
  class Cell
  {

  public:

    //!
    typedef std::list< Particle > particle_lst;
    //!
    typedef particle_lst::iterator particle_it;
    //!
    typedef particle_lst::const_iterator particle_cit;
    //!
    enum {NCORNER=8, NSIDE=6, NEDGE=12, NNEIGHBOR=26};
    
    //!
    struct Counters
    {
      typedef size_t counter_t;

      static MPI_Datatype mpi_type();
      void reset();

      //! number of counters
      enum { CLONED, MERGED, WENTOUT, CAMEIN, IN, CHEMLOAD, CPUTIME, NCOUNTS }; 
      
      Counters() 
        : cloned(0), 
        merged(0), 
        wentout(0), 
        camein(0), 
        in(0), 
        chemload(0), 
        cputime(0) 
      {}

      size_t cloned;   //!<
      size_t merged;   //!<
      size_t wentout;  //!<
      size_t camein;   //!<
      size_t in;       //!<
      size_t chemload; //!<
      size_t cputime;  //!< in milli seconds

    };

  public:
    
    //!
    //! \brief
    //!
    Cell( const mix_t& mix ) : scalars( mix ) {}

    //!
    //! \brief Get list of particles contained by the cell
    //!
    const particle_lst& contained() const { return m_contained; }

    //!
    //! \brief Modify list of particles contained by the cell
    //!
    particle_lst& contained() { return m_contained; }

    //!
    //! \brief Get list of particles entering the cell
    //!
    const particle_lst& incoming() const { return m_incoming; }

    //!
    //! \brief Modify list of particles entering cell
    //!
    particle_lst& incoming() { return m_incoming; }
    
    //!
    //! \brief Tranfser ownership of particles
    //!
    //! Claim new particles from incoming by splicing all in incoming into 
    //!   contained
    //!
    void own_incoming();
    
    //!
    //! Enforce particle cloning and maximum number density control
    //!
    void adjust_population(const size_t min_ppc, const size_t mean_ppc, const size_t max_ppc);

    //!
    //! \brief to MPI message (to be used as hindexed type of <tt>mpi_real</tt>s)
    //!
    //! \todo It'd be nice to have a separate class for MPI messaging, that
    //!       containes these vectors and when asked, create the MPI datatype. The
    //!       applications (such as cell) would know about it and implement
    //!       methods to fill or append or whatnot to the message.
    //!
    void append_to_message_real(std::vector<MPI_Aint>&, std::vector<int>&);
    
    //!
    //! Append cell count into MPI message
    //!
    void append_to_message_counts(
      std::vector<MPI_Aint>& chunk_address,
      std::vector<int>& chunk_length
      );


  public:
    
    //!
    Counters ctr;
    //!
    ScalarAverages scalars;
    //!
    real_t mass;
  
  private:
    
    //!
    particle_lst m_contained;
    //!
    particle_lst m_incoming;  

  };
  
  //! +--------------------+
  //! | Cell::own_incoming |
  //! +--------------------+
  inline void Cell::own_incoming()
  {
    ctr.camein += m_incoming.size();
    m_contained.splice( m_contained.begin(), m_incoming );
  }

  //!
  //! Transfer a single particle from one cell to another
  //!
  inline void transfer( Cell::particle_it p, Cell& old_host, Cell& new_host )
  {
    new_host.incoming().splice(
      new_host.incoming().begin(), old_host.contained(), p
      );
  }

  //! +-----------------------+
  //! | Cell::Counters::reset |
  //! +-----------------------+
  inline void Cell::Counters::reset() 
  {
    cloned   = 0;
    merged   = 0;
    wentout  = 0;
    camein   = 0;
    in       = 0;
    chemload = 0;
    cputime  = 0;
  }

  //! +--------------------------+
  //! | Cell::Counters::mpi_type |
  //! +--------------------------+
  inline MPI_Datatype Cell::Counters::mpi_type()
  {
    static mytools::mpi::Datatype type;
    static bool initialized = false;
    if (!initialized)
    {
      MPI_Type_contiguous( NCOUNTS,
        mytools::mpi::get_datatype( counter_t() ),
        type.p()
        );
      MPI_Type_commit( type.p() );
      initialized = true;
    }
    return type;
  }

} // namespace iplmcfd

#endif // IPLMCFD_CELL_HPP_