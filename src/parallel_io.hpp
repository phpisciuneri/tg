#ifndef IPLMCFD_PARALLEL_IO_HPP_
#define IPLMCFD_PARALLEL_IO_HPP_

#include "domain.hpp"
#include "subdomain.hpp"
#include "field.hpp"

#include <map>
#include <boost/scoped_ptr.hpp>

namespace iplmcfd {

  //!
  //! \class Parallel_io
  //! \brief
  //!
  class Parallel_io
  {
  public:

    //!
    //!
    //!
    Parallel_io( const Domain& d, const Subdomain& s, const Field& fields );

    //!
    //!
    //!
    ~Parallel_io();

    //!
    //!
    //!
    void output_paraview( const Field& fields, int iter, int iter_stats = 0);
    
    //!
    //! complete non-blocking collectives and close previous output
    //!
    void close();

  private:

    //!
    struct Impl;
    //!
    boost::scoped_ptr<Impl> p;

  };

} // namespace iplmcfd

#endif // IPLMCFD_PARALLEL_IO_HPP_
