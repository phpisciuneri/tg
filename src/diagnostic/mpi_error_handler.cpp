//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include "mpi_error_handler.hpp"
#include <boost/format.hpp>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

extern "C" void lmcfd_mpi_error_handler_fn( MPI_Comm* comm, int* errcode, ... )
{
   char error_string[MPI_MAX_ERROR_STRING];
   int len;
   MPI_Error_string( *errcode, error_string, &len );
   error_string[len] = '\0';
   throw std::logic_error((
      boost::format("Error in MPI call from %1%:%2%\n  %3%")
      % ::iplmcfd::mpi_error_handler::instance().filename()
      % ::iplmcfd::mpi_error_handler::instance().line()
      % error_string
      ).str());   
}



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
namespace iplmcfd {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



mpi_error_handler::mpi_error_handler()
{
   MPI_Comm_create_errhandler( lmcfd_mpi_error_handler_fn, &handler_handle_ );
   MPI_Comm_set_errhandler(MPI_COMM_WORLD, handler_handle_);
}

mpi_error_handler::~mpi_error_handler()
{
   MPI_Errhandler_free( &handler_handle_ );
}





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
} // namespace iplmcfd
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

