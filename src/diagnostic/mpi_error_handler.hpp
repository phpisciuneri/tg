#ifndef IPLMCFD_MPI_ERROR_HANDLER_HPP_
#define IPLMCFD_MPI_ERROR_HANDLER_HPP_

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#include <string>
#include <mpi.h>
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
namespace iplmcfd {
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


/** \brief Error handler and more for mpi calls
 *
 * The idea is to make MPI errors more informative by specifying at which file 
 * at what line the erroneous call is. Implementation is really simple: We have
 * this singleton error handler object that stores file and line location. When
 * an error is triggered the handler will print out the file/line info. Then
 * will call default mpi error handler. 
 */

class mpi_error_handler 
{

   //===========================================================================
   // SINGLETON
   //===========================================================================
private: 
   mpi_error_handler();
   ~mpi_error_handler();
public:
   static mpi_error_handler& instance() 
   {
      static mpi_error_handler singleton;
      return singleton;
   }  


   //===========================================================================
   // INTERFACE ----------------------------------------------------------------
public:
   void set(const std::string& filename, size_t line)
   {
      filename_ = filename;
      line_ = line;     
   }
   const std::string& filename() const { return filename_; }
   size_t line() const { return line_; }
   //===========================================================================



   //===========================================================================
   // DATA MEMBERS -------------------------------------------------------------
private:
   std::string filename_;
   size_t line_;
   MPI_Errhandler handler_handle_;
   //===========================================================================
};





//==============================================================================
// MACRO Interface for automated error handler reset
//==============================================================================
#define LMCFD_MPI_REPORT_ERROR( statement ) \
   ::lmcfd::mpi_error_handler::instance().set(__FILE__, __LINE__); \
   statement;





//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
} // namespace iplmcfd
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~



#endif // IPLMCFD_MPI_ERROR_HANDLER_HPP_
