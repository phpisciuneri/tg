#ifndef IPLMCFD_BINARY_STREAM_HPP_
#define IPLMCFD_BINARY_STREAM_HPP_


#include <boost/shared_ptr.hpp>
#include <ios>
#include <string>

namespace iplmcfd {


   /** I failed to (summon courage and convince myself the usefulness of) 
    *  properly derive a streambuf using FILE* as an underlying buffer. 
    *  This is just a stupid class used by save/load routines. 
    */
   class binary_stream
   {
      
      // factory
   public:
      typedef boost::shared_ptr<binary_stream> ptr_t;
      typedef std::ios::openmode openmode;
      typedef std::streamoff     streamoff;
      typedef std::streampos     streampos;
      typedef std::ios::seekdir  seekdir;
      
      static ptr_t fstream(const std::string&, openmode);  // return fstream implementation
      static ptr_t fstream();                              // return fstream implementation

      static ptr_t cstdio (const std::string&, openmode);   // return cstdio  FILE* implementation
      static ptr_t cstdio ();
      
      static ptr_t optimal(const std::string&, openmode);  // return optimal for a given system
      static ptr_t optimal();

      // interface 
   public: 
      virtual bool is_open() const = 0;
      virtual void open(const std::string&, openmode)  = 0;
      virtual void close() = 0;
      virtual binary_stream& write( const char* buffer, size_t sz ) = 0;
      virtual binary_stream& read( char* buffer, size_t sz ) = 0;
      virtual binary_stream& seek( streamoff, seekdir ) = 0;
      virtual binary_stream& seek( streampos ) = 0;
      virtual streampos tell() = 0;
      virtual bool fail() const = 0;
      virtual ~binary_stream(){};

   public:
      operator bool() const { return !fail(); }
   };


} // namespace iplmcfd

#endif // IPLMCFD_BINARY_STREAM_HPP_
