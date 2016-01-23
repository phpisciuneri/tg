#ifndef IPLMCFD_BINARY_FSTREAM_HPP_
#define IPLMCFD_BINARY_FSTREAM_HPP_


#include "binary_stream.hpp"
#include <fstream>


namespace iplmcfd {




//==============================================================================
// CSTDIO  ---------------------------------------------------------------------
class binary_fstream : public binary_stream
{
public:
   // interface 
   binary_fstream() {}
   binary_fstream( const std::string& name, openmode mod) : m_f(name.c_str(),mod|std::ios::binary) {m_f.exceptions(std::ios::eofbit|std::ios::failbit|std::ios::badbit);}
   ~binary_fstream() {}
   // interface 
   virtual bool is_open() const { return const_cast<std::fstream&>(m_f).is_open(); }
   virtual void open(const std::string& name, openmode mod) { m_f.open(name.c_str(), mod|std::ios::binary); }
   virtual void close() { m_f.close(); }
   virtual binary_stream& write( const char* buffer, size_t sz );
   virtual binary_stream& read( char* buffer, size_t sz );
   virtual binary_stream& seek( streamoff, seekdir );
   virtual binary_stream& seek( streampos );
   virtual streampos tell() { return m_f.tellg(); }
   virtual bool fail() const { return m_f.fail(); }

private:
   std::fstream m_f;
};






//==============================================================================
// Implementation --------------------------------------------------------------

//
inline binary_stream& binary_fstream::write( const char* buffer, size_t sz )
{
   m_f.write(buffer, sz);
   return *this;
}

//
inline binary_stream& binary_fstream::read( char* buffer, size_t sz )
{
   m_f.read(buffer, sz);
   return *this;
}

//
inline binary_stream& binary_fstream::seek( 
   streamoff offset, 
   seekdir  seekdir
   )
{
   m_f.seekg(offset, seekdir);
   return *this;
}

//
inline binary_stream& binary_fstream::seek(streampos p) 
{
   m_f.seekg(p);
   return *this;
}



} // namespace iplmcfd

#endif // IPLMCFD_BINARY_FSTREAM_HPP_
