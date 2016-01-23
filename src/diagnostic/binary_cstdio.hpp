#ifndef IPLMCFD_BINARY_CSTDIO_HPP_
#define IPLMCFD_BINARY_CSTDIO_HPP_


#include "binary_stream.hpp"
#include <cstdio>


namespace iplmcfd {




//==============================================================================
// CSTDIO  ---------------------------------------------------------------------
class binary_cstdio : public binary_stream
{
public:
   binary_cstdio() : m_fp(NULL), m_failed(true) {}
   binary_cstdio( const std::string&, openmode);
   ~binary_cstdio() { close(); }
   // interface 
   virtual bool is_open() const;
   virtual void open(const std::string&, openmode);
   virtual void close();
   virtual binary_stream& write( const char* buffer, size_t sz );
   virtual binary_stream& read( char* buffer, size_t sz );
   virtual binary_stream& seek( streamoff, seekdir );
   virtual binary_stream& seek( streampos );
   virtual streampos tell();
   virtual bool fail() const;
private:
   std::FILE* m_fp;
   bool m_failed;
};






//==============================================================================
// Implementation --------------------------------------------------------------

//
inline binary_cstdio::binary_cstdio( 
                                    const std::string& name, 
                                    openmode mod
                                    ) 
                                    : m_fp(NULL)
                                    , m_failed(false)
{
   open(name, mod);
}



//
inline bool binary_cstdio::is_open() const
{
   return m_fp!=NULL;
}

//
void binary_cstdio::open(const std::string& name, openmode mod)
{
   close(); // current will be closed even if this open fails
   std::string cmode; 
   using std::ios; 
#define isset(x) ((mod&(x))==(x))
   if (isset(ios::in))
      if (isset(ios::out))
      {
         if     (isset(ios::app))   cmode = "a+";
         else if(isset(ios::trunc)) cmode = "w+";
         else cmode = "r+";
      }
      else cmode = "r";
   else if (isset(ios::out)) cmode = "w";
   else if (isset(ios::app)) cmode = "a";
#undef isset
   else 
   {
      m_failed = true;
      m_fp = NULL;
      return;
   }
   //
   // binary!
   cmode += "b";
   m_fp = std::fopen(name.c_str(), cmode.c_str());
   m_failed = !is_open();
}


//
inline void binary_cstdio::close()
{
   if (m_fp) 
   {
      m_failed = std::fclose(m_fp)!=0;
      m_fp = NULL;
   }
}



//
inline binary_stream& binary_cstdio::write( const char* buffer, size_t sz )
{
   m_failed = std::fwrite(buffer, sizeof(char), sz, m_fp) != sz;
   return *this;
}

//
inline binary_stream& binary_cstdio::read( char* buffer, size_t sz )
{
   m_failed = std::fread(buffer, sizeof(char), sz, m_fp) != sz;
   return *this;
}

//
inline binary_stream& binary_cstdio::seek( 
   streamoff offset, 
   seekdir  sk
   )
{
   if (sk==std::ios::beg) m_failed = std::fseek(m_fp, offset, SEEK_SET)!=0;
   if (sk==std::ios::cur) m_failed = std::fseek(m_fp, offset, SEEK_CUR)!=0;
   if (sk==std::ios::end) m_failed = std::fseek(m_fp, offset, SEEK_END)!=0;
   return *this;
}

//
inline binary_stream& binary_cstdio::seek( 
   streampos pos
   )
{
   m_failed = std::fseek(m_fp, pos, SEEK_SET)!=0;
   return *this;
}

//
inline binary_stream::streampos binary_cstdio::tell() 
{
   return std::ftell( m_fp );
}


//
inline bool binary_cstdio::fail() const
{
#ifndef BOOST_MSVC 
   using std::ferror; // MSVC somehow managed not to wrap ferror in std::
#endif
   return m_failed && ferror(m_fp);
}



} // namespace iplmcfd

#endif // IPLMCFD_BINARY_CSTDIO_HPP_
