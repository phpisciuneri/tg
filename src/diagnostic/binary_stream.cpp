#include "binary_stream.hpp"
#include "binary_cstdio.hpp"
#include "binary_fstream.hpp"



namespace iplmcfd {



//==============================================================================
// Factory functions -----------------------------------------------------------
//
//
binary_stream::ptr_t binary_stream::cstdio()
{
   return ptr_t(new binary_cstdio);
}


binary_stream::ptr_t binary_stream::cstdio(
   const std::string& name, 
   openmode op
   )
{
   return ptr_t(new binary_cstdio(name, op));
}

//
//
binary_stream::ptr_t binary_stream::fstream()
{
   return ptr_t(new binary_fstream);
}

binary_stream::ptr_t binary_stream::fstream(
   const std::string& name, 
   openmode op
   )
{
   return ptr_t(new binary_fstream(name, op));
}

//
//
#ifndef LMCFD_INEFFICIENT_FSTREAM
#define use_impl fstream
#else 
#define use_impl cstdio
#endif

binary_stream::ptr_t binary_stream::optimal() { return use_impl(); }

binary_stream::ptr_t binary_stream::optimal(
   const std::string& name, 
   openmode op
   )
{
   return use_impl(name,op);
}

#undef use_impl



} // namespace iplmcfd


