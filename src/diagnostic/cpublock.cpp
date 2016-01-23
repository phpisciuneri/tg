

#include "cpublock.hpp"
#include <mpi.h>


namespace iplmcfd {



   cpublock::cpublock(size_t& t, size_t n, size_t mult)
   {
      m_t  = &t;
      m_n  = n;
      m_bt = MPI_Wtime();
      m_mult = mult;
   }


   cpublock::~cpublock()
   {
      m_bt = (MPI_Wtime() - m_bt)/m_n;
      *m_t += (size_t)(m_bt*m_mult);
   }



} // namespace iplmcfd 
