#ifndef IPLMCFD_CPUBLOCK_HPP_
#define IPLMCFD_CPUBLOCK_HPP_


#include <cstdlib> // for size_t

namespace iplmcfd {



class cpublock
{
public:
   // t is the time element that will be appended the blocks time. it is normally
   // a real quantity  but for implementation purposes I needed to make it size_t (see cell:counter)
   // and therefore, 
   // mult is what the actual CPU time will be multiplied with before truncating it into size_t.
   //
   // n is the number of units to amortize the block time (for example, the block that needs to be timed
   // is inside a loop and each iteration is pretty much homogeneous. so for efficiency the cpublock is
   // put outside the loop, but will be amortized by the loop iteration count. 
   //
   cpublock (size_t& t, size_t n, size_t mult = 1000000);
   ~cpublock();
private:
   size_t* m_t;
   size_t m_n;
   double m_bt; // block time
   size_t m_mult; 
};


} // namespace iplmcfd 


#endif // IPLMCFD_CPUBLOCK_HPP_
