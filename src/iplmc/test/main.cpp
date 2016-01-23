#include <boost/cstdlib.hpp>
#include <boost/bind.hpp>
#include <boost/test/included/prg_exec_monitor.hpp>
#include "vectorlist.hpp"
#include "sortedvector.hpp"
#include "randomshuffler.hpp"




int cpp_main(int,char*[])
{
   size_t npg   = 10;
   size_t nCell = 100000;
   size_t nParticles = nCell * npg;
   
   {
      vectorlist vl(nParticles, nCell);
      { boost::progress_timer tm;
      randomshuffler shuffler(nCell);
      vl.for_each_p( boost::bind(boost::ref(shuffler),_1) );
      std::cout << shuffler.total_moves() << std::endl;
      vl.redistribute();
      }
   }
   
   {
      sortedvector vl(nParticles, nCell);
      { boost::progress_timer tm; 
      randomshuffler shuffler(nCell);
      vl.for_each_p( boost::bind(boost::ref(shuffler),_1) );
      std::cout << shuffler.total_moves() << std::endl;
      vl.redistribute();
      }
   }   

   return boost::exit_success;

}
