#ifndef IPLMCFD_VERSION_HPP_
#define IPLMCFD_VERSION_HPP_


#include <string>

namespace iplmcfd {


   class version
   {
      // singleton
   private:
      version(); 
      ~version(){};
   public:
      static const version& instance() { 
         static version singleton;
         return singleton;
      }
            
      // methods
   public:
      static const int LEN = 16;

      std::string str() const { return m_str; }
      int         num() const { return m_num; }


   private:
      char m_str[LEN]; ///> fixed size for dumping easily
      int m_num;
   };


} // namespace iplmcfd


#endif // IPLMCFD_VERSION_HPP_
