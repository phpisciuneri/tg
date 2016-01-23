#ifndef BIN_HPP_
#define BIN_HPP_


#include <list>



template < class Element >
class bin
{



   //===========================================================================
   //  TYPEDEFS ----------------------------------------------------------------
public:
   typedef std::list< Element >   array;
   typedef typename array::iterator       iterator;
   typedef typename array::const_iterator const_iterator;
   //===========================================================================






   //===========================================================================
   // OWN    Assume ownership of all incoming elements
   //===========================================================================
   void own()
   {
      contained.splice(contained.begin(), incoming);
   }




   //===========================================================================
   // TRANSFER  Transfer an element from one bin to another
   //===========================================================================
   friend void transfer( 
      iterator p, 
      bin& from_host,  
      bin& to_host
      )
   {
      to_host.incoming.splice(
         to_host.incoming.begin(), from_host.contained, p 
         );
   }


   //===========================================================================
   // SIZE   Number of contained elements
   //===========================================================================
   size_t size() const { return contained.size(); }






   //===========================================================================
   // DATA MEMBERS -------------------------------------------------------------
public:
   array incoming;
   array contained;
   //===========================================================================



};

  

#endif //BIN_HPP_
