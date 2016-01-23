#include <boost/test/unit_test.hpp>
#include "../domain.hpp"


BOOST_AUTO_TEST_SUITE( domain )




//  +-------------------------------------------+
//  |  Custom predicate for Vector comparisons  |
//  +-------------------------------------------+
// 
// See: 
// http://www.boost.org/doc/libs/1_46_0/libs/test/doc/html/utf/testing-tools/custom-predicate.html
//
template<class T1, class T2>
boost::test_tools::predicate_result
all_equal( const T1& a1, const T2& a2 )
{
   boost::test_tools::predicate_result res = all_elements(a1 == a2 );
   if ( !res )  
      res.message() << a1 << " != " <<  a2;
   return res;
}




//  +---------------------------+
//  |  Check Methods of Domain  |
//  +---------------------------+

BOOST_AUTO_TEST_CASE( methods )
{
   using namespace iplmcfd;

   // Test stuff

   shape_vector shape(3,4,5);
   r3d dx( 1,1,1 );
   r3d xmin( 0,0,0 );
   bool_vector periods( false, false, true );
   size_t npoints = shape[0]*shape[1]*shape[2];

   //////////////////////////////////////////////////////////////////////////
   // CTOR
   BOOST_REQUIRE_NO_THROW( Domain dom(shape,dx,xmin,periods) );   

   // Now create to use: 
   Domain dom(shape,dx,xmin,periods);

   //////////////////////////////////////////////////////////////////////////
   BOOST_CHECK_EQUAL( dom.npoints(), npoints  );
   BOOST_CHECK( all_equal(dom.nx(), shape ) );


   //////////////////////////////////////////////////////////////////////////
   //  reflect
   {      
      boost::array<i3d,4> s  = { i3d(1,3,5), i3d(1,3,7), i3d(1,3,15), i3d(1,3,-1) };
      boost::array<i3d,4> rs = { i3d(1,3,0), i3d(1,3,2), i3d(1,3,0) , i3d(1,3,4)  };
      for (size_t i=0; i<s.size(); i++) 
         BOOST_CHECK( all_equal( dom.reflect(s[i]), rs[i]) );
   }


   //////////////////////////////////////////////////////////////////////////
   // subs -> index
   {   
      i3d subs[] = { i3d(0,0,0), i3d(2,0,0), i3d(0,1,0), i3d(0,1,1), i3d(0,0,5) };
      gid inds[] = { 0         , 2         , 3         , 15        , 0          };

      for( int i=0; i<sizeof(inds)/sizeof(gid); ++i )  {
         BOOST_CHECK_EQUAL( dom.index(subs[i]),  inds[i] );
         BOOST_CHECK_EQUAL( dom.sub2ind(subs[i]),  inds[i] );
      }
   }


   //////////////////////////////////////////////////////////////////////////
   // index -> subs
   {   
      i3d subs[] = { i3d(0,0,0), i3d(2,0,0), i3d(0,1,0), i3d(0,1,1), i3d(0,0,3) };
      gid inds[] = { 0         , 2         , 3         , 15        , 36         };

      for( int i=0; i<sizeof(inds)/sizeof(gid); ++i )  {
         BOOST_CHECK( all_equal(dom.subs(inds[i]), subs[i]) );
         BOOST_CHECK( all_equal(dom.ind2sub(inds[i]), subs[i]));
      }
   }

   
   //////////////////////////////////////////////////////////////////////////
   // out_of_bounds
   BOOST_CHECK( dom.out_of_bounds(npoints)  );
   BOOST_CHECK( dom.out_of_bounds(i3d(0,-1,0))  );
   BOOST_CHECK( dom.out_of_bounds(i3d(3,0,0) )  );
   BOOST_CHECK( dom.out_of_bounds(i3d(0,4,0) )  );
   BOOST_CHECK( dom.out_of_bounds(i3d(0,0,5) )  );
   BOOST_CHECK( dom.out_of_bounds(i3d(3,0,5) )  );
   
   BOOST_CHECK( !dom.out_of_bounds_periodic(i3d(0,0,5) )  );


   //////////////////////////////////////////////////////////////////////////
   // len, x, xmin, xmax ...


   i3d sub(1,2,0);
   gid I = sub[0] + sub[1]*shape[0] + sub[2]*shape[0]*shape[1];
   r3d x ( xmin + dx*sub );

   for (int i=0; i<NDIM; ++i) {
      BOOST_CHECK_EQUAL( dom.len()[i]  ,  dx[i]*(shape[i]-1) );
      BOOST_CHECK_EQUAL( dom.xmin()[i] ,  xmin[i] );
      BOOST_CHECK_EQUAL( dom.xmax()[i] ,  xmin[i]+dx[i]*(shape[i]-1) );
      BOOST_CHECK_EQUAL( dom.x(I)[i]   ,  x[i] );
      BOOST_CHECK_EQUAL( dom.x(sub)[i] ,  x[i] );

      BOOST_CHECK_EQUAL( dom.c_len()[i]  ,  dx[i]*shape[i] );
      BOOST_CHECK_EQUAL( dom.c_xmin()[i] ,  xmin[i]-dx[i]/2 );
      BOOST_CHECK_EQUAL( dom.c_xmax()[i] ,  xmin[i]-dx[i]/2 + dx[i]*shape[i] );

      BOOST_CHECK_EQUAL( dom.c_x(I)[i]  ,  x[i]-dx[i]/2 );
      BOOST_CHECK_EQUAL( dom.c_x(sub)[i],  x[i]-dx[i]/2 );

      BOOST_CHECK_EQUAL( dom.odx()[i],  1/dx[i] );
   }
   

}

BOOST_AUTO_TEST_SUITE_END()
