#ifndef IPLMCFD_PARTICLEMPITYPE_HPP_
#define IPLMCFD_PARTICLEMPITYPE_HPP_

#include "particle.hpp"
#include <mytools/mpi/opaque_objects.hpp>
#include <mpi.h>


namespace iplmcfd {

  //!
  //! \file particle_mpi_type.hpp
  //! \class ParticleMpiType
  //! \brief Partition the particle type into \e fixed and \e transient
  //!        sections.
  //!
  class ParticleMpiType
  {

  public:

    //! STATIC section contains the members of the particle struct whose 
    //! addresses are known once particle base address is known. The DYNAMIC 
    //! section contains runtime allocated member of the particle struct 
    //! (namely the phi_vector). 
    enum Sections { STATIC, DYNAMIC, NSECTION };
    //!
    typedef MPI_Aint section_addresses[NSECTION];
    //!
    typedef MPI_Datatype section_types[NSECTION];

  public:

    //!
    //!
    //!
    ParticleMpiType( const mix_t& mix )
    {
      // create a sample particle to work on
      Particle sample( mix );   

      // DYNAMIC part
      // create phi vector type (num_species+1 reals)
      MPI_Type_contiguous( 
        static_cast<int>(sample.phi.size()),  
        mpi_real, 
        types_[DYNAMIC].p()
        );
      MPI_Type_commit(types_[DYNAMIC].p());
      // store phi size for run time assertions, see get_part_addresses()
      phi_size = sample.phi.size();

      // FIXED part
      MPI_Type_contiguous( 
        Particle::fixed_num_reals(), 
        mpi_real, 
        types_[STATIC].p() 
        );
      MPI_Type_commit( types_[STATIC].p() );
    }

    //!
    //! \brief Returns for a given particle, the memory address of partitions. 
    //!
    void get_section_addresses( 
      const Particle& p, section_addresses addr ) const
    {
      BOOST_ASSERT( p.phi.size() == phi_size );
      MPI_Address( const_cast< real_t* >( p.first_fixed_member() ), 
        &addr[STATIC]  );
      MPI_Address( const_cast< real_t* >( &p.phi[0] ), &addr[DYNAMIC]  );
    }

    //!
    //! \brief Returns for a given particle, the MPI_TYPES of individual 
    //!        partitions.
    //!
    void get_section_types( section_types& dty ) const
    {
      for ( int i=0; i!=NSECTION; ++i )
        dty[i] = types_[i];
    }
    
  private:   
    
    //!
    mytools::mpi::Datatype types_[NSECTION];
    //!
    std::size_t phi_size;
    
  };

} // namespace iplmcfd

#endif // IPLMCFD_PARTICLEMPITYPE_HPP_
