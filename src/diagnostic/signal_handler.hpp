#ifndef IPLMCFD_SIGNAL_HANDLER_HPP_
#define IPLMCFD_SIGNAL_HANDLER_HPP_


namespace iplmcfd {

class signal_handler
{
   public:
      static bool interrupted() 
        #ifdef LMCFD_ENABLE_SIGNAL_HANDLING
         ;
        #else
        { return false; }
        #endif
};

} // namespace iplmcfd


#endif // IPLMCFD_SIGNAL_HANDLER_HPP_
