#include "enzo_defines.hpp"

#define DEBUG_MESSAGE            CALL f_warning(__FILE__,__LINE__)
#define ERROR_MESSAGE            CALL f_error(__FILE__,__LINE__)
#define WARNING_MESSAGE          CALL f_warning(__FILE__,__LINE__)

#ifdef CONFIG_PRECISION_SINGLE
#  define ENZO_REAL real*4
#endif

#ifdef CONFIG_PRECISION_DOUBLE
#  define ENZO_REAL real*8
#endif

#ifdef CONFIG_PRECISION_QUAD
#  define ENZO_REAL real*16
#endif