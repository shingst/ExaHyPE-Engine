// @tood Remove
#define SHM_DEBUG_PREFIX "SHM\t"

#if SHM_INVADE_DEBUG>=1
 #ifndef TBB_USE_THREADING_TOOLS
  #warning We recommend to compile with -DTBB_USE_THREADING_TOOLS if SHMInvade is compile in debug mode
 #endif
 #ifndef TBB_USE_ASSERT
  #warning We recommend to compile with -DTBB_USE_ASSERT if SHMInvade is compile in debug mode
 #endif
#endif

#define SHM_MIN_SLEEP 1
#define SHM_MAX_SLEEP 30
