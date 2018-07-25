/**
Copyright (C) 2018, Martin Schreiber and Tobias Weinzierl

All rights reserved.

Open Source License
===================

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.

2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.

3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

4. Scientific and commercial work using the software should cite the authors' corresponding papers.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef _SHMINVADE_SHMMACROS_H_
#define _SHMINVADE_SHMMACROS_H_

#define SHM_DEBUG_PREFIX    "SHM"
#define SHM_DEBUG_SEPARATOR "\t"

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

#include <string>

namespace shminvade {
  std::string getSHMDebugPrefix();
}

#endif
