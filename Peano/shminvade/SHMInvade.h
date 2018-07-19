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
#ifndef _SHMINVADE_SHMINVADE_H_
#define _SHMINVADE_SHMINVADE_H_


#include <sys/types.h>
#include <set>


namespace shminvade {
  class SHMInvade;
}


/**
 * Fundamental invasion object. If you create it, you have tell the object
 * how many cores you'd like to invade. You may also try to invade all cores.
 * Then, the object asks the strategy to identify whether there are cores
 * available and, if successful, memorises those.
 *
 * If you call retreat() or if the object is destroyed, it tells the
 * SHMController that all booked cores are not required anymore.
 */
class shminvade::SHMInvade {
  private:
    std::set<int> _occupiedCores;
  public:
    static constexpr int MaxCores = -1;

    SHMInvade(int cores);
    ~SHMInvade();
    void retreat();
};

#endif
