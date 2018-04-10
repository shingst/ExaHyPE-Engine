# SHMInvade library #


## Build SHMInvade ##

We propose to build two versions of the library. This mirrors the approach
realised by Intel with their TBB:

export TBB_INC=-I/opt/tbb/include

g++ -std=c++14 -fPIC -O2 $TBB_INC -DSHM_INVADE_DEBUG=1 -DTBB_USE_ASSERT -DTBB_USE_THREADING_TOOLS *.cpp -shared -o libshminvade_debug.so
g++ -std=c++14 -fPIC -O3 $TBB_INC -DNDEBUG -DSHM_INVADE_DEBUG=0 *.cpp -shared -o libshminvade.so
rm *.out


Please note that you should translate with 

-DTBB_USE_THREADING_TOOLS

If you want to use the libs with the Intel threading tools.



Aufgreifen im Paper:

Meinen eigenen Kontext-Kommentar aufgreifen
Muteces
Low priority for context
Manuel deletion within templates
