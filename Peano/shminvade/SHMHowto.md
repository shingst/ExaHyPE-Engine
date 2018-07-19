# SHMInvade library #


### Build SHMInvade ###

We propose to build two versions of the library. This mirrors the approach
realised by Intel with their TBB:

export TBB_INC=-I/opt/tbb/include

g++ -std=c++14 -fPIC -O2 $TBB_INC -DSHM_INVADE_DEBUG=1 -DTBB_USE_ASSERT -DTBB_USE_THREADING_TOOLS *.cpp -shared -o libshminvade_debug.so
g++ -std=c++14 -fPIC -O2 $TBB_INC -DSHM_INVADE_DEBUG=8 -DTBB_USE_ASSERT -DTBB_USE_THREADING_TOOLS *.cpp -shared -o libshminvade_debug.so
g++ -std=c++14 -fPIC -O3 $TBB_INC -DNDEBUG -DSHM_INVADE_DEBUG=0 *.cpp -shared -o libshminvade.so
rm *.o


Please note that you should translate your own code with 

-DTBB_USE_THREADING_TOOLS

too, if you want to use the libs with the Intel threading tools.


If you use a debug level bigger than 1 (currently we support 0, 1 and 2), then 
SHMInvade will dump status messages to the terminal. You can configure which 
prefix it uses by redefining the variable SHM_DEBUG_PREFIX.




### SHMInvade in Peano ###

To use Peano with SHMInvade, please ensure that you compile with two flags:

-DSharedTBB -DTBBInvade

Peano's technical architecture (tarch) itself is totally agnostic of SHMInvade. 
However, all task and loop accesses going through the namespace 
peano::datatraversal do issue SHMInvade commands. 

This means that Peano itself is SHMInvade-ready but does not administer the
usage of SHMInvade. It is up to the Peano user to realise an invasion 
strategy. This can require additional work on the user side.

Besides using a different compiler argument, you have to ensure that all 
SHMInvade files are found by the compiler and linker. I usually rely on the 
environment variables TBB_INC and TBB_SHLIB to provide all information for 
TBB. If you follow this pattern, then it is convenient to add further 

-I/my/shminvade/path/without/the/shminvade/subdirectory

-L/my/shminvade/path/with/the/shminvade/subdirectory -lshminvade

statements to these environment variables.



### SHMInvade in ExaHyPE ###

In ExaHyPE, one can switch to the invasive TBB by setting the shared memory 
mode to 

SHAREDMEM=TBBInvade

Once this is done, you have to modify/extend TBB_INC and TBB_SHLIB as 
discussed above.



### Pinning ###

SHMInvade does its own pinning through TBB Observers. For most systems, there
should be no need for manual further pinning.





Meinen eigenen Kontext-Kommentar aufgreifen
Muteces
Low priority for context
Manuel deletion within templates
