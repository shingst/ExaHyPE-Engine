##
# @file This file is part of the ExaHyPE project.
# @author ExaHyPE Group (exahype@lists.lrz.de)
#
# @section LICENSE
#
# Copyright (c) 2016  http://exahype.eu
# All rights reserved.
#
# The project has received funding from the European Union's Horizon 
# 2020 research and innovation programme under grant agreement
# No 671698. For copyrights and licensing, please consult the webpage.
#
# Released under the BSD 3 Open Source License.
# For the full license text, see LICENSE.txt
#
#
# @section DESCRIPTION
#
# Lists supported configurations
#
# Comments on the meaning of the arguments are copied from 
# https://github.com/hfp/libxsmm
# 

architectures = [ 
                    "noarch", 
                    "wsm",               # westmere
                    "snb",               # sandy bridge
                    "hsw",               # haswell
                    "knc",               # knights corner
                    "knl"                # knights landing
                ]

prefetchStrategies = [
                        "nopf",              # no prefetching at all, just 3 inputs (*A, *B, *C)
                        "pfsigonly",         # just prefetching signature, 6 inputs (*A, *B, *C, *A', *B', *C')
                        "BL2viaC",           # uses accesses to *C to prefetch *B'
                        "AL2",               # uses accesses to *A to prefetch *A'
                        "curAL2",            # prefetches current *A ahead in the kernel
                        "AL2_BL2viaC",       # combines AL2 and BL2viaC
                        "curAL2_BL2viaC",    # combines curAL2 and BL2viaC
                        "AL2jpst",           # aggressive *A' prefetch of first rows without any structure
                        "AL2jpst_BL2viaC"    # combines AL2jpst and BL2viaC
                     ]

numerics   =  [
                "linear",
                "nonlinear"
              ]
