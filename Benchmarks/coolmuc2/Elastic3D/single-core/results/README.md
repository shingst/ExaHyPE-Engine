Optimisations:
'hsw-O2-novec', 'hsw-O2-vec', 'hsw-O3-vec', 'noarch-O2-novec', 'noarch-O2-vec', 'noarch-O3-vec'
which maps to the following ICC flags:
noarch => // nothing
hsw    => -xCORE-AVX2
knl    => -xMIC-AVX512
novec  => -no-vec -no-fma
vec    => -vec -fma

Time is normalised per time step and (p+1)^3, p: polynomial order
