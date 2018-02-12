#!/usr/bin/env python

# Noise on a circle: A demo video
# Svenk, 2018-01-25

from pylab import *
ion()

nframes = 20
for i in range(nframes):
	clf()
	N = 100
	th = linspace( 0, 2*pi, N )
	noise = rand(N) * 4.
	detr = 10
	r = detr+noise
	plot( detr*cos(th), detr*sin(th), "-", label="Ideal circle")
	plot( r*cos(th), r*sin(th), "o-", label="Noisy circle")

	gca().set_aspect(1)
	ylim(-20,20)
	xlim(-20,20)
	savefig("plot%d.png" % i)

# join pngs to movie using this command:
# ffmpeg -f image2 -r 6 -i 'plot%d.png' output.webm
