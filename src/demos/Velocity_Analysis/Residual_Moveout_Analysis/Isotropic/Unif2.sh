#! /bin/sh
# shell for uniformly sampling velocity from a layered model
set -v

# building velocity model for ray tracing
nz=51 dz=50 fz=.0  labelz="Depth (m)"
nx=81 dx=50 fx=0.0  labelx="Distance (m)"
ninf=0 npmax=201 
unif2 <input >vfile  ninf=$ninf  npmax=$npmax \
	nz=$nz dz=$dz fz=$fz nx=$nx dx=$dx fx=$fx \
	v00=1500

# building velocity derivative model for ray tracing
unif2 <input >pvfile  ninf=$ninf  npmax=$npmax \
	nz=$nz dz=$dz fz=$fz nx=$nx dx=$dx fx=$fx \
	v00=1 

exit 0

