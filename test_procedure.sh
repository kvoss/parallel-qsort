#/*
# * Parallel Quicksort implementation using MPI
# *
# * Author: Krzysztof Voss [shobbo@gmail.com]
# *
# */

#!/bin/tcsh

module add openmpi/gcc

foreach np ( 1 2 4 8 16 32 )
	foreach ds ( 1 4 16 64 )
		@ dss = $ds * 1000000
		set kom = `mpicc -O2 -DTEST -DNSIZE=$dss -o p-qsort p-qsort.c`

		@ ttime = 0
		foreach i ( `seq 10` )
			set t = `mpirun -n $np p-qsort | awk '{print $4*1000}'`
			@ ttime += $t
		end
		@ ttime /= 10

		echo "numproc: $np, datasize: $ds [m], time: $ttime [ms]"
	end
end

