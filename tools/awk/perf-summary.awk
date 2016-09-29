# parses Enzo-p output and prints the following
#
#   number of solver iterations
#   time per iteration
#   number of blocks
#   number of particles
#   time per block
#   time per particle
#   time per zone
#   current memory
#   highest memory
#

BEGIN {
    done = 0;
    block_zones_real = 1;
    block_zones_ghost = 1;
    i_root_blocks = 0;
    i_root_size = 0;
    bytes_high_start = 0;
    num_blocks_start = 0;
}

/processors/{num_processors = $6; }
# Block size

/ghost_depth/ { ghosts = $6;}
/root_blocks/ {
    n = $6;
    block_zones_real = block_zones_real * n;
    block_zones_ghost = block_zones_ghost * (2*ghosts * n);
}

# Linear Solver
/iter 0000/  {time_start=$2;}
/Simulation cycle 0000/{cycle_start=$2;}
/Simulation cycle /{cycle_last=$2;}
/final iter/ {
    max_iter      = $7;
    time_stop     = $2;
    time_per_iter = (time_stop-time_start)/max_iter; 
}
/Simulation cycle/ {cycle = $5;}
/bytes-high / {
    bytes_high = $6;
    if (bytes_high_start == 0) bytes_high_start = $6;
}
/bytes-highest/ {  bytes_highest = $6; }
/num-blocks/ {
    num_blocks = $6;
    if (num_blocks_start == 0) num_blocks_start = $6;
}
/num-particles/ {
    num_particles = $NF;
    if (num_particles_start == 0) {
	num_particles_start = num_particles;
    }
}
/Mesh:root_blocks/ {
    root_blocks[i_root_blocks] = $6;
    i_root_blocks++;
}
/Mesh:root_size/ {
    root_size[i_root_size] = $6;
    i_root_size++;
}
/Adapt:max_level/ { max_level = $6 }

/Mesh:root_rank/ {root_rank = $6; }
/END ENZO/ {
    done = 1;
    time_final = $2;
}

/simulation time-usec/  {time_simulation  = $6}
/initial time-usec/  {time_initial  = $6}
/cycle time-usec/    {time_cycle    = $6}
/compute time-usec/  {time_compute  = $6};
/adapt time-usec/    {time_adapt    = $6};
/refresh time-usec/  {time_refresh  = $6};
/output time-usec/   {time_output   = $6};
/stopping time-usec/ {time_stopping = $6};

END {

    format = "%20s: %g\n";
    format3 = "%20s: %d %d %d\n";


    if (done == 1) {
	print ("Run completed");
    } else {
	print ("RUN INCOMPLETE");
    }

    print "\nSUMMARY\n"
    
    printf (format, "Total time",time_final);
    printf (format, "Cycles",cycle);
    printf (format, "Init time",cycle_start);
    printf (format, "Cycle time",(cycle_last-cycle_start)/cycle);

    printf "\nPHASES\n"

    printf (format, "Simulation",time_simulation/num_processors*0.000001);
    printf (format, "Initial",time_initial/num_processors*0.000001);
    printf (format, "Cycling",time_cycle/num_processors*0.000001);
    printf("\n");
    printf (format, "Compute",time_compute/num_processors*0.000001);
    printf (format, "Adapt",time_adapt/num_processors*0.000001);
    printf (format, "Refresh",time_refresh/num_processors*0.000001);
    printf (format, "Output",time_output/num_processors*0.000001);
    printf (format, "Stopping",time_stopping/num_processors*0.000001);
    
    printf "\nMESH HIERARCHY\n"
    
    printf (format, "rank",root_rank);
    printf (format3,"root blocks",root_blocks[0],root_blocks[1],root_blocks[2]);
    printf (format3,"root size",root_size[0],root_size[1],root_size[2]);
    printf (format3,"effective size",
	    root_size[0]*2^max_level,
	    root_size[1]*2^max_level,
	    root_size[2]*2^max_level);
    printf (format3,"block size",(root_size[0]/root_blocks[0]),
	    (root_size[1]/root_blocks[1]),
	    (root_size[2]/root_blocks[2]));
    printf (format, "max-level",max_level);
    printf (format, "num-blocks",num_blocks);
    printf (format, "blocks change",num_blocks / num_blocks_start);

    printf ("\nPARALLELISM\n");
    printf (format, "num processes", num_processors);

    printf ("\nMEMORY\n");
    printf (format, "bytes-high",bytes_high);
    printf (format, "bytes-highest",bytes_highest);
    printf (format, "memory change",bytes_high / bytes_high_start);


    printf ("\nLINEAR SOLVER");
    if (time_per_iter > 0) {
	print;
	printf (format, "time per block",time_per_iter/num_blocks);
	printf (format, "time per real zone",time_per_iter/(num_blocks*block_zones_real));
	printf (format, "time per zone",time_per_iter/(num_blocks*block_zones_ghost));
	printf (format, "solver max-iter",max_iter);
	printf (format, "time per iter",time_per_iter);
    } else {
	print " not called"
    }

    printf ("\nPARTICLES");
    
    if (num_particles > 0) {
	print
	printf (format,"num-particles",num_particles);
	printf (format, "time per particle",time_per_iter/num_particles);
	printf (format, "particles change",num_particles / num_particles_start);
    } else {
	print " not used"
    }
}

