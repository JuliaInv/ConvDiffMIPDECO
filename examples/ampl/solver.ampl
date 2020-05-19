#####################################################
# solver.ampl : solver name and its options
#####################################################
option solver bonmin;	# ... path/binary of the AMPLsolver binary
option bonmin_options "bonmin.algorithm=B-BB outlev=1 bonmin.time_limit=259200 max_cpu_time=259200";	# ... solver options
solve;
