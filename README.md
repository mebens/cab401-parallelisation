My final project for CAB401 at QUT. The task was to take a sequential, computationally-expensive program and parallelise it. I chose a provided program which computes genome similarity between pairs of bacteria in a set.

Files:

* original_raw.cpp: the originally provided program (first labelled "improved.cpp")
* original.cpp: the original program modified to work with the data files in the data folder, and a clock_gettime timing.
* modified.cpp: a modified sequential version of the program for improved performance.
* parallel.cpp: the parallel version of the program.

These were compiled using Qt Creator and MinGW 64-bit.
