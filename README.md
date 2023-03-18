# LJ Simulation

A small attempt at high-performance computing. This program simulates a behavior of a particle system intercating according to the Lennard-Jones potential.

https://user-images.githubusercontent.com/7950897/226097289-261508ea-6366-4c04-a0f2-240ee9ee508b.mp4

## Building

Program is built with CMake and a target build system of your choice. For parallel computing, it relies on OpenMP, so make sure it is avaialble on your system.

## Tools

`generate_input.py` - generates an input file for the program.

`parse_output.py` - procesesses an output file and generates a video (ffmpeg must be availeble on your PATH).
