# Falcon 9 Simulation
A continuous simulation of a rocket taking off and landing in a 2D environment
maps the rockets velocity and location at any given point in time, accounting
for thrust, gravity, and drag as a function of altitude-based changes in
air density. The rocket follows a mission launch profile, accelerating from a
ground state, attaining an apogee of 200km, falling back down, then executing
a landing burn before touchdown. Each stage is run multiple times to ensure
the mission profile is attained, adjusting parameters as necessary. The resultant
simulation is thus able to determine the optimal flight path for the rocket
which minimizes leftover fuel after landing, and thus maximizes the horizontal
velocity attained at apogee.

## Running the simulation
We initially did the simulation in using __Matlab__ we then ported it to __Python__ to analyze and improve the performance.

For running the python version, issue `python simulation.py --help` to get the availabe arguments

To run the matlab version (_stale_) go to `./matlab/` and  issue `make run`

## Team Members
Hamid & [Grant](https://github.com/gfenn)
