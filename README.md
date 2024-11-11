This repository contains a simulation of a crystal of Argon including its kinematic and thermodynamic properties.

A detailed description on how it was built is laid out in 'instructions.pdf' (in Polish).

'Argon_simulation.cpp' is the main file where all calculations are done, it can be launched using any up-to-date C++ compiler and needs to be located in the same folder as 'parameters.txt'; the latter file consists of examplary parameters for a basic simulation to be performed. Running 'Argon_simulation.cpp' will produce 2 text files -- with thermodynamics functions of state ('state.txt') and atomic positions ('xyz.txt'). Those are necessary for further analysis.

In order to visualise the simulation, 'animate.zip' needs to be unpacked. Inside, there is 'animateAr.cpp' file using OpenGL/GLUT library fully prepared for it. First, the program must be build using included 'Makefile'. Then, 2 text files -- 'params' and 'xyz' -- must be filled with what was produced as a result of 'Argon_simulation.cpp' calculations. When run, 'animateAr.cpp' will produce a visual input to inserted numbers.

Finally, 'Argon_simulation.ipynb' is a Python code containing analysis of the programme based on the points listed in 'instructions.pdf'.
