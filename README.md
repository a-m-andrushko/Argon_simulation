This repository contains a simulation of a crystal of Argon including its kinematic and thermodynamic properties.

A detailed description on how it was built is laid out in _instructions.pdf_ (in Polish).

_Argon_simulation.cpp_ is the main file where all calculations are done, it can be launched using any up-to-date C++ compiler and needs to be located in the same folder as _parameters.txt_; the latter file consists of examplary parameters for a basic simulation to be performed. Running _Argon_simulation.cpp_ will produce 2 text files -- with thermodynamics functions of state (_state.txt_) and atoms' positions (_xyz.txt_). Those are necessary for further analysis.

In order to visualise the simulation, _animate.zip_ needs to be unpacked. Inside, there is _animateAr.cpp_ file using OpenGL/GLUT library fully prepared for it. First, the program must be build using included _Makefile_. Then, 2 text files -- _params_ and _xyz_ -- must be filled with what was produced as a result of _Argon_simulation.cpp_ calculations. When run, _animateAr.cpp_ will produce a visual input to inserted numbers.

Finally, _Argon_simulation.ipynb_ is a Python code containing analysis of the programme based on the points listed in _instructions.pdf_.
