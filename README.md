simpleflow
==========

A simple OpenCL / OpenGL fluid simulation and renderer.

Warning: Requires a fairly powerful GPU!

Controls: 
* Movement: 
	* Mouse to look
	* WASD to move around
	* Q/E to go up/down
* Waves: 
	* R/F to make stronger/weaker
	* V to set power to 0
	* Y/C to rotate direction
	* X to reverse direction
* Features
	* 1/3 to decrease/increase simulation time step
	* 2 to reverse time step (does not work very well)
	* 4/5 to decrease/increase smoothing iterations
	* 6 to toggle liquid shading
* Misc:
	* P to pause everything for a second
	* Escape to quit

To make, open the .sln file in Visual Studio and compile, everything 
that is required should be included in the repository. Currently,
only windows is supported, but most of the code (everything that is
not user input handling) is platform independent and should work on
any platform with a new enough OpenGL.

![Screenshot](http://aka-san.halcy.de/share/Particle_Fluids_2014-02-04_04-07-34.png)

Particles are simulated using an OpenCL-based Smoothed Particle
Hydrodynamics implementation, colliding with a height map. The 
resulting particles are then rendered and smoothed using screen-space 
curvature smoothing, and eventually shaded and combined with the 
heightmap geometry.
