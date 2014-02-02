#ifndef __SETTINGS_H__
#define __SETTINGS_H__

// Window settings
#define WINDOW_TITLE "Particle Fluids"
#define WINDOW_WIDTH 1280
#define WINDOW_HEIGHT 720
//#define FULLSCREEN

// OpenGL Debug Mode toggle
#define DEBUG

// Particle settings
#define NUM_PARTICLES (4096*128)

// Grid size for the acceleration grid.
// Note that if this is not a power of two, things WILL break.
#define GRID_SIZE 128

// Rendering settings
#define RESOLUTION_DIVIDER 1
#define SMOOTHING_ITERATIONS 120

#endif