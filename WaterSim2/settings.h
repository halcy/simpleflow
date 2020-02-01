#ifndef __SETTINGS_H__
#define __SETTINGS_H__

// Window settings
#define WINDOW_TITLE "Particle Fluids"
#define WINDOW_WIDTH 1280
#define WINDOW_HEIGHT 720
// #define FULLSCREEN

// OpenGL Debug Mode toggle
// #define DEBUG

// Particle settings
#define NUM_PARTICLES (4096*64*4)

// Grid size for the acceleration grid.
// Note that if this is not a power of two, things WILL break.
#define GRID_SIZE 256

// Bounding box and terrain settings. Important for aligning the terrain
#define AABB_XZ 10.5f

// Simulation settings
#define TIMESTEP 0.01f
#define ITERS_PER_FRAME 1

// Rendering settings
#define RESOLUTION_DIVIDER 1
#define SMOOTHING_ITERATIONS 240

// Camera settings
#define CAM_MOVESPEED 0.2f
#define CAM_ROTSPEED 0.005f

#endif