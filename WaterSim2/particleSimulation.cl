/**
 * A SPH-based fluid simulation
 *
 * Known Problems:
 * - High density makes everything explode
 *
 * L. Diener, Jan. 2014
 **/

///////////////////////////////////////////// DEFINES /////////////////////////////////////////////

// A randomly chosen number
#define PI 3.1415f

// A small number
#define EPSILON 0.001f

// Downwards force
#define GRAVITY ((float4){0.0f, -9.81f, 0.0f, 0.0f})

// Simulation area. Should obviously be smaller than the grid
#define AABB_XZ 10.5f
#define AABB_Y 2.0f
#define AABB ((float4){AABB_XZ, AABB_Y, AABB_XZ, 10000000.0f})

// TODO play with those until things look reasonable
// Alternately, attempt to understand the derivation of the kernel 
// and set them to something that Makes Sense
#define KERNEL_RADIUS 0.1f
#define REST_DENSITY 500.0f
#define GAS_CONSTANT 1000.0f
#define VISCOSITY_CONSTANT 60.0f
#define DAMPENING 0.85f

// Grid information
#define GRID_SIZE_A ((int)128*2)
#define GRID_SIZE_XZ ((int)220*2)
#define GRID_SIZE_Y ((int)42*2)

#define GRID_HALF_XZ ((float)GRID_SIZE_XZ * KERNEL_RADIUS * 0.5f)
#define GRID_HALF_Y ((float)GRID_SIZE_Y * KERNEL_RADIUS * 0.5f)

#define GRID_COORD_XZ(val) ((floor((val + GRID_HALF_XZ) / KERNEL_RADIUS)))
#define GRID_COORD_Y(val) ((floor((val + GRID_HALF_Y) / KERNEL_RADIUS)))

#define GRID_ID(p) (int)(GRID_COORD_XZ((p).x) + GRID_COORD_Y((p).y) * GRID_SIZE_Y + GRID_COORD_XZ((p).z) * GRID_SIZE_Y * GRID_SIZE_XZ)

// Try not to access things outside
// #define GRID_FUMBLE(x) max(0, min((x), GRID_SIZE_A * GRID_SIZE_A * GRID_SIZE_A))

// For MAXIMUM PERFORMANCE, a noop
#define GRID_FUMBLE(x) (x)

// Maximum particles to check for density calculations and such
#define MAX_PARTICLES 8 // 32 for more accuracy, but 8 is good enough

// To prevent obvious numerical explosions, clamp various things.
#define RHO_CLAMP_LOW 500.0f
#define RHO_CLAMP_HIGH 1000000.0f
#define VELOCITY_CLAMP 8.0f

// Terrain settings
#define SCALE 8.0f
#define SCALE_HEIGHT 128.0f
#define SCALE_HEIGHT_SUB 4.0f

///////////////////////////////////////////// COLLISSIONS /////////////////////////////////////////////

// Vector to AABB
float4 AABBIntersect(float4 pos, float4 wallDistance) {
	float4 dp = wallDistance - fmax(pos, wallDistance);
	float4 dn = -(wallDistance + fmin(pos, -wallDistance));
	return sign(dp + dn) * fmax(-dp, dn);
}

float4 reflect(float4 v, float4 n) {
	return -2.0f * dot(v, n) * n + v;
}

float heightmap(float u, float v, const __global float* heightData) {
	int heightU = ((u + AABB_XZ) / (AABB_XZ * 2.0f * SCALE)) * 4096.0f;
	int heightV = ((v + AABB_XZ) / (AABB_XZ * 2.0f * SCALE)) * 4096.0f;
	int heightIdx = heightU + heightV * 4096;
	heightIdx = max(0, min(heightIdx, 4096 * 4096));
	return heightData[heightIdx] * SCALE_HEIGHT - 2.0f - SCALE_HEIGHT_SUB;
	//return (sin(u) + 1.0f) - 6.0f;
}

float4 grad(float u, float v, const __global float* heightData) {
	float off = (AABB_XZ * 2.0f * SCALE) / 4096.0f;
	float dup = heightmap(u + off, v, heightData);
	float dun = heightmap(u - off, v, heightData);
	float dvp = heightmap(u, v + off, heightData);
	float dvn = heightmap(u, v - off, heightData);
	float4 du = (float4)(u + off, dup, v, 0.0f) - (float4)(u - off, dun, v, 0.0f);
	float4 dv = (float4)(u, dvp, v + off, 0.0f) - (float4)(u, dvn, v - off, 0.0f);
	return cross(du, dv);
}

void collide(float4 p0, float4* pv, float4* vv, const __global float* heightData) {
	float4 p = *pv;
	float4 v = *vv;

	float height = heightmap(p.x, p.z, heightData);
	if(p.y < height) {
		// Find colission with heightmap by marching
		float4 pA = p0;
		float4 pB = p;
		float4 pC = pA + 0.5f * (pB - pA);
		while(length(pB - pA) > EPSILON) {
			float pCh = heightmap(pC.x, pC.z, heightData);
			if(pC.y < pCh) {
				pB = pC;
			}
			else {
				pA = pC;
			}
			pC = pA + 0.5f * (pB - pA);
		}
		float4 pI = pA;

		// Determine normal and rest length
		float rest = length(p - pI);
		float4 norm = normalize(grad(pI.x, pI.z, heightData));
		v = reflect(v, norm) * DAMPENING;
		p = p + rest * normalize(v);
		height = heightmap(p.x, p.z, heightData);
		if(p.y < height) {
			p.y = height + EPSILON;
		}
	}

	*pv = p;
	*vv = v;
}

///////////////////////////////////////////// SPH FORCES /////////////////////////////////////////////

// Density contribution of one particle at a position
float rhoSingle(float4 posI, float4 posJ) {
	float4 posDiff = posI - posJ;
	float kernelNorm = 315.0f / (64.0f * PI * pow(KERNEL_RADIUS, 9.0f));
	float kernelVal = pow(pow(KERNEL_RADIUS, 2.0f) - dot(posDiff, posDiff), 3.0f);
	return kernelNorm * kernelVal;
}

// Density at a position
float rho(float4 pos, const float4* particles, int numParticles) {
	float rhoV = 0.0f;
	for(int j = 0; j < numParticles; j++) {
		rhoV += rhoSingle(pos, particles[j]);
	}
	return fmax(RHO_CLAMP_LOW, fmin(rhoV, RHO_CLAMP_HIGH));
}

// Pressure for a given density
float pressure(float rhoV) {
	return GAS_CONSTANT * (rhoV - REST_DENSITY);
}

// Normalized pressure for a given density
float normPressure(float rhoV) {
	return pressure(rhoV) / pow(rhoV, 2.0f);
}

// Pressure gradient contribution of one particle at a position with given normalized pressures
float4 pressureDSingle(float4 posI, float4 posJ, float normPressureI, float normPressureJ) {
	float4 posDiff = posI - posJ;
	float pressureSum = normPressureI + normPressureJ;
	float kernelNorm = -45.0f / (PI * pow(KERNEL_RADIUS, 6.0f));
	float kernelVal = pow(KERNEL_RADIUS - dot(posDiff, posDiff), 2.0f);
	float4 pressureDirection = normalize(posDiff);
	return pressureSum * kernelNorm * kernelVal * pressureDirection;
}

// Pressure gradient at a position with given normalized pressure
float4 pressureD(
	float4 pos, 
	float normPressureI, 
	const float4* particles, 
	const float4* particleData, 
	int numParticles
) {
	float4 pressureDV = (float4){0.0f, 0.0f, 0.0f, 0.0f};
	for(int j = 0; j < numParticles; j++) {
		float rhoJ = particleData[j].x;
		float normPressureJ = particleData[j].y;
		pressureDV += pressureDSingle(pos, particles[j], normPressureI, normPressureJ);
	}
	return pressureDV;
}

// Viscous term contribution of one particle at a position with given densitiy and velocities
float4 viscousTermSingle(float4 posI, float4 posJ, float rhoJ, float4 vI, float4 vJ) {
	float4 posDiff = posI - posJ;
	float4 normVelDiff = (vJ - vI) / rhoJ;
	float kernelNorm = 45.0f / (PI * pow(KERNEL_RADIUS, 6.0f));
	float kernelVal = KERNEL_RADIUS - dot(posDiff, posDiff);
	return normVelDiff * kernelNorm * kernelVal;
}

// Viscous term at a position with given density and velocity
float4 viscousTerm(
	float4 pos, 
	float rhoI, 
	float4 vI, 
	const float4* particles, 
	const float4* particleVelocity, 
	const float4* particleData, 
	int numParticles
) {
	float4 viscousV = (float4){0.0f, 0.0f, 0.0f, 0.0f};
	for(int j = 0; j < numParticles; j++) {
		float rhoJ = particleData[j].x;
		viscousV += viscousTermSingle(pos, particles[j], rhoJ, vI, particleVelocity[j]);
	}
	return (VISCOSITY_CONSTANT / rhoI) * viscousV;
}

// Force for a given density, pressure gradient and viscous term
float4 sphForce(float rhoV, float4 pressureDV, float4 viscousV) {
	return viscousV - pressureDV / rhoV;
}

// Force acting on some particle
float4 getForce(
	float4 particle, 
	float4 velocity, 
	float4 data, 
	const float4* particles, 
	const float4* particleVelocity, 
	const float4* particleData,
	float timeTotal,
	int numParticles,
	float4 windDirection,
	float windPower
) {
	float4 force = GRAVITY;

	// "Wind"
	if(particle.y > -AABB_Y + 0.5f) {
		float windRunner = dot(windDirection, particle);
		float windAffect = min(1.0f, max(0.0f, particle.y + 0.5f + AABB_Y));
		float shove = (sin(timeTotal + windRunner * 2.0f) + 1.0f) * windAffect;
		force.y += sin(timeTotal + windRunner) * 4.0f * windPower;
		force -= windDirection * sin(timeTotal + windRunner) * 2.0f * windPower;
	}

	float rhoV = data.x;
	float normPressureV = data.y;
	float4 pressureDV = pressureD(particle, normPressureV, particles, particleData, numParticles);
	float4 viscousV = viscousTerm(particle, rhoV, velocity, particles, particleVelocity, particleData, numParticles);
	force += sphForce(rhoV, pressureDV, viscousV);
	return force;
}

///////////////////////////////////////////// GRID ACCESS /////////////////////////////////////////////

// Find particles in radius around a position
int getRelevantParticles(
	float4 pos,
	const __global float4* particlesIn, 
	const __global float4* particleVelocityIn,
	const __global float4* particleDataIn,
	float4* particles, 
	float4* particleVelocity,
	float4* particleData,
	const __global int* gridSize,
	const __global int* gridSum,
	const __global int* cellSelect
) {
	// Determine relevant cells
	int centralId = GRID_ID(pos);
	int cells[27];
	for(int j = 0; j < 27; j++) {
		int i = cellSelect[j];
		int xd = (i % 3) - 1;
		int yd = (((i / 3) % 3) - 1) * GRID_SIZE_Y;
		int zd = ((i / 9) - 1) * GRID_SIZE_Y * GRID_SIZE_XZ;
		cells[j] = GRID_FUMBLE(centralId + xd + yd + zd); 
	}

	// Grab particles within radius
	int totalParticles = 0;
	for(int i = 0; i < 27; i++) {
		int particleCount = gridSize[cells[i]];
		int particleStart = gridSum[cells[i]] - particleCount;
		for(int j = 0; j < particleCount; j++) {
			float4 particle = particlesIn[particleStart + j];
			float4 particleDiff = pos - particle;
			if((totalParticles < MAX_PARTICLES) && (dot(particleDiff, particleDiff) < KERNEL_RADIUS * KERNEL_RADIUS)) {
				particles[totalParticles] = particle;
				particleVelocity[totalParticles] = particleVelocityIn[particleStart + j];
				particleData[totalParticles] = particleDataIn[particleStart + j];
				totalParticles++;
			}
		}
	}

	// Done
	return totalParticles;
}

// Same as above, but positions only
int getRelevantParticlePositions(
	float4 pos,
	const __global float4* particlesIn, 
	float4* particles,
	const __global int* gridSize,
	const __global int* gridSum,
	const __global int* cellSelect
) {
	// Determine relevant cells
	int centralId = GRID_ID(pos);
	int cells[27];
	for(int j = 0; j < 27; j++) {
		int i = cellSelect[j];
		int xd = (i % 3) - 1;
		int yd = (((i / 3) % 3) - 1) * GRID_SIZE_Y;
		int zd = ((i / 9) - 1) * GRID_SIZE_Y * GRID_SIZE_XZ;
		cells[j] = GRID_FUMBLE(centralId + xd + yd + zd); 
	}

	// Grab particles within radius
	int totalParticles = 0;
	for(int i = 0; i < 27; i++) {
		int particleCount = gridSize[cells[i]];
		int particleStart = gridSum[cells[i]] - particleCount;
		for(int j = 0; j < particleCount; j++) {
			float4 particle = particlesIn[particleStart + j];
			float4 particleDiff = pos - particle;
			if((totalParticles < MAX_PARTICLES) && (dot(particleDiff, particleDiff) < KERNEL_RADIUS * KERNEL_RADIUS)) {
				particles[totalParticles] = particle;
				totalParticles++;
			}
		}
	}

	// Done
	return totalParticles;
}

///////////////////////////////////////////// INTEGRATION PRECALC /////////////////////////////////////////////

// Precalculates data for the simulation step
// Caching density / normalized pressure is a massive time / memory access saver later on
__kernel void CalculateParticleData(
	const __global float4* particles,
	__global float4* particleData,
	const __global int* gridSize,
	const __global int* gridSum,
	const __global int* cellSelect,
	const int numParticles
) {
	int gid = get_global_id(0);
	if(gid < numParticles) {
		float4 particle = particles[gid];
		float4 particlesClose[MAX_PARTICLES];
		int particlesCloseCount = getRelevantParticlePositions(
			particle,
			particles,
			particlesClose,
			gridSize,
			gridSum,
			cellSelect
		);
		float rhoV = rho(particle, particlesClose, particlesCloseCount);
		float normPressureV = normPressure(rhoV);

		// Data for simulation
		particleData[gid].x = rhoV;
		particleData[gid].y = normPressureV;
	}
}

///////////////////////////////////////////// INTEGRATION /////////////////////////////////////////////

// Calculates new positions for the particles. The actual "simulation"*
// (*) Accurate if the time step and the radius are small enough!**
// (**) Accurate if you do not think too hard about what "small enough" means!
__kernel void IntegratePosition(
	const __global float4* particlesIn, 
	__global float4* particlesOut,
	const __global float4* particleVelocityIn, 
	__global float4* particleVelocityOut,
	const __global float4* particleData,
	const __global int* gridSize,
	const __global int* gridSum,
	const __global int* cellSelect,
	const float dT,
	const float timeTotal,
	const int numParticles,
	const __global float* heightData,
	const float4 windDir,
	const float windPower
) {
	int gid = get_global_id(0);
	if(gid < numParticles) {
		// Read initial position and velocity
		float4 x0 = particlesIn[gid];
		float4 v0 = particleVelocityIn[gid];
		float4 d0 = particleData[gid];

		// Find relevant neighbours
		float4 particlesClose[MAX_PARTICLES];
		float4 particleVelocityClose[MAX_PARTICLES];
		float4 particleDataClose[MAX_PARTICLES];
		int particlesCloseCount = getRelevantParticles(
			x0,
			particlesIn,
			particleVelocityIn,
			particleData,
			particlesClose, 
			particleVelocityClose,
			particleDataClose,
			gridSize,
			gridSum,
			cellSelect
		);

		// Calculate initial force
		float4 F0 = getForce(
			x0, 
			v0, 
			d0, 
			particlesClose, 
			particleVelocityClose, 
			particleDataClose, 
			timeTotal, 
			particlesCloseCount,
			windDir,
			windPower
		);

		// "Calculate" acceleration
		float4 a0 = F0;
	
		// Integrate position
		float4 x1 = x0 + v0 * dT + 0.5f * a0 * dT * dT;
	
		// "Calculate" acceleration after position integration
		// Currently set to acceleration before position integration.
		// This is not technically correct, but not TERRIBLY bad.*
		// (*) I hope
		float4 a1 = F0;

		// Integrate velocity
		float4 v1 = v0 + 0.5f * (a0 + a1) * dT;

		// Collide with terrain and a box
		float4 intersectTest = AABBIntersect(x1, AABB);
		collide(x0, &x1, &v1, heightData);
		if(length(intersectTest) != 0.0f) {
			x1 += intersectTest + normalize(intersectTest) * 0.01f;
			v1 *= -normalize(fabs(intersectTest)) * DAMPENING;
		}
		
		 // Catastrophe prevention
		v1 = fabs(v1) > VELOCITY_CLAMP ? VELOCITY_CLAMP * sign(v1) : v1;

		// Write back
		x1.w = length(v1);
		particlesOut[gid] = x1;
		particleVelocityOut[gid] = v1;
	}
}

///////////////////////////////////////////// GRID PRECALC /////////////////////////////////////////////

// Precalculates data for the grid assignment
// Volatile and atomic operations may(tm) not be the fastest, and
// the grid data structure eats GPU memory for breakfast, but 
// conceptually this is very neat
__kernel void CalculateGridData(
	const __global float4* particles,
	__global int* particleOffset,
	volatile __global int* gridCounts,
	const int numParticles
) {
	int gid = get_global_id(0);
	if(gid < numParticles) {
		// Grid information: Cell ID and cell-relative address
		float4 particle = particles[gid];
		int gridId = GRID_FUMBLE(GRID_ID(particle));
		particleOffset[gid] = atomic_add(&gridCounts[gridId], 1);
	}
}


// Clears (zeros) grid data structure so it can be used again
__kernel void ClearGrid(__global int* grid, const int n) {
	int gid = get_global_id(0);
	if(gid < n) {
		grid[gid] = 0;
	}
}

// Reorders particles according to computed indices in the grid
__kernel void ReorderByGrid(
	const __global float4* particlesIn, 
	__global float4* particlesOut,
	const __global float4* particleVelocityIn, 
	__global float4* particleVelocityOut,
	const __global int* particleOffset,
	const __global int* gridSize, 
	const __global int* gridSum,
	const int numParticles
) {
	int gid = get_global_id(0);
	if(gid < numParticles) {
		// Figure out where to copy to
		float4 particle = particlesIn[gid];
		int gridCell = GRID_FUMBLE(GRID_ID(particle));
		int gridOffset = particleOffset[gid];
		int gridCellStart = gridSum[gridCell] - gridSize[gridCell];
		int newIndex = gridCellStart + gridOffset;

		// Copy
		particle.w = 0.0f;
		particlesOut[newIndex] = particle;
		particleVelocityOut[newIndex] = particleVelocityIn[gid];
	}
}

// A prefix sum, inclusive.
// I am getting so much value out of this prefix sum it's incredible
// Prefix sum for President 2016
__kernel void PrefixSum(const __global int* inArray, __global int* outArray, const int n, const int offset) {
	int gid = get_global_id(0);
	if(gid < n) {
		if(gid < offset) {
			outArray[gid] = inArray[gid];
		}
		else {
			outArray[gid] = inArray[gid] + inArray[gid - offset];
		}
	}
}
