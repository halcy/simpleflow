/**
 * A SPH-based fluid simulation that is sadly not particularly fast nor good.
 * L. Diener, Jan. 2014
 **/

///////////////////////////////////////////// DEFINES /////////////////////////////////////////////

// A randomly chosen number
#define PI 3.1415f

// Downwards force
#define GRAVITY ((float4){0.0f, -9.81f, 0.0f, 0.0f})

// Simulation area. Should obviously be smaller than the grid
#define AABB ((float4){3.0f, 3.0f, 3.0f, 10000000.0f})

// TODO play with those until things look reasonable
// Alternately, attempt to understand the derivation of the kernel 
// and set them to something that Makes Sense
#define KERNEL_RADIUS 0.1f
#define REST_DENSITY 500.0f
#define GAS_CONSTANT 1000.0f
#define VISCOSITY_CONSTANT 50.0f
#define DAMPENING 0.7f

// Grid information
#define GRID_SIZE 128
#define GRID_HALF (floor((float)GRID_SIZE*KERNEL_RADIUS*0.5f))
#define GRID_COORD(val) ((floor((val + GRID_HALF) / KERNEL_RADIUS)))
#define GRID_ID(p) (GRID_COORD((p).x) + GRID_COORD((p).y) * GRID_SIZE + GRID_COORD((p).z) * GRID_SIZE * GRID_SIZE)

// Maximum particles to check for density calculations and such
#define MAX_PARTICLES 32

// To prevent obvious numerical explosions, clamp velocity.
#define VELOCITY_CLAMP 4.0f

///////////////////////////////////////////// COLLISSIONS /////////////////////////////////////////////

// Vector to AABB
float4 AABBIntersect(float4 pos, float4 wallDistance) {
	float4 dp = wallDistance - max(pos, wallDistance);
	float4 dn = -(wallDistance + min(pos, -wallDistance));
	return sign(dp + dn) * max(-dp, dn);
}

///////////////////////////////////////////// SPH FORCES /////////////////////////////////////////////

// Density contribution of one particle at a position
float rhoSingle(float4 posI, float4 posJ) {
	float4 posDiff = posI - posJ;
	if(dot(posDiff, posDiff) > KERNEL_RADIUS * KERNEL_RADIUS) return 0.0f;
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
	return rhoV;
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
	if(dot(posDiff, posDiff) > KERNEL_RADIUS * KERNEL_RADIUS) return (float4){0.0f, 0.0f, 0.0f, 0.0f};
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
	if(dot(posDiff, posDiff) > KERNEL_RADIUS * KERNEL_RADIUS) return (float4){0.0f, 0.0f, 0.0f, 0.0f};
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
	int numParticles
) {
	float4 force = GRAVITY;

	// "Wind"
	if(particle.x < 2.0f) {
		float wind_affect = min(1.0f, max(0.0f, particle.y + 2.9f));
		float shove = (sin(timeTotal + particle.x * 2.0f) + 1.0f) * wind_affect;
		force.y += sin(-timeTotal + particle.x) * 7.0f;
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
	const __global int* gridSum
) {
	// Determine relevant cells
	const float4 o = (float4){KERNEL_RADIUS, -KERNEL_RADIUS, 0.0f, 0.0f};
	const int cells[27] = {
		GRID_ID(pos + o.yxxw),
		GRID_ID(pos + o.yxyw),
		GRID_ID(pos + o.yxzw),
		GRID_ID(pos + o.yyxw),
		GRID_ID(pos + o.yyyw),
		GRID_ID(pos + o.yyzw),
		GRID_ID(pos + o.yzxw),
		GRID_ID(pos + o.yzyw),
		GRID_ID(pos + o.yzzw),

		GRID_ID(pos + o.zxxw),
		GRID_ID(pos + o.zxyw),
		GRID_ID(pos + o.zxzw),
		GRID_ID(pos + o.zyxw),
		GRID_ID(pos + o.zyyw),
		GRID_ID(pos + o.zyzw),
		GRID_ID(pos + o.zzxw),
		GRID_ID(pos + o.zzyw),
		GRID_ID(pos + o.zzzw),

		GRID_ID(pos + o.xxxw),
		GRID_ID(pos + o.xxyw),
		GRID_ID(pos + o.xxzw),
		GRID_ID(pos + o.xyxw),
		GRID_ID(pos + o.xyyw),
		GRID_ID(pos + o.xyzw),
		GRID_ID(pos + o.xzxw),
		GRID_ID(pos + o.xzyw),
		GRID_ID(pos + o.xzzw),
	};

	// Grab particles within radius
	int totalParticles = 0;
	for(int i = 0; i < 27; i++) {
		int particleCount = gridSize[cells[i]];
		int particleStart = gridSum[cells[i]] - particleCount;
		for(int j = 0; j < particleCount; j++) {
			float4 particle = particlesIn[particleStart + j];
			float4 particleDiff = pos - particle;
			if(totalParticles < MAX_PARTICLES && (dot(particleDiff, particleDiff) < KERNEL_RADIUS * KERNEL_RADIUS)) {
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
	const __global int* gridSum
) {
	// Determine relevant cells
	const float4 o = (float4){KERNEL_RADIUS, -KERNEL_RADIUS, 0.0f, 0.0f};
	const int cells[27] = {
		GRID_ID(pos + o.yxxw),
		GRID_ID(pos + o.yxyw),
		GRID_ID(pos + o.yxzw),
		GRID_ID(pos + o.yyxw),
		GRID_ID(pos + o.yyyw),
		GRID_ID(pos + o.yyzw),
		GRID_ID(pos + o.yzxw),
		GRID_ID(pos + o.yzyw),
		GRID_ID(pos + o.yzzw),

		GRID_ID(pos + o.zxxw),
		GRID_ID(pos + o.zxyw),
		GRID_ID(pos + o.zxzw),
		GRID_ID(pos + o.zyxw),
		GRID_ID(pos + o.zyyw),
		GRID_ID(pos + o.zyzw),
		GRID_ID(pos + o.zzxw),
		GRID_ID(pos + o.zzyw),
		GRID_ID(pos + o.zzzw),

		GRID_ID(pos + o.xxxw),
		GRID_ID(pos + o.xxyw),
		GRID_ID(pos + o.xxzw),
		GRID_ID(pos + o.xyxw),
		GRID_ID(pos + o.xyyw),
		GRID_ID(pos + o.xyzw),
		GRID_ID(pos + o.xzxw),
		GRID_ID(pos + o.xzyw),
		GRID_ID(pos + o.xzzw),
	};

	// Grab particles within radius
	int totalParticles = 0;
	for(int i = 0; i < 27; i++) {
		int particleCount = gridSize[cells[i]];
		int particleStart = gridSum[cells[i]] - particleCount;
		for(int j = 0; j < particleCount; j++) {
			float4 particle = particlesIn[particleStart + j];
			float4 particleDiff = pos - particle;
			if(totalParticles < MAX_PARTICLES && (dot(particleDiff, particleDiff) < KERNEL_RADIUS * KERNEL_RADIUS)) {
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
			gridSum
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
	const float dT,
	const float timeTotal,
	const int numParticles
) {
	int gid = get_global_id(0);
	if(gid < 65536) {
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
			gridSum
		);

		// Calculate initial force
		float4 F0 = getForce(x0, v0, d0, particlesClose, particleVelocityClose, particleDataClose, timeTotal, particlesCloseCount);

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

		// Collide with a box
		float4 intersectTest = AABBIntersect(x1, AABB);
		if(length(intersectTest) != 0.0f) {
			x1 += intersectTest + normalize(intersectTest) * 0.01f;
			v1 *= -normalize(fabs(intersectTest)) * DAMPENING;
		}
		
		 // Catastrophe prevention
		v1 = fabs(v1) > VELOCITY_CLAMP ? VELOCITY_CLAMP * sign(v1) : v1;

		// Write back
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
	__global float4* particleData,
	__global int* gridCounts,
	const int numParticles
) {
	int gid = get_global_id(0);
	if(gid < numParticles) {
		// Grid information: Cell ID and cell-relative address
		float4 particle = particles[gid];
		int gridId = GRID_ID(particle);
		particleData[gid].z = gridId;
		particleData[gid].w = atomic_inc(&gridCounts[gridId]);
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
	const __global float4* particleDataIn,
	__global float4* particleDataOut,
	const __global int* gridSize, 
	const __global int* gridSum,
	const int numParticles
) {
	int gid = get_global_id(0);
	if(gid < 65536) {
		// Figure out where to copy to
		float4 particleInfo = particleDataIn[gid];
		int gridCell = particleInfo.z;
		int gridOffset = particleInfo.w;
		int gridCellStart = gridSum[gridCell] - gridSize[gridCell];
		int newIndex = gridCellStart + gridOffset;

		// Copy
		particlesOut[newIndex] = particlesIn[gid];
		particleVelocityOut[newIndex] = particleVelocityIn[gid];
		particleDataOut[newIndex] = particleDataIn[gid];
	}
}

// A prefix sum, inclusive.
// I am getting so much value out of this prefix sum it's incredible
// Prefix sum for President 2016
__kernel void PrefixSum(const __global int* inArray, __global int* outArray, int n, int offset) {
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
