/**
 * A particle fluid renderer.
 */

// Settings
#include "settings.h"

// "THIS FUNCTION IS UNSAFE"
#pragma warning(disable: 4996)

// Include files with neat things.
#include "glhelpers.h"
#include "Vector.h"
#include "pgmloader.h"

// Include CL/GL interop tools.
#include "clgl.h"

// Data like shaders and such
#include "data.h"

// Vertex array object
GLuint vertexArray;

// Time variable
float time;

// Forward declaration
void draw();
void update();
void updateI(int);
void handleKeypress(unsigned char k, int x, int y);

// Particles
typedef struct particle_t {
	float x;
	float y;
	float z;
	float w;
} particle;

// [-0.5, 0.5] RNG
float centeredUnitRand() {
	return (((float)rand()/(float)RAND_MAX) - 0.5f);
}

//////////////////////////////// INIT FUNCTIONS ////////////////////////////////

// Make a window for doing OpenGL
void makeWindow(int argc, char** argv) {
	glutInit(&argc, argv);
	glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH);
	glutInitContextVersion(4, 2);
	#ifdef DEBUG
		glutInitContextFlags(GLUT_CORE_PROFILE | GLUT_DEBUG);
	#else
		glutInitContextFlags(GLUT_COMPATIBILITY_PROFILE);
	#endif

	// Open window / full screen
	#ifdef FULLSCREEN
		char modeString[255] = "";
		sprintf(modeString, "%dx%d:24", WINDOW_WIDTH, WINDOW_HEIGHT);
		glutGameModeString(modeString);
		glutEnterGameMode();
	#else
		glutInitWindowSize(WINDOW_WIDTH, WINDOW_HEIGHT);
		glutCreateWindow(WINDOW_TITLE);
	#endif

	// Use GLEW, and make sure it imports EVERYTHING
	glewExperimental = GL_TRUE;
	glewInit();

	#ifdef DEBUG
		registerGlDebugLogger(GL_DEBUG_SEVERITY_MEDIUM);
	#endif

	// Set up OpenGL features
	glEnable(GL_DEPTH_TEST);
	glClearDepth(1.0f);
	glDepthFunc(GL_LEQUAL);
	glEnable(GL_CULL_FACE);
	glCullFace(GL_BACK);
}

// Functions that glut calls for us
void setupGlutCallbacks() {
	glutDisplayFunc(draw);
	glutTimerFunc(15, updateI, 0);
	glutKeyboardFunc(handleKeypress);
}

// Initialize buffers, textures, etc.
void initObjects() {
	// Put a bunch of particles into space
	particle* particleData = (particle*)malloc(sizeof(particle) * NUM_PARTICLES);
	GLuint* particleElements = (GLuint*)malloc(sizeof(GLuint) * NUM_PARTICLES);

	srand(666);
	for(int i = 0; i < NUM_PARTICLES; i++) {
		particleData[i].x = centeredUnitRand() * 3.0f * 1.99f;
		particleData[i].y = centeredUnitRand() * 3.0f * 1.99f;
		particleData[i].z = centeredUnitRand() * 3.0f * 1.99f;
		particleData[i].w = 0.0f;
		particleElements[i] = i;
	}
	for(int i = 0; i < 2; i++) {
		particles.vertexBuffer[i] = makeBO(
			GL_ARRAY_BUFFER,
			particleData,
			sizeof(particle) * NUM_PARTICLES,
			GL_DYNAMIC_DRAW
		);
	}
	particles.elementBuffer = makeBO(
		GL_ELEMENT_ARRAY_BUFFER,
		particleElements,
		sizeof(GLuint) * NUM_PARTICLES,
		GL_STATIC_DRAW
	);

	free(particleData);
	free(particleElements);

	// Prepare a screen quad to render postprocessed things.
	Vector quadData[] = {
		{-1.0f, -1.0f, 0.0f},
		{ 1.0f, -1.0f, 0.0f},
		{ 1.0f,  1.0f, 0.0f},
		{-1.0f,  1.0f, 0.0f}
	};
	GLuint quadElements[] = {0, 1, 3, 1, 2, 3};

	screenQuad.vertexBuffer = makeBO(
		GL_ARRAY_BUFFER,
		quadData,
		sizeof(Vector) * 4,
		GL_STATIC_DRAW
	);
	screenQuad.elementBuffer = makeBO(
		GL_ELEMENT_ARRAY_BUFFER,
		quadElements,
		sizeof(GLuint) * 6,
		GL_STATIC_DRAW
	);

	// Load textures.
	terrain.envTexture = loadTexture("skymap_b.tga");

	// Create a VAO and bind it
	glGenVertexArrays(1, &vertexArray);
	glBindVertexArray(vertexArray);

	// Prepare a lot of zero'd data
	particle* zeroData = (particle*)malloc(sizeof(particle) * NUM_PARTICLES);
	cl_int* zeroGrid = (cl_int*)malloc(sizeof(cl_int) * GRID_SIZE * GRID_SIZE * GRID_SIZE);
	for(int i = 0; i < NUM_PARTICLES; i++) {
		zeroData[i].x = 0.0f;
		zeroData[i].y = 0.0f;
		zeroData[i].z = 0.0f;
		zeroData[i].w = 0.0f;
	}
	for(int i = 0; i < GRID_SIZE * GRID_SIZE * GRID_SIZE; i++) {
		zeroGrid[i] = 0;
	}

	// Share some buffers with OpenCL and create some more
	for(int i = 0; i < 2; i++) {
		particles.particleBuffer[i] = sharedBuffer(particles.vertexBuffer[i], CL_MEM_READ_WRITE);
		particles.velocityBuffer[i] = clCreateBuffer(
			clContext(), 
			CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
			NUM_PARTICLES * sizeof(cl_float) * 4, 
			zeroData, 
			NULL
		);
		particles.gridBuffer[i] = clCreateBuffer(
			clContext(), 
			CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
			sizeof(cl_int) * GRID_SIZE * GRID_SIZE * GRID_SIZE, 
			zeroGrid, 
			NULL
		);
	}
	particles.dataBuffer = clCreateBuffer(
		clContext(), 
		CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
		NUM_PARTICLES * sizeof(cl_float) * 4, 
		zeroData, 
		NULL
	);
	particles.offsetBuffer = clCreateBuffer(
		clContext(), 
		CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
		NUM_PARTICLES * sizeof(cl_int), 
		zeroData, 
		NULL
	);
	particles.gridSizeBuffer = clCreateBuffer(
		clContext(), 
		CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
		sizeof(cl_int) * GRID_SIZE * GRID_SIZE * GRID_SIZE, 
		zeroGrid, 
		NULL
	);
	for(int i = 0; i < 27; i++) {
		zeroGrid[i] = i;
	}
	particles.cellSelectBuffer = clCreateBuffer(
		clContext(), 
		CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
		sizeof(cl_int) * 27, 
		zeroGrid, 
		NULL
	);
	free(zeroData);
	free(zeroGrid);

	// Load terrain
	terrain.heightData = loadPGM("grand_canyon.pgm", 4096, 4096);
	particles.terrainBuffer = clCreateBuffer(
		clContext(), 
		CL_MEM_READ_ONLY | CL_MEM_COPY_HOST_PTR, 
		sizeof(cl_float) * 4096 * 4096, 
		terrain.heightData, 
		NULL
	);
	
	// Make terrain texture
	terrain.heightTexture = genFloatTexture(terrain.heightData, 4096, 4096);

	// Make terrain geometry
	PaddedVector* terrainVertices = (PaddedVector*)malloc(sizeof(PaddedVector) * 512 * 512);
	for(int x = 0; x < 4096; x += 8) {
		for(int y = 0; y < 4096; y += 8) {
			int xx = x / 8;
			int yy = y / 8;
			PaddedVector v = PadVector(MakeVector(
				xx * (1.0f / 512.0f) * (AABB_XZ * 2.0f) - AABB_XZ, 
				terrain.heightData[(x / 8) + 4096 * (y / 8)] * 128.0f - 3.0f - 3.0f, 
				yy * (1.0f / 512.0f) * (AABB_XZ * 2.0f) - AABB_XZ
			));
			v.pad = 1.0f;
			terrainVertices[xx + 512 * yy] = v;
		}
	}
	terrain.vertexBuffer = makeBO(
		GL_ARRAY_BUFFER,
		terrainVertices,
		sizeof(PaddedVector) * 512 * 512,
		GL_STATIC_DRAW
	);
	free(terrainVertices);

	GLuint* terrainElements = (GLuint*)malloc(sizeof(GLuint) * 512 * 512 * 6);
	int quadIndex = 0;
	for(int x = 0; x < 511; x++) {
		for(int y = 0; y < 511; y++) {
			terrainElements[quadIndex * 6 + 0] = (x + 0) + (y + 0) * 512;
			terrainElements[quadIndex * 6 + 1] = (x + 1) + (y + 1) * 512;
			terrainElements[quadIndex * 6 + 2] = (x + 1) + (y + 0) * 512;
			terrainElements[quadIndex * 6 + 3] = (x + 0) + (y + 0) * 512;
			terrainElements[quadIndex * 6 + 4] = (x + 0) + (y + 1) * 512;
			terrainElements[quadIndex * 6 + 5] = (x + 1) + (y + 1) * 512;
			quadIndex++;
		}
	}
	terrain.elementBuffer = makeBO(
		GL_ELEMENT_ARRAY_BUFFER,
		terrainElements,
		sizeof(GLuint) * 512 * 512 * 6,
		GL_STATIC_DRAW
	);

	free(terrainElements);
}

// Create an OpenCL program and kernels
void initPrograms() {
	openCLKernels.program = clProgramFromFile("particleSimulation.cl");
	openCLKernels.gridClearKernel = clCreateKernel(openCLKernels.program, "ClearGrid", NULL);
	openCLKernels.gridKernel = clCreateKernel(openCLKernels.program, "CalculateGridData", NULL);
	openCLKernels.prefixSumKernel = clCreateKernel(openCLKernels.program, "PrefixSum", NULL);
	openCLKernels.gridReorderKernel = clCreateKernel(openCLKernels.program, "ReorderByGrid", NULL);
	openCLKernels.dataKernel = clCreateKernel(openCLKernels.program, "CalculateParticleData", NULL);
	openCLKernels.simulationKernel = clCreateKernel(openCLKernels.program, "IntegratePosition", NULL);
}

// Set up an FBO to render to
void initFramebuffers() {
	// Depth buffer
    GLuint depthBuffer;
    glGenRenderbuffers(1, &depthBuffer);
    glBindRenderbuffer(GL_RENDERBUFFER, depthBuffer);
    glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, WINDOW_WIDTH, WINDOW_HEIGHT);

	// FBO for background
	glGenFramebuffers(1, &framebuffers.backgroundFBO);
    glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.backgroundFBO);
    framebuffers.backgroundTexture = makeTextureBuffer(WINDOW_WIDTH, WINDOW_HEIGHT, GL_RGBA, GL_RGBA32F);
    glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framebuffers.backgroundTexture, 0);

	// Attach depth
    glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBuffer);

	// Depth buffer, low-res.
	GLuint depthBufferLowres;
    glGenRenderbuffers(1, &depthBufferLowres);
    glBindRenderbuffer(GL_RENDERBUFFER, depthBufferLowres);
	glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER);

	// FBOs for particle depth
	for(int i = 0; i < 2; i++) {
		glGenFramebuffers(1, &framebuffers.particleFBO[i]);
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.particleFBO[i]);

		framebuffers.particleTexture[i] = makeTextureBuffer(WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER, GL_RED, GL_R32F);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framebuffers.particleTexture[i], 0);

		// Attach depth
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBufferLowres);
	}

	// FBOs for particle velocity
	glGenFramebuffers(1, &framebuffers.velocityFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.velocityFBO);

	framebuffers.particleVelocityTexture = makeTextureBuffer(WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER, GL_RED, GL_R32F);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framebuffers.particleVelocityTexture, 0);

	// Attach depth
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBufferLowres);

	// FBOs for particle thickness
	for(int i = 0; i < 2; i++) {
		glGenFramebuffers(1, &framebuffers.particleThicknessFBO[i]);
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.particleThicknessFBO[i]);

		framebuffers.particleThicknessTexture[i] = makeTextureBuffer(WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER, GL_RED, GL_R32F);
		glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framebuffers.particleThicknessTexture[i], 0);

		// Attach depth
		glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBufferLowres);
	}

    glBindFramebuffer(GL_FRAMEBUFFER, 0);

	// FBO for final particle color
	glGenFramebuffers(1, &framebuffers.particleColorFBO);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.particleColorFBO);

	framebuffers.particleColorTexture = makeTextureBuffer(WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER, GL_RGBA, GL_RGBA32F);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, framebuffers.particleColorTexture, 0);

	// Attach depth
	glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, depthBufferLowres);
}

// Initialize shaders.
void initShaders() {
	// Load a simple Vertex/Fragment shader
	GLuint vertexShader = loadShader(GL_VERTEX_SHADER, "simple.vert");
	GLuint fragmentShader = loadShader(GL_FRAGMENT_SHADER, "simple.frag");
	objectShader.shaderProgram = makeShaderProgram(vertexShader, fragmentShader);

	// Get locations of attributes and uniforms used inside.
	objectShader.vertexPosition = glGetAttribLocation(objectShader.shaderProgram, "vertex");

	objectShader.modelviewMatrix = glGetUniformLocation(objectShader.shaderProgram, "modelview");
	objectShader.normalviewMatrix = glGetUniformLocation(objectShader.shaderProgram, "normalview");
	objectShader.projectionMatrix = glGetUniformLocation(objectShader.shaderProgram, "projection");
	objectShader.terrainTexture = glGetUniformLocation(objectShader.shaderProgram, "terrainTexture");

	// Bind output variables
	glBindFragDataLocation(objectShader.shaderProgram, 0, "outColor");

	// Create particle depth rendering shader
	vertexShader = loadShader(GL_VERTEX_SHADER, "particles.vert");
	fragmentShader = loadShader(GL_FRAGMENT_SHADER, "particledepth.frag");
	particleShader.shaderProgram = makeShaderProgram(vertexShader, fragmentShader);
	
	// Get locations of attributes and uniforms used inside.
	particleShader.vertexPosition = glGetAttribLocation(particleShader.shaderProgram, "vertex");

	particleShader.modelviewMatrix = glGetUniformLocation(particleShader.shaderProgram, "modelview");
	particleShader.projectionMatrix = glGetUniformLocation(particleShader.shaderProgram, "projection");
	particleShader.screenSize = glGetUniformLocation(particleShader.shaderProgram, "screenSize");
	particleShader.terrainTexture = glGetUniformLocation(particleShader.shaderProgram, "terrainTexture");

	// Bind output variables
	glBindFragDataLocation(particleShader.shaderProgram, 0, "particleDepth");

	// Create particle thickness rendering shader
	vertexShader = loadShader(GL_VERTEX_SHADER, "particles.vert");
	fragmentShader = loadShader(GL_FRAGMENT_SHADER, "particlethickness.frag");
	particleThicknessShader.shaderProgram = makeShaderProgram(vertexShader, fragmentShader);

	// Get locations of attributes and uniforms used inside.
	particleThicknessShader.vertexPosition = glGetAttribLocation(particleThicknessShader.shaderProgram, "vertex");

	particleThicknessShader.modelviewMatrix = glGetUniformLocation(particleThicknessShader.shaderProgram, "modelview");
	particleThicknessShader.projectionMatrix = glGetUniformLocation(particleThicknessShader.shaderProgram, "projection");
	particleThicknessShader.screenSize = glGetUniformLocation(particleThicknessShader.shaderProgram, "screenSize");

	particleThicknessShader.terrainTexture = glGetUniformLocation(particleThicknessShader.shaderProgram, "terrainTexture");

	// Bind output variables
	glBindFragDataLocation(particleThicknessShader.shaderProgram, 0, "particleThickness");

	// Create particle thickness rendering shader
	vertexShader = loadShader(GL_VERTEX_SHADER, "particles.vert");
	fragmentShader = loadShader(GL_FRAGMENT_SHADER, "particlevelocity.frag");
	particleVelocityShader.shaderProgram = makeShaderProgram(vertexShader, fragmentShader);

	// Get locations of attributes and uniforms used inside.
	particleVelocityShader.vertexPosition = glGetAttribLocation(particleVelocityShader.shaderProgram, "vertex");

	particleVelocityShader.modelviewMatrix = glGetUniformLocation(particleVelocityShader.shaderProgram, "modelview");
	particleVelocityShader.projectionMatrix = glGetUniformLocation(particleVelocityShader.shaderProgram, "projection");
	particleVelocityShader.screenSize = glGetUniformLocation(particleVelocityShader.shaderProgram, "screenSize");

	// Bind output variables
	glBindFragDataLocation(particleVelocityShader.shaderProgram, 0, "velocityMap");

	// Create curvature flow shader
	vertexShader = loadShader(GL_VERTEX_SHADER, "quad.vert");
	fragmentShader = loadShader(GL_FRAGMENT_SHADER, "curvatureflow.frag");
	curvatureFlowShader.shaderProgram = makeShaderProgram(vertexShader, fragmentShader);

	// Attributes / uniforms
	curvatureFlowShader.vertexPosition = glGetAttribLocation(curvatureFlowShader.shaderProgram, "vertex");
	curvatureFlowShader.particleTexture = glGetUniformLocation(curvatureFlowShader.shaderProgram, "particleTexture");
	curvatureFlowShader.projectionMatrix = glGetUniformLocation(curvatureFlowShader.shaderProgram, "projection");
	curvatureFlowShader.screenSize = glGetUniformLocation(curvatureFlowShader.shaderProgram, "screenSize");

	// Output
	glBindFragDataLocation(curvatureFlowShader.shaderProgram, 0, "outDepth");

	// Create liquid shading shader
	vertexShader = loadShader(GL_VERTEX_SHADER, "quad.vert");
	fragmentShader = loadShader(GL_FRAGMENT_SHADER, "liquidshade.frag");
	liquidShadeShader.shaderProgram = makeShaderProgram(vertexShader, fragmentShader);

	// Attributes / uniforms
	liquidShadeShader.vertexPosition = glGetAttribLocation(liquidShadeShader.shaderProgram, "vertex");
	liquidShadeShader.particleTexture = glGetUniformLocation(liquidShadeShader.shaderProgram, "particleTexture");
	liquidShadeShader.particleThicknessTexture = glGetUniformLocation(liquidShadeShader.shaderProgram, "particleThicknessTexture");
	liquidShadeShader.environmentTexture = glGetUniformLocation(liquidShadeShader.shaderProgram, "environmentTexture");
	liquidShadeShader.velocityTexture = glGetUniformLocation(liquidShadeShader.shaderProgram, "velocityTexture");
	liquidShadeShader.projectionMatrix = glGetUniformLocation(liquidShadeShader.shaderProgram, "projection");
	liquidShadeShader.modelviewMatrix = glGetUniformLocation(liquidShadeShader.shaderProgram, "modelview");
	liquidShadeShader.screenSize = glGetUniformLocation(liquidShadeShader.shaderProgram, "screenSize");

	// Output
	glBindFragDataLocation(liquidShadeShader.shaderProgram, 0, "outColor");

	// Create composition shader
	vertexShader = loadShader(GL_VERTEX_SHADER, "quad.vert");
	fragmentShader = loadShader(GL_FRAGMENT_SHADER, "compose.frag");
	compositionShader.shaderProgram = makeShaderProgram(vertexShader, fragmentShader);

	// Attributes / uniforms
	compositionShader.vertexPosition = glGetAttribLocation(compositionShader.shaderProgram, "vertex");
	compositionShader.backgroundTexture = glGetUniformLocation(compositionShader.shaderProgram, "backgroundTexture");
	compositionShader.terrainTexture = glGetUniformLocation(compositionShader.shaderProgram, "terrainTexture");
	compositionShader.particleTexture = glGetUniformLocation(compositionShader.shaderProgram, "particleTexture");
	compositionShader.modelviewMatrix = glGetUniformLocation(compositionShader.shaderProgram, "modelview");

	// Output
	glBindFragDataLocation(compositionShader.shaderProgram, 0, "outColor");
}

//////////////////////// STATUS UPDATE FUNCTIONS ////////////////////////////////////////

// Update status variables.
// Called every 15ms, unless the PC is too slow.
void updateI(int) { update(); }

void update() {
	// Update things here.
	time += 0.02f;

	// Redraw screen now, please, and call again in 15ms.
	glutPostRedisplay();
	glutTimerFunc(15, updateI, 0);

	
	// Input is win32 only
	HWND window = GetActiveWindow();
	POINT p;
	GetCursorPos(&p);
	ScreenToClient(window, &p);
	glutWarpPointer(WINDOW_WIDTH / 2, WINDOW_HEIGHT / 2);

	float angleX = (p.x - (WINDOW_WIDTH / 2)) * CAM_ROTSPEED;
	Quaternion rotX = RotationQuaternion(-angleX, MakeVector(0, 1, 0));
	camera.front = TransformVector(RotationMatrixFromQuaternion(rotX), camera.front);

	camera.elevation += (p.y - (WINDOW_HEIGHT / 2)) * CAM_ROTSPEED;
	camera.elevation = max(-0.9f, min(camera.elevation, 0.9f));

	if(GetAsyncKeyState('W') != 0) {
		camera.pos = VectorAdd(camera.pos, VectorMul(camera.front, CAM_MOVESPEED));
	}

	if(GetAsyncKeyState('S') != 0) {
		camera.pos = VectorAdd(camera.pos, VectorMul(camera.front, -CAM_MOVESPEED));
	}

	if(GetAsyncKeyState('E') != 0) {
		camera.pos = VectorAdd(camera.pos, VectorMul(camera.up, CAM_MOVESPEED));
	}

	if(GetAsyncKeyState('Q') != 0) {
		camera.pos = VectorAdd(camera.pos, VectorMul(camera.up, -CAM_MOVESPEED));
	}

	if(GetAsyncKeyState('D') != 0) {
		camera.pos = VectorAdd(camera.pos, VectorMul(VectorCross(camera.front, camera.up), CAM_MOVESPEED));
	}

	if(GetAsyncKeyState('A') != 0) {
		camera.pos = VectorAdd(camera.pos, VectorMul(VectorCross(camera.front, camera.up), -CAM_MOVESPEED));
	}
}


// Shuffling helper
void shuffle(int *array, size_t n) {
	if (n > 1) {
		size_t i;
		for (i = 0; i < n - 1; i++) {
			size_t j = i + rand() / (RAND_MAX / (n - i) + 1);
			int t = array[j];
			array[j] = array[i];
			array[i] = t;
		}
	}
}

// Draw the scene to the screen
void draw() {
	/////////////////////// PART 1: SIMULATION /////////////////////////////////

	// Grab buffers for OpenCL
	acquireGLBuffer(particles.particleBuffer[particles.currentBuffer]);
	acquireGLBuffer(particles.particleBuffer[1 - particles.currentBuffer]);

	// Prepare to run some kernels
	cl_int numParticles = NUM_PARTICLES;
	cl_int gridElements = GRID_SIZE * GRID_SIZE * GRID_SIZE;

	cl_uint workSize[3] = {numParticles, 0, 0};
	cl_uint gridWorkSize[3] = {gridElements, 0, 0};
	cl_uint workgroupSize[3] = {1024, 0, 0};

	// Clear grid
	clSetKernelArg(openCLKernels.gridClearKernel, 0, sizeof(cl_mem), &particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.gridClearKernel, 1, sizeof(cl_int), &gridElements);
	clRunKernel(openCLKernels.gridClearKernel, gridWorkSize, workgroupSize);

	// Compute grid positions
	clSetKernelArg(openCLKernels.gridKernel, 0, sizeof(cl_mem), &particles.particleBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridKernel, 1, sizeof(cl_mem), &particles.offsetBuffer);
	clSetKernelArg(openCLKernels.gridKernel, 2, sizeof(cl_mem), &particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.gridKernel, 3, sizeof(cl_int), &numParticles);
	clRunKernel(openCLKernels.gridKernel, workSize, workgroupSize);

	// Compute prefix sum for grid
	clSetKernelArg(openCLKernels.prefixSumKernel, 2, sizeof(cl_uint), (void*)&gridElements);

	int pingpong = 0;
	for(cl_int offset = 1; offset <= gridElements; offset *= 2) {
		if(offset == 1) {
			clSetKernelArg(openCLKernels.prefixSumKernel, 0, sizeof(cl_mem), (void*)&particles.gridSizeBuffer);
		}
		else {
			clSetKernelArg(openCLKernels.prefixSumKernel, 0, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[0] : particles.gridBuffer[1]));
		}
		clSetKernelArg(openCLKernels.prefixSumKernel, 1, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[1] : particles.gridBuffer[0]));
		clSetKernelArg(openCLKernels.prefixSumKernel, 3, sizeof(cl_int), (void*)&offset);
		clRunKernel(openCLKernels.prefixSumKernel, gridWorkSize, workgroupSize);
		pingpong = 1 - pingpong;
	}

	// Reorganize particles 
	clSetKernelArg(openCLKernels.gridReorderKernel, 0, sizeof(cl_mem), &particles.particleBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 1, sizeof(cl_mem), &particles.particleBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 2, sizeof(cl_mem), &particles.velocityBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 3, sizeof(cl_mem), &particles.velocityBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 4, sizeof(cl_mem), &particles.offsetBuffer);
	clSetKernelArg(openCLKernels.gridReorderKernel, 5, sizeof(cl_mem), (void*)&particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.gridReorderKernel, 6, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[0] : particles.gridBuffer[1]));
	clSetKernelArg(openCLKernels.gridReorderKernel, 7, sizeof(cl_int), &numParticles);
	clRunKernel(openCLKernels.gridReorderKernel, workSize, workgroupSize);

	particle* testData = (particle*)malloc(sizeof(cl_float) * numParticles * 4);

	// Swap particle buffers
	particles.currentBuffer = 1 - particles.currentBuffer;

	// Send new cell select buffer
	int cellSelect[27];
	for(int i = 0; i < 27; i++) {
		cellSelect[i] = i;
	}
	shuffle(cellSelect, 27);
	clEnqueueWriteBuffer(clCommandQueue(), particles.cellSelectBuffer, true, 0, 27 * sizeof(cl_int), cellSelect, 0, 0, 0);
	clFinish(clCommandQueue());

	// Recalculate densities and normalized pressure derivatives
	clSetKernelArg(openCLKernels.dataKernel, 0, sizeof(cl_mem), &particles.particleBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.dataKernel, 1, sizeof(cl_mem), &particles.dataBuffer);

	clSetKernelArg(openCLKernels.dataKernel, 2, sizeof(cl_mem), (void*)&particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.dataKernel, 3, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[1] : particles.gridBuffer[0]));

	clSetKernelArg(openCLKernels.dataKernel, 4, sizeof(cl_mem), (void*)&particles.cellSelectBuffer);

	clSetKernelArg(openCLKernels.dataKernel, 5, sizeof(cl_int), &numParticles);
	clRunKernel(openCLKernels.dataKernel, workSize, workgroupSize);

	// Send new cell select buffer
	cellSelect[27];
	for(int i = 0; i < 27; i++) {
		cellSelect[i] = i;
	}
	shuffle(cellSelect, 27);
	clEnqueueWriteBuffer(clCommandQueue(), particles.cellSelectBuffer, true, 0, 27 * sizeof(cl_int), cellSelect, 0, 0, 0);
	clFinish(clCommandQueue());

	// Integrate position
	float dT = TIMESTEP;
	clSetKernelArg(openCLKernels.simulationKernel, 0, sizeof(cl_mem), &particles.particleBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.simulationKernel, 1, sizeof(cl_mem), &particles.particleBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.simulationKernel, 2, sizeof(cl_mem), &particles.velocityBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.simulationKernel, 3, sizeof(cl_mem), &particles.velocityBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.simulationKernel, 4, sizeof(cl_mem), &particles.dataBuffer);
	
	clSetKernelArg(openCLKernels.simulationKernel, 5, sizeof(cl_mem), (void*)&particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.simulationKernel, 6, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[1] : particles.gridBuffer[0]));

	clSetKernelArg(openCLKernels.simulationKernel, 7, sizeof(cl_mem), (void*)&particles.cellSelectBuffer);

	clSetKernelArg(openCLKernels.simulationKernel, 8, sizeof(cl_float), &dT);
	clSetKernelArg(openCLKernels.simulationKernel, 9, sizeof(cl_float), &time);
	clSetKernelArg(openCLKernels.simulationKernel, 10, sizeof(cl_int), &numParticles);

	clSetKernelArg(openCLKernels.simulationKernel, 11, sizeof(cl_mem), &particles.terrainBuffer);

	PaddedVector windDir = PadVector(
		TransformVector(YAxisRotationMatrix(particles.windAngle), MakeVector(1.0f, 0.0f, 0.0f))
	);
	clSetKernelArg(openCLKernels.simulationKernel, 12, sizeof(cl_float) * 4, &windDir);
	clSetKernelArg(openCLKernels.simulationKernel, 13, sizeof(cl_float), &particles.windPower);

	clRunKernel(openCLKernels.simulationKernel, workSize, workgroupSize);

	// Release buffers back to OpenGL
	releaseGLBuffer(particles.particleBuffer[particles.currentBuffer]);
	releaseGLBuffer(particles.particleBuffer[1 - particles.currentBuffer]);

	// Swap particle buffers
	particles.currentBuffer = 1 - particles.currentBuffer;

	//////////////////////// PART 2: RENDERIING ////////////////////////////////

	// Clear everything first thing.
	glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
	
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.backgroundFBO);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	
	// Activate shader.
	glUseProgram(objectShader.shaderProgram);

	// Set projection
	Matrix projection = PerspectiveMatrix(
		45.0f,
		(float)WINDOW_WIDTH/(float)WINDOW_HEIGHT,
		0.01f,
		100.0f
	);
	MatrixAsUniform(objectShader.projectionMatrix, projection);
	
	// Vertices
	glBindBuffer(GL_ARRAY_BUFFER, terrain.vertexBuffer);
	glVertexAttribPointer(
		objectShader.vertexPosition,
		4,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 4,
		(void*)0
	);
	glEnableVertexAttribArray(objectShader.vertexPosition);
	
	// Set modelview according to camera
	Matrix rot = lookatMatrix(camera.pos, VectorAdd(camera.pos, camera.front), camera.up);
	rot = MatrixMul(RotationMatrix(camera.elevation, MakeVector(1, 0, 0)), rot);
	Matrix trans = TranslationMatrix(-camera.pos.x, -camera.pos.y, -camera.pos.z);
	Matrix modelview =  MatrixMul(rot, trans);
	MatrixAsUniform(objectShader.modelviewMatrix, modelview);

	// Normal view matrix - inverse transpose of modelview.
	Matrix normalview = MatrixTranspose(FastMatrixInverse(modelview));
	MatrixAsUniform(objectShader.normalviewMatrix, normalview);

	// Set heightmap texture
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, terrain.heightTexture );
	glUniform1i(objectShader.terrainTexture, 0);

	// Turn off culling
	glDisable(GL_CULL_FACE);

	// Send element buffer to GPU and draw.
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, terrain.elementBuffer);

	glDrawElements(
		GL_TRIANGLES,
		512 * 512 * 6,
		GL_UNSIGNED_INT,
		(void*)0
	);
	
	// Turn culling back on
	glEnable(GL_CULL_FACE);

	// Switch to low-res viewport
	glViewport(0, 0, WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER);

	// Low-res projection matrix
	Matrix projectionLowres = PerspectiveMatrix(
		45.0f,
		(float)(WINDOW_WIDTH / RESOLUTION_DIVIDER) / (float)(WINDOW_HEIGHT / RESOLUTION_DIVIDER),
		0.01f,
		100.0f
	);

	// Activate particle depth FBO
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.particleFBO[0]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Activate shader
	glUseProgram(particleShader.shaderProgram);

	// Set textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, framebuffers.backgroundTexture);
	glUniform1i(particleShader.terrainTexture, 0);

	// Send uniforms
	MatrixAsUniform(particleShader.modelviewMatrix, modelview);
	MatrixAsUniform(particleShader.projectionMatrix, projectionLowres);
	glUniform2f(particleShader.screenSize, WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER);

	// Bind new buffer and set up arrtibutes
	glBindBuffer(GL_ARRAY_BUFFER, particles.vertexBuffer[particles.currentBuffer]);
	glVertexAttribPointer(
		particleShader.vertexPosition,
		4,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 4,
		(void*)0
	);
	glEnableVertexAttribArray(particleShader.vertexPosition);

	// Bind element buffer and draw particles
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, particles.elementBuffer);
	glDrawElements(
		GL_POINTS,
		NUM_PARTICLES,
		GL_UNSIGNED_INT,
		(void*)0
	);

	// Activate particle thickness FBO
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.particleThicknessFBO[0]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Activate shader
	glUseProgram(particleThicknessShader.shaderProgram);

	// Send uniforms
	MatrixAsUniform(particleThicknessShader.modelviewMatrix, modelview);
	MatrixAsUniform(particleThicknessShader.projectionMatrix, projectionLowres);
	glUniform2f(particleThicknessShader.screenSize, WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER);

	// Set textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, framebuffers.backgroundTexture);
	glUniform1i(particleThicknessShader.terrainTexture, 0);

	// Bind new buffer and set up arrtibutes
	glBindBuffer(GL_ARRAY_BUFFER, particles.vertexBuffer[particles.currentBuffer]);
	glVertexAttribPointer(
		particleThicknessShader.vertexPosition,
		4,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 4,
		(void*)0
	);
	glEnableVertexAttribArray(particleThicknessShader.vertexPosition);

	// Enable additive blending and disable depth test
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE);
	glDisable(GL_DEPTH_TEST);

	// Bind element buffer and draw particles, this time rendering thickness map
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, particles.elementBuffer);
	glDrawElements(
		GL_POINTS,
		NUM_PARTICLES,
		GL_UNSIGNED_INT,
		(void*)0
	);
	
	// Turn blending back off and depth test back on
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

	// Activate particle velocity FBO
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.velocityFBO);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Activate shader
	glUseProgram(particleVelocityShader.shaderProgram);

	// Send uniforms
	MatrixAsUniform(particleVelocityShader.modelviewMatrix, modelview);
	MatrixAsUniform(particleVelocityShader.projectionMatrix, projectionLowres);
	glUniform2f(particleVelocityShader.screenSize, WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER);

	// Bind new buffer and set up arrtibutes
	glBindBuffer(GL_ARRAY_BUFFER, particles.vertexBuffer[particles.currentBuffer]);
	glVertexAttribPointer(
		particleVelocityShader.vertexPosition,
		4,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 4,
		(void*)0
	);
	glEnableVertexAttribArray(particleVelocityShader.vertexPosition);

	// Bind element buffer and draw particles, this time rendering velocity map
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, particles.elementBuffer);
	glDrawElements(
		GL_POINTS,
		NUM_PARTICLES,
		GL_UNSIGNED_INT,
		(void*)0
	);

	// Curvature flow smoothing begins
	glUseProgram(curvatureFlowShader.shaderProgram);

	// Send uniforms
	glUniform1i(curvatureFlowShader.particleTexture, 0);

	glUniform2f(curvatureFlowShader.screenSize, WINDOW_WIDTH / RESOLUTION_DIVIDER, WINDOW_HEIGHT / RESOLUTION_DIVIDER);
	MatrixAsUniform(curvatureFlowShader.projectionMatrix, projectionLowres);

	// Prepare state
	glActiveTexture(GL_TEXTURE0);
	glBindBuffer(GL_ARRAY_BUFFER, screenQuad.vertexBuffer);
	glVertexAttribPointer(
		curvatureFlowShader.vertexPosition,
		3,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 3,
		(void*)0
	);
	glEnableVertexAttribArray(curvatureFlowShader.vertexPosition);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, screenQuad.elementBuffer);

	// Smoothing loop
	glDisable(GL_DEPTH_TEST);
	pingpong = 0;
	for(int i = 0; i < SMOOTHING_ITERATIONS; i++) {
		// Bind no FBO
		glBindFramebuffer(GL_FRAMEBUFFER, 0);

		// Bind texture
		glBindTexture(GL_TEXTURE_2D, framebuffers.particleTexture[pingpong]);

		// Activate proper FBO and clear
		glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.particleFBO[1 - pingpong]);

		// Draw a quad
		glDrawElements(
			GL_TRIANGLES,
			6,
			GL_UNSIGNED_INT,
			(void*)0
		);

		// Switch buffers
		pingpong = 1 - pingpong;
	}
	glEnable(GL_DEPTH_TEST);

	// Activate particle color FBO
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.particleColorFBO);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Liquid shading shader
	glUseProgram(liquidShadeShader.shaderProgram);

	// Bind and set textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, framebuffers.particleTexture[0]);
	glUniform1i(liquidShadeShader.particleTexture, 0);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, framebuffers.particleThicknessTexture[0]);
	glUniform1i(liquidShadeShader.particleThicknessTexture, 1);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, terrain.envTexture);
	glUniform1i(liquidShadeShader.environmentTexture, 2);

	glActiveTexture(GL_TEXTURE3);
	glBindTexture(GL_TEXTURE_2D, framebuffers.particleVelocityTexture);
	glUniform1i(liquidShadeShader.velocityTexture, 3);

	// Send uniforms
	glUniform2f(liquidShadeShader.screenSize, WINDOW_WIDTH, WINDOW_HEIGHT);
	MatrixAsUniform(liquidShadeShader.modelviewMatrix, modelview);
	MatrixAsUniform(liquidShadeShader.projectionMatrix, projection);

	// Draw a quad
    glDisable(GL_DEPTH_TEST);
	glBindBuffer(GL_ARRAY_BUFFER, screenQuad.vertexBuffer);
	glVertexAttribPointer(
		liquidShadeShader.vertexPosition,
		3,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 3,
		(void*)0
	);
	glEnableVertexAttribArray(liquidShadeShader.vertexPosition);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, screenQuad.elementBuffer);
	glDrawElements(
		GL_TRIANGLES,
		6,
		GL_UNSIGNED_INT,
		(void*)0
	);
	glEnable(GL_DEPTH_TEST);

	// Switch back to full-res viewport
	glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

	// Deactivate FBOs
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Compose shader
	glUseProgram(compositionShader.shaderProgram);

	// Bind and set textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, terrain.envTexture);
	glUniform1i(compositionShader.backgroundTexture, 0);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, framebuffers.particleColorTexture);
	glUniform1i(compositionShader.particleTexture, 1);

	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, framebuffers.backgroundTexture);
	glUniform1i(compositionShader.terrainTexture, 2);

	// Send uniforms
	MatrixAsUniform(compositionShader.modelviewMatrix, modelview);
	
	// Draw a quad
    glDisable(GL_DEPTH_TEST);
	glBindBuffer(GL_ARRAY_BUFFER, screenQuad.vertexBuffer);
	glVertexAttribPointer(
		compositionShader.vertexPosition,
		3,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 3,
		(void*)0
	);
	glEnableVertexAttribArray(compositionShader.vertexPosition);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, screenQuad.elementBuffer);
	glDrawElements(
		GL_TRIANGLES,
		6,
		GL_UNSIGNED_INT,
		(void*)0
	);
	glEnable(GL_DEPTH_TEST);

	// Switch drawing area and displayed area.
	glutSwapBuffers();
}

// Key press handler
// Mostly here to allow us to break on escape
void handleKeypress(unsigned char k, int x, int y) {
	switch(k) {
		case 27: // Escape -> die.
			exit(0);
		break;

		case 'y':
			particles.windAngle = particles.windAngle - 0.05f;
		break;

		case 'x':
			particles.windAngle = particles.windAngle + 3.14f;
		break;

		case 'c':
			particles.windAngle = particles.windAngle + 0.05f;
		break;

		case 'r':
			particles.windPower = max(0.0f, min(particles.windPower + 0.1f, 2.0f));
		break;

		case 'f':
			particles.windPower = max(0.0f, min(particles.windPower - 0.1f, 2.0f));
		break;

		case 'v':
			particles.windPower = 0.0f;
		break;

		case 'p': // Neat for debugging. Wait a second on 'p'.
			Sleep(1000);
		break;
	}
}

// Set things up and run
int main(int argc, char** argv) {
	makeWindow(argc, argv);
	acquireSharedOpenCLContext();
	initObjects();
	initPrograms();
	initFramebuffers();
	initShaders();
	setupGlutCallbacks();

	camera.pos = MakeVector(0, 0, -12);
	camera.up = MakeVector(0, 1, 0);
	camera.front = MakeVector(0, 0, 1);
	
	glutSetCursor(GLUT_CURSOR_NONE);
	glutMainLoop();

	releaseSharedOpenCLContext();
}