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

// Include CL/GL interop tools.
#include "clgl.h"

// Buffer objects
struct {
	GLuint vertexBuffer;
	GLuint elementBuffer;
} cube;

struct {
	int currentBuffer;
	GLuint vertexBuffer[2];
	
	cl_mem particleBuffer[2];
	cl_mem velocityBuffer[2];
	cl_mem dataBuffer[2];
	cl_mem gridSizeBuffer;
	cl_mem gridBuffer[2];

	GLuint elementBuffer;
} particles;

struct {
	GLuint vertexBuffer;
	GLuint elementBuffer;
} screenQuad;

struct {
	GLuint backgroundFBO;
	GLuint backgroundTexture;

	GLuint particleFBO[2];
	GLuint particleTexture[2];

	GLuint particleThicknessFBO[2];
	GLuint particleThicknessTexture[2];

	GLuint particleColorFBO;
	GLuint particleColorTexture;
} framebuffers;

// A bunch of shaders
struct {
	GLuint shaderProgram;

	GLuint vertexPosition;
	GLuint vertexTexcoords;
	GLuint vertexNormal;

	GLuint modelviewMatrix;
	GLuint normalviewMatrix;
	GLuint projectionMatrix;

	GLuint texture;
} objectShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint modelviewMatrix;
	GLuint projectionMatrix;
	GLuint screenSize;
} particleShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint modelviewMatrix;
	GLuint projectionMatrix;
	GLuint screenSize;
} particleThicknessShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint projectionMatrix;
	GLuint screenSize;

	GLuint particleTexture;
} curvatureFlowShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint projectionMatrix;
	GLuint screenSize;

	GLuint environmentTexture;
	GLuint particleTexture;
	GLuint particleThicknessTexture;
} liquidShadeShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint backgroundTexture;
	GLuint particleTexture;
} compositionShader;

// A bunch of kernels
struct {
	cl_program program;

	cl_kernel gridClearKernel;
	cl_kernel gridKernel;
	cl_kernel prefixSumKernel;
	cl_kernel gridReorderKernel;
	cl_kernel dataKernel;
	cl_kernel simulationKernel;
} openCLKernels;


// Textures
GLuint defaultTextureData;

// Vertex array object
GLuint vertexArray;

// Status variable
float angle;

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

#include "cubeverts.h"

// [-0.5, 0.5] RNG
float centeredUnitRand() {
	return (((float)rand()/(float)RAND_MAX) - 0.5f);
}

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
}

// Initialize buffers, textures, etc.
void initObjects() {
	// Create a cube VBO
	cube.vertexBuffer = makeBO(
		GL_ARRAY_BUFFER,
		cubeVertices,
		sizeof(cubeVertices),
		GL_STATIC_DRAW
	);
	cube.elementBuffer = makeBO(
		GL_ELEMENT_ARRAY_BUFFER,
		cubeElements,
		sizeof(cubeElements),
		GL_STATIC_DRAW
	);

	// Put a bunch of particles into space
	particle* particleData = (particle*)malloc(sizeof(particle) * NUM_PARTICLES);
	GLushort* particleElements = (GLushort*)malloc(sizeof(GLushort) * NUM_PARTICLES);

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
		sizeof(GLushort) * NUM_PARTICLES,
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
	GLushort quadElements[] = {0, 1, 3, 1, 2, 3};

	screenQuad.vertexBuffer = makeBO(
		GL_ARRAY_BUFFER,
		quadData,
		sizeof(Vector) * 4,
		GL_STATIC_DRAW
	);
	screenQuad.elementBuffer = makeBO(
		GL_ELEMENT_ARRAY_BUFFER,
		quadElements,
		sizeof(GLushort) * 6,
		GL_STATIC_DRAW
	);

	// Load textures.
	defaultTextureData = loadTexture("skymap_b.tga");

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
		particles.dataBuffer[i] = clCreateBuffer(
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
	particles.gridSizeBuffer = clCreateBuffer(
		clContext(), 
		CL_MEM_READ_WRITE | CL_MEM_COPY_HOST_PTR, 
		sizeof(cl_int) * GRID_SIZE * GRID_SIZE * GRID_SIZE, 
		zeroGrid, 
		NULL
	);
	
	free(zeroData);

	// Create an OpenCL program and kernel
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
	objectShader.vertexNormal = glGetAttribLocation(objectShader.shaderProgram, "vnormal");
	objectShader.vertexTexcoords = glGetAttribLocation(objectShader.shaderProgram, "texcoords");

	objectShader.modelviewMatrix = glGetUniformLocation(objectShader.shaderProgram, "modelview");
	objectShader.normalviewMatrix = glGetUniformLocation(objectShader.shaderProgram, "normalview");
	objectShader.projectionMatrix = glGetUniformLocation(objectShader.shaderProgram, "projection");
	objectShader.texture = glGetUniformLocation(objectShader.shaderProgram, "texture");

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

	// Bind output variables
	glBindFragDataLocation(particleThicknessShader.shaderProgram, 0, "particleThickness");

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
	liquidShadeShader.projectionMatrix = glGetUniformLocation(liquidShadeShader.shaderProgram, "projection");
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
	compositionShader.particleTexture = glGetUniformLocation(compositionShader.shaderProgram, "particleTexture");

	// Output
	glBindFragDataLocation(compositionShader.shaderProgram, 0, "outColor");
}

// Update status variables.
// Called every 15ms, unless the PC is too slow.
void updateI(int) { update(); }

void update() {
	// Update things here.
	angle += 0.02f;

	// Redraw screen now, please, and call again in 15ms.
	glutPostRedisplay();
	glutTimerFunc(15, updateI, 0);
	glutKeyboardFunc(handleKeypress);
}

// Draw the scene to the screen
void draw() {
	// Grab buffers for OpenCL
	acquireGLBuffer(particles.particleBuffer[particles.currentBuffer]);
	acquireGLBuffer(particles.particleBuffer[1 - particles.currentBuffer]);

	// Prepare to run some kernels
	size_t workSize[3] = {NUM_PARTICLES, 0, 0};
	size_t workgroupSize[3] = {1024, 0, 0};
	int numParticles = NUM_PARTICLES;
	int gridElements = GRID_SIZE * GRID_SIZE * GRID_SIZE;

	// Clear grid
	cl_uint clearGlobalWorkSize[] = {gridElements, 0, 0};
	cl_uint clearLocalWorkSize[] = {min(workgroupSize[0], clearGlobalWorkSize[0]), 0, 0};
	clSetKernelArg(openCLKernels.gridClearKernel, 0, sizeof(cl_mem), &particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.gridClearKernel, 1, sizeof(cl_int), &gridElements);
	clRunKernel(openCLKernels.gridClearKernel, clearGlobalWorkSize, clearLocalWorkSize);

	// Compute grid positions
	clSetKernelArg(openCLKernels.gridKernel, 0, sizeof(cl_mem), &particles.particleBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridKernel, 1, sizeof(cl_mem), &particles.dataBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridKernel, 2, sizeof(cl_mem), &particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.gridKernel, 3, sizeof(cl_int), &numParticles);
	clRunKernel(openCLKernels.gridKernel, workSize, workgroupSize);

	// Compute prefix sum for grid
	clSetKernelArg(openCLKernels.prefixSumKernel, 2, sizeof(cl_uint), (void*)&gridElements);

	int pingpong = 0;
	for(cl_uint offset = 1; offset <= gridElements; offset *= 2) {
		cl_uint globalWorkSize[] = {gridElements, 0, 0};
		cl_uint localWorkSize[] = {min(workgroupSize[0], globalWorkSize[0]), 0, 0};
		if(offset == 1) {
			clSetKernelArg(openCLKernels.prefixSumKernel, 0, sizeof(cl_mem), (void*)&particles.gridSizeBuffer);
		}
		else {
			clSetKernelArg(openCLKernels.prefixSumKernel, 0, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[0] : particles.gridBuffer[1]));
		}
		clSetKernelArg(openCLKernels.prefixSumKernel, 1, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[1] : particles.gridBuffer[0]));
		clSetKernelArg(openCLKernels.prefixSumKernel, 3, sizeof(cl_int), (void*)&offset);
		clRunKernel(openCLKernels.prefixSumKernel, globalWorkSize, localWorkSize);
		pingpong = 1 - pingpong;
	}

	// Reorganize particles
	clSetKernelArg(openCLKernels.gridReorderKernel, 0, sizeof(cl_mem), &particles.particleBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 1, sizeof(cl_mem), &particles.particleBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 2, sizeof(cl_mem), &particles.velocityBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 3, sizeof(cl_mem), &particles.velocityBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 4, sizeof(cl_mem), &particles.dataBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 5, sizeof(cl_mem), &particles.dataBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.gridReorderKernel, 6, sizeof(cl_mem), (void*)&particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.gridReorderKernel, 7, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[0] : particles.gridBuffer[1]));
	clSetKernelArg(openCLKernels.gridReorderKernel, 8, sizeof(cl_int), &numParticles);
	clRunKernel(openCLKernels.gridReorderKernel, workSize, workgroupSize);

	// Swap particle buffers
	particles.currentBuffer = 1 - particles.currentBuffer;

	// Recalculate densities and normalized pressure derivatives
	clSetKernelArg(openCLKernels.dataKernel, 0, sizeof(cl_mem), &particles.particleBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.dataKernel, 1, sizeof(cl_mem), &particles.dataBuffer[particles.currentBuffer]);

	clSetKernelArg(openCLKernels.dataKernel, 2, sizeof(cl_mem), (void*)&particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.dataKernel, 3, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[1] : particles.gridBuffer[0]));

	clSetKernelArg(openCLKernels.dataKernel, 4, sizeof(cl_int), &numParticles);
	clRunKernel(openCLKernels.dataKernel, workSize, workgroupSize);

	// Integrate position
	float dT = 0.01f;
	clSetKernelArg(openCLKernels.simulationKernel, 0, sizeof(cl_mem), &particles.particleBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.simulationKernel, 1, sizeof(cl_mem), &particles.particleBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.simulationKernel, 2, sizeof(cl_mem), &particles.velocityBuffer[particles.currentBuffer]);
	clSetKernelArg(openCLKernels.simulationKernel, 3, sizeof(cl_mem), &particles.velocityBuffer[1 - particles.currentBuffer]);
	clSetKernelArg(openCLKernels.simulationKernel, 4, sizeof(cl_mem), &particles.dataBuffer[particles.currentBuffer]);
	
	clSetKernelArg(openCLKernels.simulationKernel, 5, sizeof(cl_mem), (void*)&particles.gridSizeBuffer);
	clSetKernelArg(openCLKernels.simulationKernel, 6, sizeof(cl_mem), (void*)&(pingpong == 0 ? particles.gridBuffer[1] : particles.gridBuffer[0]));

	clSetKernelArg(openCLKernels.simulationKernel, 7, sizeof(cl_float), &dT);
	clSetKernelArg(openCLKernels.simulationKernel, 8, sizeof(cl_float), &angle);
	clSetKernelArg(openCLKernels.simulationKernel, 9, sizeof(cl_int), &numParticles);
	clRunKernel(openCLKernels.simulationKernel, workSize, workgroupSize);

	// Release buffers back to OpenGL
	releaseGLBuffer(particles.particleBuffer[particles.currentBuffer]);
	releaseGLBuffer(particles.particleBuffer[1 - particles.currentBuffer]);

	// Swap particle buffers
	particles.currentBuffer = 1 - particles.currentBuffer;

	// Clear everything first thing.
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

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
	glBindBuffer(GL_ARRAY_BUFFER, cube.vertexBuffer);
	glVertexAttribPointer(
		objectShader.vertexPosition,
		4,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 9,
		(void*)0
	);
	glEnableVertexAttribArray(objectShader.vertexPosition);

	// Normals
	glVertexAttribPointer(
		objectShader.vertexNormal,
		3,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 9,
		(void*)(sizeof(GLfloat) * 4)
	);
	glEnableVertexAttribArray(objectShader.vertexNormal);
	
	// Texture coordinates
	glVertexAttribPointer(
		objectShader.vertexTexcoords,
		2,
		GL_FLOAT,
		GL_FALSE,
		sizeof(GLfloat) * 9,
		(void*)(sizeof(GLfloat) * 7)
	);
	glEnableVertexAttribArray(objectShader.vertexTexcoords);
	
	// Set modelview - the objects rotation and translation.
	/*Matrix rot = XAxisRotationMatrix(angle);
	rot = MatrixMul(rot, ZAxisRotationMatrix(angle / 2.0f));
	Matrix trans = TranslationMatrix(0.0f, 0.0f, -5.0f);
	Matrix modelview = MatrixMul(trans, rot);
	MatrixAsUniform(objectShader.modelviewMatrix, modelview);*/

	Matrix rot = XAxisRotationMatrix(0.5f);
	rot = MatrixMul(rot, YAxisRotationMatrix(3.14f / 4.0f));
	Matrix trans = TranslationMatrix(0.0f, 2.0f, -12.0f);
	Matrix modelview = MatrixMul(trans, rot);
	MatrixAsUniform(objectShader.modelviewMatrix, modelview);

	// Normal view matrix - inverse transpose of modelview.
	Matrix normalview = MatrixTranspose(FastMatrixInverse(modelview));
	MatrixAsUniform(objectShader.normalviewMatrix, normalview);

	// Textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, defaultTextureData);
	glUniform1i(objectShader.texture, 0);

	// Send element buffer to GPU and draw.
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube.elementBuffer);

	glDrawElements(
		GL_TRIANGLES,
		36,
		GL_UNSIGNED_SHORT,
		(void*)0
	);
	
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
	glBindFramebuffer(GL_FRAMEBUFFER, framebuffers.particleFBO[0]);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Activate shader
	glUseProgram(particleShader.shaderProgram);

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
		GL_UNSIGNED_SHORT,
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
		GL_UNSIGNED_SHORT,
		(void*)0
	);
	
	// Turn blending back off and depth test back on
	glDisable(GL_BLEND);
	glEnable(GL_DEPTH_TEST);

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
			GL_UNSIGNED_SHORT,
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
	glBindTexture(GL_TEXTURE_2D, defaultTextureData);
	glUniform1i(liquidShadeShader.environmentTexture, 2);

	// Send uniforms
	glUniform2f(liquidShadeShader.screenSize, WINDOW_WIDTH, WINDOW_HEIGHT);
	MatrixAsUniform(liquidShadeShader.projectionMatrix, projection);

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
	glEnableVertexAttribArray(liquidShadeShader.vertexPosition);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, screenQuad.elementBuffer);
	glDrawElements(
		GL_TRIANGLES,
		6,
		GL_UNSIGNED_SHORT,
		(void*)0
	);
	glEnable(GL_DEPTH_TEST);

	// Switch back to full-res viewport
	glViewport(0, 0, WINDOW_WIDTH, WINDOW_HEIGHT);

	// Deactivate FBOs
	glBindFramebuffer(GL_FRAMEBUFFER, 0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Liquid shading shader
	glUseProgram(compositionShader.shaderProgram);

	// Bind and set textures
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, defaultTextureData);
	glUniform1i(compositionShader.backgroundTexture, 0);

	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, framebuffers.particleColorTexture);
	glUniform1i(compositionShader.particleTexture, 1);

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
		GL_UNSIGNED_SHORT,
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
	initFramebuffers();
	initShaders();
	setupGlutCallbacks();
	glutMainLoop();
	releaseSharedOpenCLContext();
}