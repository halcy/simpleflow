#ifndef __DATA_H__
#define __DATA_H__

// Buffer objects for OpenGL and OpenCL
struct {
	int currentBuffer;
	GLuint vertexBuffer[2];
	
	cl_mem particleBuffer[2];
	cl_mem velocityBuffer[2];
	cl_mem dataBuffer;
	cl_mem gridSizeBuffer;
	cl_mem gridBuffer[2];
	cl_mem offsetBuffer;
	cl_mem cellSelectBuffer;

	cl_mem terrainBuffer;

	GLuint elementBuffer;

	float windAngle;
	float windPower;
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

	GLuint velocityFBO;
	GLuint particleVelocityTexture;

	GLuint particleThicknessFBO[2];
	GLuint particleThicknessTexture[2];

	GLuint particleColorFBO;
	GLuint particleColorTexture;
} framebuffers;

// A bunch of shaders
struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint modelviewMatrix;
	GLuint normalviewMatrix;
	GLuint projectionMatrix;

	GLuint terrainTexture;
} objectShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint modelviewMatrix;
	GLuint projectionMatrix;
	GLuint screenSize;

	GLuint terrainTexture;
} particleShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint modelviewMatrix;
	GLuint projectionMatrix;
	GLuint screenSize;

	GLuint terrainTexture;
} particleThicknessShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint modelviewMatrix;
	GLuint projectionMatrix;
	GLuint screenSize;
} particleVelocityShader;

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
	GLuint modelviewMatrix;
	GLuint screenSize;

	GLuint environmentTexture;
	GLuint particleTexture;
	GLuint particleThicknessTexture;
	GLuint velocityTexture;
} liquidShadeShader;

struct {
	GLuint shaderProgram;

	GLuint vertexPosition;

	GLuint modelviewMatrix;

	GLuint backgroundTexture;
	GLuint particleTexture;
	GLuint terrainTexture;
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

// Environment info
struct {
	Vector pos;
	Vector front;
	Vector up;
	float elevation;
} camera;

// A terrain
struct {
	float* heightData;
	GLuint heightTexture;

	GLuint vertexBuffer;
	GLuint elementBuffer;

	GLuint envTexture;
} terrain;

#endif