/**
 * OpenGL / OpenCL interop code
 */

#include "clgl.h"
#include "glhelpers.h"

// "THIS FUNCTION IS UNSAFE"
#pragma warning(disable: 4996)

struct {
	cl_platform_id platformId;
	cl_device_id deviceId;
	cl_context context;
	cl_command_queue queue;
} openCLState;


void __stdcall openCLError(const char *errinfo, const void *private_info, size_t cb, void *user_data) {
	printf("OpenCL error: %s\n", errinfo);
	getc(stdin);
}

void acquireSharedOpenCLContext() {
	// Error code
	cl_int clError;

	// Get platform ID(s)
	cl_uint platform_count;
	clGetPlatformIDs(0, NULL, &platform_count);
	cl_platform_id* platform_ids = new cl_platform_id[platform_count];
	clGetPlatformIDs(platform_count, platform_ids, NULL);

	//get a reference to the first available GPU device
	for(unsigned int i = 0; i < platform_count; i++) {
		cl_uint device_count = 0;
		clError = clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_GPU, 0, NULL, &device_count);
		if(clError == CL_DEVICE_NOT_FOUND || device_count == 0) {
			continue;
		}

		cl_device_id devices[1];
		clGetDeviceIDs(platform_ids[i], CL_DEVICE_TYPE_GPU, 1, devices, NULL);
		openCLState.platformId = platform_ids[i];
		openCLState.deviceId = devices[0];
		break;
	}

	// Create a new OpenCL context on the selected device which supports sharing with OpenGL

   // Define OS-specific context properties and create the OpenCL context
    #if defined (__APPLE__) // Mac OS
        CGLContextObj kCGLContext = CGLGetCurrentContext();
        CGLShareGroupObj kCGLShareGroup = CGLGetShareGroup(kCGLContext);
        cl_context_properties props[] = 
        {
            CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE, (cl_context_properties)kCGLShareGroup, 
            0 
        };
        g_CLContext = clCreateContext(props, 0,0, NULL, NULL, &clError);
    #else
        #ifndef _WIN32 // Lunix
            cl_context_properties props[] = 
            {
                CL_GL_CONTEXT_KHR, (cl_context_properties)glXGetCurrentContext(), 
                CL_GLX_DISPLAY_KHR, (cl_context_properties)glXGetCurrentDisplay(), 
                CL_CONTEXT_PLATFORM, (cl_context_properties)openCLState.platformId, 
                0
            };
            g_CLContext = clCreateContext(props, 1, &openCLState.deviceId, NULL, NULL, &clError);
        #else // Win32
            cl_context_properties props[] = 
            {
                CL_GL_CONTEXT_KHR, (cl_context_properties)wglGetCurrentContext(), 
                CL_WGL_HDC_KHR, (cl_context_properties)wglGetCurrentDC(), 
				CL_CONTEXT_PLATFORM, (cl_context_properties)openCLState.platformId, 
                0
            };
			openCLState.context = clCreateContext(props, 1, &openCLState.deviceId, &openCLError, NULL, &clError);
        #endif
    #endif

	// Create command queue
	openCLState.queue = clCreateCommandQueue(openCLState.context, openCLState.deviceId, 0, &clError);

	// If debug info is desired, print some useful information
	#ifdef DEBUG
		cl_ulong mem_size;
		clGetDeviceInfo (openCLState.deviceId, CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
		printf("Global memory: %u kb\n", (unsigned int)(mem_size/1024));
	
		clGetDeviceInfo (openCLState.deviceId, CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL);
		printf("Local memory: %u kb\n", (unsigned int)(mem_size/1024));
	
		clGetDeviceInfo (openCLState.deviceId, CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE, sizeof(mem_size), &mem_size, NULL);
		printf("Constant memory: %u kb\n",  (unsigned int)(mem_size/1024));
	
		size_t workgroup_sizes[3];
		clGetDeviceInfo (openCLState.deviceId, CL_DEVICE_MAX_WORK_ITEM_SIZES, sizeof(size_t) * 3, &workgroup_sizes, NULL);
		printf("Workgroup sizes max: (%u, %u, %u)\n", workgroup_sizes[0], workgroup_sizes[1], workgroup_sizes[2]);

		clGetDeviceInfo (openCLState.deviceId, 0x4003, sizeof(mem_size), &mem_size, NULL);
		printf("(NVIDIA) Warp size: %u\n",  (unsigned int)(mem_size));
	#endif

	delete platform_ids;
}

void releaseSharedOpenCLContext() {
	SAFE_RELEASE(clReleaseCommandQueue, openCLState.queue);
	SAFE_RELEASE(clReleaseContext, openCLState.context);
}

cl_mem sharedBuffer(GLuint buffer, cl_mem_flags accessFlags) {
	return clCreateFromGLBuffer(openCLState.context, accessFlags, buffer, NULL);
}

void acquireGLBuffer(cl_mem buffer) {
	clEnqueueAcquireGLObjects(openCLState.queue, 1, &buffer, 0, NULL, NULL);
}

void releaseGLBuffer(cl_mem buffer) {
	clEnqueueReleaseGLObjects(openCLState.queue, 1, &buffer, 0, NULL, NULL);
	clFinish(openCLState.queue);
}

cl_command_queue clCommandQueue() {
	return openCLState.queue;
}

cl_context clContext() {
	return openCLState.context;
}

cl_program clProgramFromFile(char* fileName) {
	char *programSrc = loadFile(fileName);
	size_t programSrcLength = strlen(programSrc);
	cl_program program = clCreateProgramWithSource(openCLState.context, 1, (const char**)&programSrc, &programSrcLength, NULL);
	//char* opts = "-Werror -cl-opt-disable";
	//char* opts = "-Werror";
	char* opts = "-cl-mad-enable -cl-no-signed-zeros -cl-unsafe-math-optimizations -cl-finite-math-only -cl-fast-relaxed-math -cl-strict-aliasing";
	if(clBuildProgram(program, 1, &openCLState.deviceId, opts, NULL, NULL) != CL_SUCCESS) {
		size_t logSize;
		clGetProgramBuildInfo(program, openCLState.deviceId, CL_PROGRAM_BUILD_LOG, 0, NULL, &logSize);

		char* buildLog = (char*)malloc(sizeof(char) * (logSize + 1));
		clGetProgramBuildInfo(program, openCLState.deviceId, CL_PROGRAM_BUILD_LOG, logSize, buildLog, NULL);
		buildLog[logSize] = '\0';

		printf("Errors building program %s:\n%s\n\n", fileName, buildLog);
		free(buildLog);

		getc(stdin);
		exit(1);
	}

	return program;
}

size_t adjustWorkSize(size_t workSize, size_t workgroupSize) {
	size_t r = workSize % workgroupSize;
	if(r == 0) {
		return workSize;
	} 
	else {
		return workSize + workgroupSize - r;
	}
}

static const char* errorToString(cl_int error) {
	#define CL_ERROR(x) case (x): return #x;
	switch(error) {
		CL_ERROR(CL_SUCCESS);
		CL_ERROR(CL_DEVICE_NOT_FOUND);
		CL_ERROR(CL_DEVICE_NOT_AVAILABLE);
		CL_ERROR(CL_COMPILER_NOT_AVAILABLE);
		CL_ERROR(CL_MEM_OBJECT_ALLOCATION_FAILURE);
		CL_ERROR(CL_OUT_OF_RESOURCES);
		CL_ERROR(CL_OUT_OF_HOST_MEMORY);
		CL_ERROR(CL_PROFILING_INFO_NOT_AVAILABLE);
		CL_ERROR(CL_MEM_COPY_OVERLAP);
		CL_ERROR(CL_IMAGE_FORMAT_MISMATCH);
		CL_ERROR(CL_IMAGE_FORMAT_NOT_SUPPORTED);
		CL_ERROR(CL_BUILD_PROGRAM_FAILURE);
		CL_ERROR(CL_MAP_FAILURE);
		CL_ERROR(CL_INVALID_VALUE);
		CL_ERROR(CL_INVALID_DEVICE_TYPE);
		CL_ERROR(CL_INVALID_PLATFORM);
		CL_ERROR(CL_INVALID_DEVICE);
		CL_ERROR(CL_INVALID_CONTEXT);
		CL_ERROR(CL_INVALID_QUEUE_PROPERTIES);
		CL_ERROR(CL_INVALID_COMMAND_QUEUE);
		CL_ERROR(CL_INVALID_HOST_PTR);
		CL_ERROR(CL_INVALID_MEM_OBJECT);
		CL_ERROR(CL_INVALID_IMAGE_FORMAT_DESCRIPTOR);
		CL_ERROR(CL_INVALID_IMAGE_SIZE);
		CL_ERROR(CL_INVALID_SAMPLER);
		CL_ERROR(CL_INVALID_BINARY);
		CL_ERROR(CL_INVALID_BUILD_OPTIONS);
		CL_ERROR(CL_INVALID_PROGRAM);
		CL_ERROR(CL_INVALID_PROGRAM_EXECUTABLE);
		CL_ERROR(CL_INVALID_KERNEL_NAME);
		CL_ERROR(CL_INVALID_KERNEL_DEFINITION);
		CL_ERROR(CL_INVALID_KERNEL);
		CL_ERROR(CL_INVALID_ARG_INDEX);
		CL_ERROR(CL_INVALID_ARG_VALUE);
		CL_ERROR(CL_INVALID_ARG_SIZE);
		CL_ERROR(CL_INVALID_KERNEL_ARGS);
		CL_ERROR(CL_INVALID_WORK_DIMENSION);
		CL_ERROR(CL_INVALID_WORK_GROUP_SIZE);
		CL_ERROR(CL_INVALID_WORK_ITEM_SIZE);
		CL_ERROR(CL_INVALID_GLOBAL_OFFSET);
		CL_ERROR(CL_INVALID_EVENT_WAIT_LIST);
		CL_ERROR(CL_INVALID_EVENT);
		CL_ERROR(CL_INVALID_OPERATION);
		CL_ERROR(CL_INVALID_GL_OBJECT);
		CL_ERROR(CL_INVALID_BUFFER_SIZE);
		CL_ERROR(CL_INVALID_MIP_LEVEL);
		default:
			return "Unknown error code";
	}
	#undef CL_ERROR
}

void clRunKernel(cl_kernel kernel, const size_t minWorkSize[3], const size_t workgroupSize[3]) {
	cl_int dimensions = 3;
	if(minWorkSize[2] == 0) {
		dimensions = 2;
	}
	if(minWorkSize[1] == 0) {
		dimensions = 1;
	}

	size_t workSize[3];
	/*for(int i = 0; i < dimensions; i++) {
		workSize[i] = adjustWorkSize(minWorkSize[i], workgroupSize[i]);
	}*/

	cl_int clErr = clEnqueueNDRangeKernel(openCLState.queue, kernel, dimensions, NULL, minWorkSize, workgroupSize, 0, NULL, NULL);
	if(clErr != CL_SUCCESS) {
		printf("Failed to run kernel: %s\n", errorToString(clErr));
		fgetc(stdin);
	}
	clErr = clFinish(openCLState.queue);
	if(clErr != CL_SUCCESS) {
		printf("Failed to finish kernel: %s\n", errorToString(clErr));
		fgetc(stdin);
	}
}
