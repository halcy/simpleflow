#include <stdio.h>  
#include <stdlib.h>  
#include <stdint.h>
#include <math.h>

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

// "THIS FUNCTION IS UNSAFE"
#pragma warning(disable: 4996)

uint8_t* readFileBytes(const char *name)  {  
	FILE *fl = fopen(name, "rb");  
	fseek(fl, 0, SEEK_END);  
	size_t len = ftell(fl);
	uint8_t *ret = (uint8_t*)malloc(len * sizeof(uint8_t));
	fseek(fl, 0, SEEK_SET);
    fread(ret, sizeof(uint8_t), len, fl);
	fclose(fl);
	return ret;  
}  

float* loadPGM(const char* fileName, int w, int h) {
		uint8_t* inputData = readFileBytes(fileName);

        // Remove header
		size_t pos = 0;
		while(inputData[pos] != '\n') {
			pos++;
		}
		pos++;

        // Remove width
		while(inputData[pos] != '\n') {
			pos++;
		}
		pos++;

        // Remove height
        while(inputData[pos] != '\n') {
			pos++;
		}
		pos++;

        // Remove rest
        while(inputData[pos] != '\n') {
			pos++;
		}
		pos++;
        printf("%d\n", pos);

        float* dataFloats = (float*)malloc(sizeof(float) * w * h);
        float minv = 10000000000.0f;
        float maxv = 0.0f;
		
        for(size_t i = 0; i < w * h; i++) {
            dataFloats[i] = (float)((inputData[2 * i + pos] << 8) + inputData[2 * i + 1 + pos]);
            minv = min(minv, dataFloats[i]);
            maxv = max(maxv, dataFloats[i]);
        }

        // Normalize and fit into range we want
        for(int i = 0; i < w * h; i++) {
            dataFloats[i] = ((dataFloats[i] - minv) / (maxv - minv)) * 0.07;
        }
		
		free(inputData);

		return dataFloats;
}