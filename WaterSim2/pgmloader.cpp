#include <stdio.h>  
#include <stdlib.h>  
#include <stdint.h>
#include <math.h>

#define min(a,b) ((a)<(b)?(a):(b))
#define max(a,b) ((a)>(b)?(a):(b))

// "THIS FUNCTION IS UNSAFE"
#pragma warning(disable: 4996)

uint8_t* readFileBytes(const char *name)  {  
	FILE *fl = fopen(name, "r");  
	fseek(fl, 0, SEEK_END);  
	long len = ftell(fl);  
	uint8_t *ret = (uint8_t*)malloc(len);  
	fseek(fl, 0, SEEK_SET);  
	fread(ret, 1, len, fl);  
	fclose(fl);  
	return ret;  
}  

float* loadPGM(const char* fileName, int w, int h) {
		uint8_t* inputData = readFileBytes(fileName);

        // Remove header
		int pos = 0;
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

        float* dataFloats = (float*)malloc(sizeof(float) * w * h);
        float minv = 100000000000000.0f;
        float maxv = 0.0f;
		
        for(int i = 0; i < w * h; i++) {
            uint8_t pixelData[2];
            pixelData[0] = inputData[2*i+1 + pos];
            pixelData[1] = inputData[2*i + pos];
            uint16_t* pixelDataProper = (uint16_t*)pixelData;
            dataFloats[i] = (float)pixelDataProper[0];
            minv = min(minv, dataFloats[i]);
            maxv = max(maxv, dataFloats[i]);
        }

        // Normalize
        for(int i = 0; i < w * h; i++) {
            dataFloats[i] = (dataFloats[i] - minv) / (maxv - minv);
        }
		
		free(inputData);

		return dataFloats;
}