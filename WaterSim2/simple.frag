#version 420 core

// Simple fragment shader.
// Does texturing and phong shading.

// From vertex shader
in vec3 objectPos;
in vec3 worldPos;

// Textures
uniform sampler2D terrainTexture;

// Output
out vec4 outColor;

float heightmap(vec2 coords) {
	//return texture(terrainTexture, coords) * 128.0f - 4.0f;
	return texture(terrainTexture, coords) * 32.0f - 1.0f;
}

vec2 heightmapCoords(vec2 coordsIn) {
	return ((coordsIn + vec2(10.5f)) / (10.5f * 2.0f)) / 8.0f;
}

vec3 grad(vec2 pos) {
	float off = 1.0f / 4096f;
	float yx = heightmap(pos + vec2(off, 0.0f));
	float yz = heightmap(pos + vec2(0.0f, off));
	float yxn = heightmap(pos - vec2(off, 0.0f));
	float yzn = heightmap(pos - vec2(0.0f, off));
	vec3 vx = normalize(vec3(pos.x + off, yx, pos.y) - vec3(pos.x - off, yxn, pos.y));
	vec3 vz = normalize(vec3(pos.x, yx, pos.y + off) - vec3(pos.x, yzn, pos.y - off));
	return cross(vx, vz);
}

void main() {
	vec2 heightcoords = heightmapCoords(objectPos.xz);
	float heightcol = heightmap(heightcoords);

	vec3 normal = normalize(grad(heightcoords));
	normal.y = -normal.y;
	
	// Sun is at infinity
	vec3 lightDir = vec3(1.0f, 0.0f, -1.0f);
	float lambert = max(0.0f, dot(normalize(lightDir), normal));
	float light = (lambert + 0.6f) * 0.5f;
	
	outColor = vec4(light, light, light, gl_FragCoord.z);
}
