#version 420 core

// Particle velocity fragment shader

// Parameters from the vertex shader
in vec3 eyespacePos;
in float eyespaceRadius;
in float velocity;

// Uniforms
uniform mat4 projection;
uniform vec2 screenSize;

// Output
out float velocityMap;

void main() {
	vec3 normal;

	// See where we are inside the point sprite
	normal.xy = (gl_PointCoord - 0.5f) * 2.0f;
	float dist = length(normal);
	
	// Outside sphere? Discard.
	if(dist > 1.0f) {
		discard;
	}

	velocityMap = velocity;
}
