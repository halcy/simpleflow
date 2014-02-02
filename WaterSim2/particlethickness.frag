#version 420 core

// Particle depth fragment shader

// Parameters from the vertex shader
in vec3 eyespacePos;
in float eyespaceRadius;

// Uniforms
uniform mat4 projection;
uniform vec2 screenSize;

// Output
out float particleThickness;

void main() {
	vec3 normal;

	// See where we are inside the point sprite
	normal.xy = (gl_PointCoord - 0.5f) * 2.0f;
	float dist = length(normal);
	
	// Outside sphere? Discard.
	if(dist > 1.0f) {
		discard;
	}

	// Additive-Blend to get a sort-of thickness.
	// TODO maybe a constant value isn't great here...
	// particleThickness = 1.0f;
	particleThickness = 1.0f - length(normal);
}
