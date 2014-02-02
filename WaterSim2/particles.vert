#version 420 core

// Particle vertex shader

// Attributes: Position only, since these are simply particles
in vec4 vertex;

// Uniforms
uniform mat4 projection;
uniform mat4 modelview;
uniform vec2 screenSize;

// Parameters passed to the fragment shader.
out vec3 eyespacePos;
out float eyespaceRadius;

void main() {
	// Transform
	vec4 eyespacePos4 = vec4(vertex.xyz, 1.0f) * modelview;
	eyespacePos = eyespacePos4.xyz;
	eyespaceRadius = 1.0f / (-eyespacePos.z * 4.0f * (1.0f / screenSize.y));
	vec4 clipspacePos = eyespacePos4 * projection;

	// Set up variables for rasterizer
	gl_Position = clipspacePos;
	gl_PointSize = eyespaceRadius;
}
