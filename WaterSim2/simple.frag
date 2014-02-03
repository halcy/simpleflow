#version 420 core

// Simple fragment shader.
// Does texturing and phong shading.

// From vertex shader
in vec3 objectPos;

// Textures
uniform sampler2D texture;

// Output
out vec4 outColor;

void main() {
	float heightcol = (objectPos.y + 5.0f) / 8.0f;
	outColor = vec4(heightcol, heightcol, heightcol, gl_FragCoord.z);
}
