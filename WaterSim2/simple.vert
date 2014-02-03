#version 420 core

// Simple vertex shader.
// Transforms and projects vertices and calculates parameters for lighting.

// Attributes: Position, normal, texture coordinates
in vec4 vertex;

// Same for the whole model or scene: Projection and Modelview matrices
uniform mat4 projection;
uniform mat4 modelview;
uniform mat4 normalview;

// To fragment shader
out vec3 objectPos;

void main() {
	// Transform the vertex according to modelview
	vec4 viewVertex;	
	viewVertex = vertex * modelview;

	// Project and send to the fragment shader
	gl_Position = viewVertex * projection;
	objectPos = vertex.xyz;
}
