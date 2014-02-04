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
out vec3 worldPos;

void main() {
	objectPos = vertex.xyz;

	// Transform the vertex according to modelview
	vec4 viewVertex;	
	viewVertex = vertex * modelview;
	worldPos = viewVertex.xyz;

	// Project and send to the fragment shader
	gl_Position = viewVertex * projection;
}
