#version 420 core

// Simple vertex shader.
// Transforms and projects vertices and calculates parameters for lighting.

// Attributes: Position, normal, texture coordinates
in vec4 vertex;
in vec3 vnormal;
in vec2 texcoords;

// Same for the whole model or scene: Projection and Modelview matrices
uniform mat4 projection;
uniform mat4 modelview;
uniform mat4 normalview;

// Parameters passed to the fragment shader.
out vec2 texcoord;
out vec3 eye;
out vec3 light;
out vec3 normal;

void main() {
	vec3 lightPos = vec3(-3.0, 3.0, 0.0);

	// Texture coordinates are passed through
	texcoord = texcoords;

	// Transform the vertex according to modelview
	vec4 viewVertex;	
	viewVertex = vertex * modelview;

	// Calculate lighting parameters for the fragment shader
	eye = normalize(viewVertex.xyz);
	light = normalize(lightPos - viewVertex.xyz);
	normal = normalize((vec4(vnormal, 0.0) * normalview).xyz);

	// Project and send to the fragment shader
	gl_Position = viewVertex * projection;
	gl_PointSize = 1.0f / (-viewVertex.z * 0.001f);
}
