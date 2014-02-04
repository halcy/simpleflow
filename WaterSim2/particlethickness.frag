#version 420 core

// Particle thickness fragment shader

// Parameters from the vertex shader
in vec3 eyespacePos;
in float eyespaceRadius;

// Uniforms
uniform mat4 projection;
uniform vec2 screenSize;

// Textures
uniform sampler2D terrainTexture;

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
	// Blur-ish effect is achieved by including length here.
	particleThickness = 1.0f - length(normal);

	// Depth needs to be calculated so we can check against terrain
	// Set up rest of normal
	normal.z = sqrt(1.0f - dist);
	normal.y = -normal.y;
	normal = normalize(normal);

	// Calculate fragment position in eye space, project to find depth
	vec4 fragPos = vec4(eyespacePos + normal * eyespaceRadius / screenSize.y, 1.0);
	vec4 clipspacePos = fragPos * projection;

	// Set up output
	float far = gl_DepthRange.far; 
	float near = gl_DepthRange.near;
	float deviceDepth = clipspacePos.z / clipspacePos.w;
	float fragDepth = (((far - near) * deviceDepth) + near + far) / 2.0;
	gl_FragDepth = fragDepth;

	if(fragDepth > texture(terrainTexture, gl_FragCoord.xy / screenSize).w) {
		discard;
	}
}
