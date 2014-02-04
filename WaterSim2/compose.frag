#version 420 core

// Composes different textures into one

// Parameters from the vertex shader
in vec2 coords;

// Textures
uniform sampler2D backgroundTexture;
uniform sampler2D terrainTexture;
uniform sampler2D particleTexture;

// Uniforms
uniform mat4 modelview;

// Output
out vec4 outColor;

vec2 spheremap(vec3 dir) {
	float m = 2.0f * sqrt(dir.x * dir.x + dir.y * dir.y + (dir.z + 1.0f) * (dir.z + 1.0f));
	return vec2(dir.x / m + 0.5f, dir.y / m + 0.5f);
}

void main() {
	vec3 dir = normalize(vec3((coords.xy - vec2(0.5f)) * 2.0f, -1.0f));
	dir = mat3(modelview) * dir;
	vec2 mapCoords = spheremap(dir);

	vec4 terrainColor = texture(terrainTexture, coords);
	vec4 backgroundColor = texture(backgroundTexture, mapCoords);
	vec4 particleColor = texture(particleTexture, coords);

	float blend = terrainColor.a == 1.0f ? 0.0f : 1.0f;
	backgroundColor = blend * terrainColor + (1.0f - blend) * backgroundColor;
	outColor = particleColor * particleColor.w + backgroundColor * (1.0f - particleColor.w);
}
