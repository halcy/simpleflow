#version 420 core

// Composes different textures into one

// Parameters from the vertex shader
in vec2 coords;

// Textures
uniform sampler2D backgroundTexture;
uniform sampler2D particleTexture;

// Output
out vec4 outColor;

void main() {
	vec3 dir = normalize(vec3(coords.xy - vec2(0.5f), -2.0f));
	vec4 backgroundColor = texture(backgroundTexture, abs(dir.xy + vec2(0.5f)));
	vec4 particleColor = texture(particleTexture, coords);
	outColor = particleColor * particleColor.w + backgroundColor * (1.0f - particleColor.w);
}
