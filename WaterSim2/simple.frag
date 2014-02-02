#version 420 core

// Simple fragment shader.
// Does texturing and phong shading.

// Parameters from the vertex shader
in vec2 texcoord;
in vec3 eye;
in vec3 light;
in vec3 normal;

// Textures
uniform sampler2D texture;

// Output
out vec4 outColor;

void main() {
	vec3 lightAmbient = vec3(0.2, 0.2, 0.2);
	vec3 lightDiffuse = vec3(0.5, 0.5, 0.5);
	vec3 lightSpecular = vec3(1.0, 1.0, 1.0);

	vec3 materialColor = texture2D(texture, texcoord).rgb;
	vec3 materialSpecular = vec3(1.0, 1.0, 1.0);
	float materialShinyness = 10.0;

	// Normal
	vec3 normalNormal = normalize(normal);

	// Ambient lighting
	vec3 color = lightAmbient * materialColor;

	// Cosine of angle between normal and vector light-vertex
	float lambertTerm = dot(normalNormal, light);

	// Avoid darkening parts and drawing specular highlights wrong
	if(lambertTerm > 0.0) {
	
		// Diffuse lighting
		color += lightDiffuse * materialColor * lambertTerm;

		// Specular highlights
		vec3 specDir = normalize(reflect(light, -normalNormal));
		float specular = max(
			0.0,
			pow(dot(specDir, eye), materialShinyness)
		);
		color += lightSpecular *
			materialSpecular *
			specular;
	}

	outColor = vec4(color, 1.0f);
}
