#version 420 core

// Shades liquids

// Parameters from the vertex shader
in vec2 coords;

// Textures
uniform sampler2D environmentTexture;
uniform sampler2D particleTexture;
uniform sampler2D particleThicknessTexture;
uniform sampler2D velocityTexture;

// Uniforms
uniform vec2 screenSize;
uniform mat4 projection;
uniform mat4 modelview;

// Output
out vec4 outColor;

float noise(vec2 coord) {
    return fract(sin(dot(coord.xy, vec2(12.9898, 78.233))) * 43758.5453);
}

vec2 spheremap(vec3 dir) {
	float m = 2.0f * sqrt(dir.x * dir.x + dir.y * dir.y + (dir.z + 1.0f) * (dir.z + 1.0f));
	return vec2(dir.x / m + 0.5f, dir.y / m + 0.5f);
}

vec3 eyespacePos(vec2 pos) {
	float depth = texture(particleTexture, pos);
	pos = (pos - vec2(0.5f)) * 2.0f;
	return(depth * vec3(-pos.x * projection[0][0], -pos.y * projection[1][1], 1.0f));
}

// Compute eye-space normal. Adapted from PySPH.
vec3 eyespaceNormal(vec2 pos) {
	// Width of one pixel
	vec2 dx = vec2(1.0f / screenSize.x, 0.0f);
	vec2 dy = vec2(0.0f, 1.0f / screenSize.y);

	// Central z
	float zc =  texture(particleTexture, pos);

	// Derivatives of z
	// For shading, one-sided only-the-one-that-works version
	float zdxp = texture(particleTexture, pos + dx);
	float zdxn = texture(particleTexture, pos - dx);
	float zdx = (zdxp == 0.0f) ? (zdxn == 0.0f ? 0.0f : (zc - zdxn)) : (zdxp - zc);

	float zdyp = texture(particleTexture, pos + dy);
	float zdyn = texture(particleTexture, pos - dy);
	float zdy = (zdyp == 0.0f) ? (zdyn == 0.0f ? 0.0f : (zc - zdyn)) : (zdyp - zc);

	// Projection inversion
	float cx = 2.0f / (screenSize.x * -projection[0][0]);
	float cy = 2.0f / (screenSize.y * -projection[1][1]);

	// Screenspace coordinates
	float sx = floor(pos.x * (screenSize.x - 1.0f));
	float sy = floor(pos.y * (screenSize.y - 1.0f));
	float wx = (screenSize.x - 2.0f * sx) / (screenSize.x * projection[0][0]);
	float wy = (screenSize.y - 2.0f * sy) / (screenSize.y * projection[1][1]);

	// Eyespace position derivatives
	vec3 pdx = normalize(vec3(cx * zc + wx * zdx, wy * zdx, zdx));
	vec3 pdy = normalize(vec3(wx * zdy, cy * zc + wy * zdy, zdy));

	return normalize(cross(pdx, pdy));
}

// Calculate fresnel coefficient
// Schlicks approximation is for lamers
float fresnel(float rr1, float rr2, vec3 n, vec3 d) {
	float r = rr1 / rr2;
	float theta1 = dot(n, -d);
	float theta2 = sqrt(1.0f - r * r * (1.0f - theta1 * theta1));

	// Figure out what the Fresnel equations say about what happens next
	float rs = (rr1 * theta1 - rr2 * theta2) / (rr1 * theta1 + rr2 * theta2);
	rs = rs * rs;
	float rp = (rr1 * theta2 - rr2 * theta1) / (rr1 * theta2 + rr2 * theta1);
	rp = rp * rp;

	return((rs + rp) / 2.0f);
}

void main() {
	float particleDepth = texture(particleTexture, coords);
	float particleThickness = texture(particleThicknessTexture, coords);
	float velocity = texture(velocityTexture, coords);

	vec3 normal = eyespaceNormal(coords);
	normal = normal * inverse(mat3(modelview));
	normal.xz = normal.zx;

	vec3 lightDir = vec3(1.0f, 1.0f, -1.0f);

	if(particleDepth == 0.0f) {
		outColor = vec4(0.0f);
	}
	else {
		vec3 pos = eyespacePos(coords);
		pos = (vec4(pos, 1.0f) * inverse(modelview)).xyz;

		float thickness = vec4(particleThickness) / 10.0f;

		float lambert = max(0.0f, dot(normalize(lightDir), normal));
		
		vec3 fromEye = normalize(pos);
		
		fromEye.xz = -fromEye.xz;
		vec3 reflectedEye = normalize(reflect(fromEye, normal));
		float specular = clamp(fresnel(1.0f, 1.5f, normal, fromEye), 0.0f, 0.4f);

		// De-specularize fast particles
		//specular = max(0.0f, specular - (velocity / 15.0f));
		
		vec4 environmentColor = texture(environmentTexture, spheremap(reflectedEye));
		vec4 particleColor = exp(-vec4(0.6f, 0.2f, 0.05f, 3.0f) * thickness);
		particleColor.w = clamp(1.0f - particleColor.w, 0.0f, 1.0f);
		particleColor.rgb = (lambert + 0.2f) * particleColor.rgb * (1.0f - specular) + specular * environmentColor.rgb;
		
		// Oil
		/*particleColor.rgb = specular * environmentColor.rgb;
		particleColor.w = 1.0f;*/

		// Add some superfake foam colouring
		//particleColor += lambert * vec4(velocity / 15.0f);
		
		outColor = particleColor;
		//outColor = environmentColor;
		//outColor = vec4(lambert, lambert, lambert, 1.0f);
	}
}
