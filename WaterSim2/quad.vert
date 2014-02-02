#version 420 core

// Screen-sized-quad vertex shader

in vec4 vertex;
out vec2 coords;

void main() {
	coords = (vertex.xy + 1.0f) / 2.0f;
	gl_Position = vertex;
}
