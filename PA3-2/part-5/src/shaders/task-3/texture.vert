attribute vec3 position;
attribute vec3 normal;
attribute vec2 uv;

uniform mat4 modelMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

varying vec2 fUv;

void main() {
  fUv = uv;
  gl_Position = projectionMatrix * modelViewMatrix * vec4(position, 1.0);
}