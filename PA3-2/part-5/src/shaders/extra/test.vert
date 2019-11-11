#define M_PI 3.1415926535897932384626433832795

attribute vec3 position;
attribute vec3 normal;
attribute vec2 uv;

uniform mat4 modelMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;
uniform float time;

varying vec2 fUv;

void main() {
    vec3 offset = position;
    float t = time / 1000.;

    offset += normal * (0.5 * sin(t * 2. * M_PI) * sin(uv.x * 8. * M_PI) + 0.5);
    offset += normal * (0.5 * cos(t * uv.y * M_PI) + 0.5);

    fUv = uv;

    gl_Position = projectionMatrix * modelViewMatrix * vec4(offset, 1.0);
}