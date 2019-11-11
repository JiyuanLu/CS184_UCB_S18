attribute vec3 position;
attribute vec3 normal;

uniform mat4 modelMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

varying vec3 fPosition;
varying vec3 fNormal;

void main() {
	// TODO: Part 5.1
	vec4 fp = modelMatrix * vec4(position, 1.0);
	vec4 fn = modelMatrix * vec4(normal, 0.0);

    fPosition = vec3(fp.x, fp.y, fp.z);
    fNormal = vec3(fn.x, fn.y, fn.z);
    gl_Position = projectionMatrix*(modelViewMatrix*fp);
}