attribute vec3 position;
attribute vec3 normal;

uniform mat4 modelMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 projectionMatrix;

varying vec3 fPosition;
varying vec3 fNormal;

vec3 vPosition;

void main() {
	// TODO: Part 5.1
    //fPosition = vec3(0.0);
    fPosition = (modelMatrix * vec4(position, 1.0)).xyz;			//world space position
    //fNormal = vec3(0.0);
    fNormal = (modelMatrix * vec4(normal, 1.0)).xyz;				//world space normal vector
    vPosition = (modelViewMatrix * vec4(position, 1.0)).xyz;
    gl_Position = (projectionMatrix * vec4(vPosition, 1.0));	//screen space position
    //gl_Position = vec4(vec3(0.0), 1.0);
}