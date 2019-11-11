attribute vec3 position;
attribute vec3 normal;
attribute vec3 tangent;
attribute vec2 uv;

uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 modelMatrix;
uniform sampler2D textureDisplacement;
uniform vec2 textureDimension;

varying vec3 fPosition;
varying vec3 fNormal;

void main() {
    vec3 offset = position;
    vec3 b = cross(normal, tangent);
    mat3 tbn = mat3(tangent, b, normal);

    float w = textureDimension[0];
    float h = textureDimension[1];
    vec2 t = vec2(uv.x + 1.0/w, uv.y);
    vec4 sap1 = texture2D(textureDisplacement, t) -texture2D(textureDisplacement, uv);
    vec4 sap2 = texture2D(textureDisplacement, t)-texture2D(textureDisplacement, uv);

    // TODO: Compute displaced vertices
    float heightScaling = 0.8;

    // TODO: Compute displaced normals
    float normalScaling = 1.0;

	float du = heightScaling * normalScaling * sap1.x;
    float dv = heightScaling * normalScaling * sap2.x;
    
    fPosition = vec3(0.0);
    fNormal = tbn * vec3(-du, -dv, 1);
    gl_Position = projectionMatrix * modelViewMatrix * vec4(offset, 1.0);
    
}