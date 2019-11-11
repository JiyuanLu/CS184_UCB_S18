attribute vec3 position;
attribute vec3 normal;		//n vector
attribute vec3 tangent;		//t vector
attribute vec2 uv;

uniform mat4 projectionMatrix;
uniform mat4 modelViewMatrix;
uniform mat4 modelMatrix;
uniform sampler2D textureDisplacement;
uniform vec2 textureDimension;

varying vec3 fPosition;
varying vec3 fNormal;

vec3 b;						//b vector		
mat3 TBN;
float du;
float dv;
float w;
float h;
float u;
float v;
float hu;
float hv;
float huv;
vec3 no;
vec3 nd;

// h(u, v) = texture2D(textureDisplacement, uv);

void main() {
    vec3 offset = position;

    // TODO: Compute displaced vertices
    float heightScaling = 0.8;		//kh

    // TODO: Compute displaced normals
    float normalScaling = 1.0;		//kn


    vec3 n = normal;
    vec3 t = tangent;
    vec3 b = cross(n, t);
    TBN = mat3(t, b, n);

    w = textureDimension.x;
    h = textureDimension.y;
    u = uv.x;
    v = uv.y;
    hu = texture2D(textureDisplacement, vec2(u + 1.0/w, v)).r;
    hv = texture2D(textureDisplacement, vec2(u, v + 1.0/h)).r;
    huv = texture2D(textureDisplacement, uv).r;
    du = heightScaling * normalScaling * (hu - huv);
    dv = heightScaling * normalScaling * (hv - huv);
    no = vec3(-du, -dv, 1.0);
    nd = TBN * no;
    fNormal = (modelMatrix * vec4(nd, 1.0)).xyz;

    fPosition = (modelMatrix * vec4(position, 1.0)).xyz;
    gl_Position = projectionMatrix * modelViewMatrix * vec4(offset, 1.0);
}