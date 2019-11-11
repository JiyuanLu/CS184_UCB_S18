precision highp float;

uniform vec3 lPosition;
uniform vec3 lIntensity;

varying vec3 fPosition;
varying vec3 fNormal;

float r;
float kd = 1.0;
vec3 l;
vec3 ld;
vec3 fIntensity;
float cosine_term;

void main() {
    // TODO: Part 5.1
    //gl_FragColor = vec4(vec3(0.0), 1.0);
    r = distance(lPosition, fPosition);
    l = lPosition - fPosition;
    l = normalize(l);
    fIntensity = lIntensity / r / r;
    cosine_term = max(0.0,  dot(l, fNormal));
    ld = kd * fIntensity * cosine_term;
    gl_FragColor = vec4(ld, 1.0);

}