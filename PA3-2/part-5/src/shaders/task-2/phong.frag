precision highp float;

uniform vec3 cameraPosition;
uniform vec3 lPosition;
uniform vec3 lIntensity;

varying vec3 fPosition;
varying vec3 fNormal;

//parameters that you should decide
float ka = 0.5;
float kd = 0.8;
float ks = 0.8;
vec3 aIntensity = vec3(0.4, 0.4, 0.4);
float p = 128.0;

//intermediate results for ambience lighting
vec3 la;

//intermediate results for diffuse lighting
float r;
vec3 l;
vec3 fIntensity;
float cosine_term;
vec3 ld;

//intermediate results for specular lighting
vec3 v;
vec3 h;		//half vector between  l and v
vec3 ls;
float cosine_power;

vec3 L;

void main() {
    // TODO: Part 5.2

    //ambience
    la = ka * aIntensity;

    //diffuse
    r = distance(lPosition, fPosition);
    l = lPosition - fPosition;
    l = normalize(l);
    fIntensity = lIntensity / r / r;
    cosine_term = max(0.0,  dot(l, fNormal));
    ld = kd * fIntensity * cosine_term;

    //specular
    v = cameraPosition - fPosition;
    v = normalize(v);
    h = (l + v) / length(l + v);
    cosine_power = pow(max(0.0, dot(h, fNormal)), p);
    ls = ks * fIntensity * cosine_power;

    //total light
    L = la + ld + ls;
    //L = la;
    //L = ld;
    //L = ls;
    gl_FragColor = vec4(L, 1.0);
}