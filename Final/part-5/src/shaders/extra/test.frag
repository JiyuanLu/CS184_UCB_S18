precision highp float;

uniform float time;
uniform vec3 cameraPosition;
uniform vec3 lPosition;
uniform vec3 lIntensity;

varying vec3 fPosition;
varying vec3 fNormal;
varying vec2 fUv;
uniform vec2 textureDimension;
uniform sampler2D texture;

void main() {
	float t = time / 10000.;
    float ka = 0.5;
    float kd = 0.5;
    float ks = 0.1;
    float p = 5.0;
    float U = mod(1.0+fUv.x- t, 1.0);
    vec4 textcolor = texture2D(texture, vec2(U, fUv.y));
 //   vec3 Ia = vec3(0.1,0.1,0.1);
    vec3 Ia = vec3(textcolor.x, textcolor.y, textcolor.z);
    vec3 d = lPosition - fPosition;
    float r = 1.0/dot(d,d);
    vec3 l = normalize(d);
    vec3 v = cameraPosition - fPosition;
    v = normalize(d);
    vec3 h = l + v;
    h = normalize(h);
    float t1 = dot(fNormal, l);
    float t2 = dot(fNormal, h);
    
    vec3 color = ka * Ia;
    
    if (t1 > 0.0) 
    	color +=  kd*r*t1*lIntensity;
    if (t2 > 0.0)
    	color = color + ks * r * pow(t2,p) * lIntensity;
    gl_FragColor = vec4(color, 1.0);

    //gl_FragColor = texture2D(texture, fUv);
}