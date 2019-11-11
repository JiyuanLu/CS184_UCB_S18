precision highp float;

uniform vec3 cameraPosition;
uniform vec3 lPosition;
uniform vec3 lIntensity;

varying vec3 fPosition;
varying vec3 fNormal;

void main() {
    // TODO: Part 5.2
    
    float ka = 0.5;
    float kd = 0.5;
    float ks = 0.5;
    float p = 15.0;
    vec3 Ia = vec3(0.3,0.3,0.3);
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
}