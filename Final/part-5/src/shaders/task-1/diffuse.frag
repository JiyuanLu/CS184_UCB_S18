precision highp float;

uniform vec3 lPosition;
uniform vec3 lIntensity;

varying vec3 fPosition;
varying vec3 fNormal;

void main() {
    // TODO: Part 5.1
    
    vec3 d = lPosition - fPosition;
    float r = 1.0/dot(d,d);
    vec3 l = normalize(d);
    vec3 color = r * lIntensity;
    float nl = dot(l, fNormal);
    if (nl > 0.0)
	 	color = nl * color;
    else
    	color = vec3(0.0,0.0,0.0); 
    gl_FragColor = vec4(color, 1.0);

}