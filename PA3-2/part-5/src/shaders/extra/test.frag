precision highp float;

uniform float time;

varying vec2 fUv;

void main() {
    float t = time / 1000.;
    float r = fUv.y;
    float g = (cos(t) * 0.25 + 0.5);
    gl_FragColor = vec4(r, g, cos(t) * sin(t), 1.);
}