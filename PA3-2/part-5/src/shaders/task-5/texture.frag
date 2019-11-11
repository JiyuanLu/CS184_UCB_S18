precision highp float;

uniform sampler2D texture;
varying vec2 fUv;

void main() {
	 gl_FragColor = texture2D(texture, fUv);
}