#ifdef GL_ES
precision mediump float;
precision mediump int;
#endif
varying vec4 vertColor;
void main() {
  vec2 coord = gl_PointCoord - vec2(0.5);
  float r = length(coord) * 2.0;
  float a = 1.0 - smoothstep(0.0, 1.0, r);
  gl_FragColor = vec4(vertColor.rgb, vertColor.a * a);
}