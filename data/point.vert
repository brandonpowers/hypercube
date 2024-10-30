uniform mat4 transform;
uniform float pointSize;
attribute vec4 vertex;
attribute float size;
attribute vec4 color;
varying vec4 vertColor;
void main() {
  vec4 pos = transform * vertex;
  gl_Position = pos;
  gl_PointSize = size * pointSize * (1.0 / -pos.z);
  vertColor = color;
}