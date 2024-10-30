uniform mat4 transform;
uniform float pointSize;

attribute vec4 position;
attribute vec4 color;
attribute float size;

varying vec4 vertColor;

void main() {
    gl_Position = transform * position;
    gl_PointSize = size * pointSize * (1.0 / -gl_Position.z);
    vertColor = color;
}