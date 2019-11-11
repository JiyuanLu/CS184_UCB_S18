import * as THREE from 'three';
import bindComputeTangents from './three-compute-tangents-poly';
import bindObjLoader from './three-obj-loader-poly';
import bindTeapotBufferGeometry from './three-teapot-buffer-geometry-poly';

bindComputeTangents(THREE);
bindObjLoader(THREE);
bindTeapotBufferGeometry(THREE);

export default THREE;