import THREE from '../lib/three';
import Renderer from '../lib/renderer';

import vert from '../shaders/extra/test.vert';
import frag from '../shaders/extra/test.frag';

import texture from '../textures/global.png';

export default class extends Renderer {
  initScene() {
    if (!this.checkShader(vert, frag)) {
      this.setErrorScene();
      return;
    }
    this.setLight();

    this.uniforms['texture'] = {
      type: "t",
      value: new THREE.TextureLoader().load(texture)
    };
    this.uniforms['textureDimension'] = {
      type: 'vec2',
      value: new THREE.Vector2(100.0, 100.0)
    };    

    const geometry = new THREE.SphereGeometry(5, 32, 32);
    const material = this.createShaderMaterial(vert, frag);
    const cube = new THREE.Mesh(geometry, material);
    this.scene.add(cube);
  }

  update(dt) {
    if (!this.focussed) {
      this.updateCamera(dt / 12000);
    }
  }
}