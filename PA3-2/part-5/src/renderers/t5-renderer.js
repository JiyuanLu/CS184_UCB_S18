import THREE from '../lib/three';
import Renderer from '../lib/renderer';

import vert from '../shaders/task-5/texture.vert';
import frag from '../shaders/task-5/texture.frag';

import texture from '../textures/earth.png';

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

    const geometry = new THREE.TeapotBufferGeometry(4, 32, 32);
    const material = this.createShaderMaterial(vert, frag);
    const object = new THREE.Mesh(geometry, material);
    this.scene.add(object);
  }

  update(dt) {
    if (!this.focussed) {
      this.updateCamera(dt / 12000);
    }
  }
}
