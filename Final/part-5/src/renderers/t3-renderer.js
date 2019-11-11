import THREE from '../lib/three';
import Renderer from '../lib/renderer';

import vert from '../shaders/task-3/texture.vert';
import frag from '../shaders/task-3/texture.frag';

import texture from '../textures/displacement2.png';

export default class extends Renderer {
  initScene() {
    if (!this.checkShader(vert, frag)) {
      this.setErrorScene();
      return;
    }

    this.uniforms['texture'] = {
      type: "t",
      value: new THREE.TextureLoader().load(texture)
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
