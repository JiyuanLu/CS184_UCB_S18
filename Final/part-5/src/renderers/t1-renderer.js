import THREE from '../lib/three';
import Renderer from '../lib/renderer';

// shader imports
import vert from '../shaders/task-1/diffuse.vert';
import frag from '../shaders/task-1/diffuse.frag';

export default class extends Renderer {
  initScene() {
    if (!this.checkShader(vert, frag)) {
      this.setErrorScene();
      return;
    }

    this.setLight();

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
