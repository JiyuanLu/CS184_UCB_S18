import THREE from '../lib/three';

class Renderer {
  destroy() {
    this.detachListeners();
    this.view = null;
    this.scene = null;
    this.camera = null;
    this.renderer = null;
  }

  setView(view) {
    this.view = view;
    this.width = view.width;
    this.height = view.height;
  }

  setStats(stats) {
    this.stats = stats;
  }

  attachListeners() {
    window.addEventListener   ('resize',      this.onresize);
    this.view.addEventListener('keydown',     this.onkeydown);
    this.view.addEventListener('keyup',       this.onkeyup);
    window.addEventListener('mousedown',   this.onmousedown);
    window.addEventListener('mouseup',     this.onmouseup);
    window.addEventListener('mousemove',   this.onmousemove);
  }

  detachListeners() {
    window.removeEventListener('resize',      this.onresize);
    this.view.addEventListener('keydown',     this.onkeydown);
    this.view.addEventListener('keyup',       this.onkeyup);
    window.addEventListener('mousedown',   this.onmousedown);
    window.addEventListener('mouseup',     this.onmouseup);
    window.addEventListener('mousemove',   this.onmousemove);
  }

  onresize = () => {
    this.width = this.view.width;
    this.height = this.view.height;
    this.renderer.setViewport(0, 0, this.width, this.height);
  };

  onkeydown = (e) => { };
  onkeydown = (e) => { };
  onkeyup = (e) => { };

  onmousedown = (e) => {
    if (!this.view.contains(e.target)) {
      this.focussed = false;
      return;
    }
    this.focussed = true;
    this.inputState.mouse[e.button] = true;
    if (e.button === 0)
      this.inputState.lastMouse = [e.clientX, e.clientY];
  };

  onmouseup = (e) => {
    this.inputState.mouse[e.button] = false;
    if (e.button === 0)
      this.inputState.lastMouse = null;
  };

  onmousemove = (e) => {
    if (!this.inputState.mouse[0]) {
      return;
    }

    const x = e.clientX;
    const y = e.clientY;

    const dx = (x - this.inputState.lastMouse[0]) / this.view.width;
    const dy = (this.inputState.lastMouse[1] - y) / this.view.height;
    this.inputState.lastMouse = [x, y];
    this.updateCamera(dx, dy);
  };

  init() {
    this.scene = new THREE.Scene();
    this.camera = new THREE.PerspectiveCamera(75, this.width / this.height, 0.1, 1000);
    this.cameraState = {
      theta: 1.8,
      phi: 0,
      distance: 15,
    };
    this.updateCamera();

    this.renderer = new THREE.WebGLRenderer({
      canvas: this.view,
      antialias: true
    });
    this.renderer.setClearColor(new THREE.Color(0.078, 0.09, 0.11));
    this.uniforms = {
      time: { type: 'f', value: 0 }
    };
    this.inputState = {
      key: {},
      mouse: {},
      lastMouse: null
    };
    this.focussed = false;
    this.stats.setMode(0);
    this.initScene();
    this.attachListeners();
  }

  setErrorScene() {
    const sphereBuffer = new THREE.SphereBufferGeometry(5, 32, 32);
    const wireframe = new THREE.WireframeGeometry(sphereBuffer);
    const object = new THREE.LineSegments(wireframe);
    object.material.depthTest = false;
    this.scene.add(object);
  }

  setLight(x=25, y=15, z=25, r=1, g=1, b=1) {
    const lGeometry = new THREE.BoxGeometry(1, 1, 1);
    const lMaterial = new THREE.MeshBasicMaterial({ color: new THREE.Color(r, g, b) });
    const light = new THREE.Mesh(lGeometry, lMaterial);
    light.position.set(x, y, z);
    this.scene.add(light);

    this.uniforms['lPosition'] = {
      t: 'vec3',
      value: light.position
    };

    this.uniforms['lIntensity'] = {
      t: 'vec3',
      value: new THREE.Vector3(r * 1000, g * 1000, b * 1000)
    }
  }

  begin() {
    this.sysLast = 0;
    this.onFrame(this.sysLast);
  }

  onFrame(t) {
    const dt = t - this.sysLast;
    this.stats.begin();
    if (dt >= 1000 / 30) {
      this.uniforms.time.value = t;
      this.update(dt);
      this.render();
      this.sysLast = t;
      this.stats.end();
    }
    requestAnimationFrame(t => this.onFrame(t));
  }

  createShaderMaterial(vertexShader, fragmentShader) {
    return new THREE.RawShaderMaterial( {
      uniforms: this.uniforms,
      vertexShader,
      fragmentShader
    });
  }

  updateCamera(dx=0, dy=0) {
    let { theta, phi, distance } = this.cameraState;

    theta = theta - dy * Math.PI;
    phi = (phi + dx * 2 * Math.PI) % (2 * Math.PI);

    if (theta > 3.14) {
      theta = 3.14;
    } else if (theta < .001) {
      theta = 0.001;
    }
    this.cameraState.theta = theta;
    this.cameraState.phi = phi;

    this.camera.position.y = distance * -Math.cos(theta);
    this.camera.position.x = distance * Math.cos(phi) * Math.sin(theta);
    this.camera.position.z = distance * Math.sin(phi) * Math.sin(theta);
    this.camera.lookAt(new THREE.Vector3(0, 0, 0));
  }

  render() {
    this.renderer.render(this.scene, this.camera);
  }

  checkShader(vertSrc, fragSrc) {
    const gl = this.renderer.getContext();
    const vert = Renderer.compileShader(gl, vertSrc, gl.VERTEX_SHADER);
    if (vert === null) {
      return false;
    }
    const frag = Renderer.compileShader(gl, fragSrc, gl.FRAGMENT_SHADER);
    if (frag === null) {
      return false;
    }
    const prog = Renderer.createProgram(gl, vert, frag);
    if (prog === null) {
      return false;
    }
    gl.deleteShader(vert);
    gl.deleteShader(frag);
    gl.deleteProgram(prog);
    return true;
  }

  static compileShader(gl, shaderSource, shaderType) {
    const shader = gl.createShader(shaderType);
    gl.shaderSource(shader, shaderSource);
    gl.compileShader(shader);
    const success = gl.getShaderParameter(shader, gl.COMPILE_STATUS);
    if (!success) {
      console.log("could not compile shader:" + gl.getShaderInfoLog(shader));
      gl.deleteShader(shader);
      return null;
    }
    return shader;
  }

  static createProgram(gl, vertexShader, fragmentShader) {
    const program = gl.createProgram();
    gl.attachShader(program, vertexShader);
    gl.attachShader(program, fragmentShader);
    gl.linkProgram(program);
    const success = gl.getProgramParameter(program, gl.LINK_STATUS);
    if (!success) {
      console.log("program filed to link:" + gl.getProgramInfoLog(program));
      gl.deleteProgram(program);
      return null;
    }
    return program;
  };

  initScene() { }
  update(dt) { }
}

export default Renderer;