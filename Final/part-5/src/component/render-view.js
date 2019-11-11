import React from 'react';
import * as Stats from 'stats-js';
import './render-view.css';

class RenderView extends React.Component {
  constructor(props) {
    super(props);
    this.state = {
      renderer: new props.renderer()
    }
  }

  render() {
    const {width, height} = this.props;
    return (
      <div className='container'>
        <div className='render-view-container' ref={$ => this.container = $}>
          <canvas
            className='render-view'
            width={width || 680}
            height={height || 400}
            ref={$ => this.$ = $}
          />
        </div>
      </div>
    );
  }

  componentDidMount() {
    const {renderer} = this.state;
    renderer.setView(this.$);

    const stats = new Stats();
    this.container.appendChild(stats.domElement);
    renderer.setStats(stats);

    renderer.init();
    renderer.begin();
  }

  componentWillUnmount() {
    const {renderer} = this.state;
    renderer.destroy();
  }
}

export default RenderView;