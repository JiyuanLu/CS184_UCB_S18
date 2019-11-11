import * as Renderers from './renderers';

export default [
  {
    title: '1. Diffuse Shading',
    renderer: Renderers.T1Renderer
  },
  {
    title: '2. Blinn-Phong Shading',
    renderer: Renderers.T2Renderer
  },
  {
    title: '3. Texture Mapping',
    renderer: Renderers.T3Renderer
  },
  {
    title: '4. Bump/Displacement Mapping',
    subsections: [
      {
        title: '4.1 Bump Mapping',
        renderer: Renderers.T4_1Renderer
      },
      {
        title: '4.2 Displacement Mapping',
        renderer: Renderers.T4_2Renderer
      }
    ]
  },
  {
    title: 'Add any extra shaders here!',
    renderer: Renderers.T5Renderer
  }
]