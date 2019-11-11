import React from 'react';
import RenderView from './render-view';
import './section.css';

function SubSection(props) {
  const {model} = props;
  const {title, renderer} = model;
  return (
    <div className='subsection'>
      <h3>{title}</h3>
      { renderer ? <RenderView renderer={renderer}/> : null }
    </div>
  );
}

function Section(props) {
  const {model} = props;
  const {title, renderer, subsections} = model;

  const subs = subsections && subsections.map(s => <SubSection key={s.title} model={s}/>);

  return (
    <div className='section'>
      <h2>{title}</h2>
      { renderer ? <RenderView renderer={renderer}/> : null }
      {subs ? <div className='subsections'> {subs} </div> : null}
    </div>
  )
}

export default Section;