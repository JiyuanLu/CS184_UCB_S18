import React  from 'react';
import Section from './component/section';
import Model from './model.js';
import './App.css';

class App extends React.Component {
  render() {
    const sections = Model.map(m => {
      return <Section key={m.title} model={m}/>
    });
    return (
      <div className="App">
        <div className='header'>
          <h1>Project 3-2 Part 5 - Shaders</h1>
        </div>
        {sections}
      </div>
    );
  }
}

export default App;
