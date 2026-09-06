/* Display saved scores only. Selecting a point never changes its AAM mapping. */
function scorePlotGroups(assemblies) {
  const groups = new Map();
  assemblies.forEach((a, index) => {
    const s = a.score;
    const retention = s.set_atom_retention;
    const changes = s.broken_bonds + s.formed_bonds;
    if (!Number.isFinite(retention) || !Number.isFinite(changes)) return;
    const key = JSON.stringify([changes, retention]);
    if (!groups.has(key)) groups.set(key, {key, changes, retention, members: []});
    groups.get(key).members.push({a, index});
  });
  for (const group of groups.values()) {
    group.members.sort((x, y) =>
      x.a.score.matched_fragment_count - y.a.score.matched_fragment_count ||
      x.a.precursors.length - y.a.precursors.length || x.index - y.index);
  }
  return [...groups.values()];
}

function renderScorePlot(assemblies, selected, onSelect) {
  const svg = document.getElementById('scorePlot');
  const choices = document.getElementById('scorePlotChoices');
  svg.replaceChildren();
  choices.replaceChildren();
  const groups = scorePlotGroups(assemblies);
  const ns = 'http://www.w3.org/2000/svg';
  const element = (tag, attrs, text) => {
    const node = document.createElementNS(ns, tag);
    for (const [key, value] of Object.entries(attrs)) node.setAttribute(key, value);
    if (text !== undefined) node.textContent = text;
    svg.appendChild(node);
    return node;
  };
  const width = 720, height = 180, left = 60, right = 20, top = 12, bottom = 35;
  svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
  const xmax = Math.max(1, ...groups.map(g => g.changes)) + 1;
  const x = value => left + value / xmax * (width - left - right);
  const y = value => height - bottom - value * (height - top - bottom);
  for (let tick = 0; tick <= 4; tick++) {
    const value = tick / 4;
    element('line', {x1:left, x2:width-right, y1:y(value), y2:y(value), stroke:'#e2e8f0'});
    element('text', {x:left-8, y:y(value)+4, 'text-anchor':'end', 'font-size':11}, `${100*value}%`);
    const cost = Math.round(xmax * value);
    element('text', {x:x(cost), y:height-bottom+16, 'text-anchor':'middle', 'font-size':11}, cost);
  }
  element('text', {x:width/2, y:height-2, 'text-anchor':'middle', 'font-size':12}, 'Structural changes: breaking + forming → fewer is better');
  element('text', {x:13, y:height/2, transform:`rotate(-90 13 ${height/2})`, 'text-anchor':'middle', 'font-size':12}, 'R retention');
  for (const group of groups) {
    const active = group.members.some(m => m.index === selected);
    const blind = group.members.some(m => !m.a.ground_truth);
    const validation = group.members.some(m => m.a.ground_truth);
    const dot = element('circle', {cx:x(group.changes), cy:y(group.retention), r:active?9:7,
      fill:blind?'#2684ff':'#00a896', stroke:active?'#111827':validation?'#00a896':'white',
      'stroke-width':active?3:2, tabindex:0, role:'button', 'data-score-key':group.key});
    const description = `${(100*group.retention).toFixed(1)}% retention; ${group.changes} changes; ${group.members.length} alternatives. Click to open; fewer fragments first.`;
    dot.setAttribute('aria-label', description);
    const title = document.createElementNS(ns, 'title');
    title.textContent = description;
    dot.appendChild(title);
    dot.style.cursor = 'pointer';
    dot.onclick = () => onSelect(group.members[0].index);
    dot.onkeydown = event => {if(event.key === 'Enter' || event.key === ' ') {event.preventDefault();dot.onclick();}};
    if (group.members.length > 1) {
      element('text', {x:x(group.changes)+11, y:y(group.retention)-9, 'font-size':11, 'pointer-events':'none'}, `×${group.members.length}`);
    }
    if (active) {
      const caption = document.createElement('span');
      caption.textContent = `At this point (${group.members.length}): `;
      choices.appendChild(caption);
      for (const {a, index} of group.members) {
        const button = document.createElement('button');
        button.textContent = `${a.ground_truth?'Validation':'Alternative '+a.rank} · ${a.pattern} · ${a.score.matched_fragment_count} fragments`;
        button.className = index === selected ? 'chosen' : '';
        button.onclick = () => onSelect(index);
        choices.appendChild(button);
      }
    }
  }
}
