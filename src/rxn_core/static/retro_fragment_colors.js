/* Fragment membership and colors come from saved matched partitions. */
function fragmentMode() {
  return document.getElementById('colorBy').value === 'fragment' &&
    !document.getElementById('symmetry').checked;
}

function fragmentLegend(assembly) {
  const legend = document.getElementById('fragmentLegend');
  legend.replaceChildren();
  const note = document.createElement('span');
  note.textContent = fragmentMode()
    ? 'Matched fragments · same color in R and P · colors local to this assembly'
    : 'Precursor colors · symmetry domains are alternatives, not fragment assignments';
  legend.appendChild(note);
  if (!fragmentMode()) return;
  assembly.fragments.forEach(fragment => {
    const chip = document.createElement('span');
    chip.className = 'fragment-chip';
    const swatch = document.createElement('i');
    swatch.style.backgroundColor = fragment.color;
    chip.appendChild(swatch);
    chip.appendChild(document.createTextNode(fragment.label + ' · ' +
      fragment.source_atoms.length + ' atoms' +
      (fragment.occupations.length > 1 ? ' ×' + fragment.occupations.length : '')));
    chip.title = 'Source atoms: ' + fragment.source_atoms.join(', ') + '\n' +
      fragment.occupations.map(o => 'Copy ' + o.copy + ' → P atoms: ' +
        o.mapping.map(pair => pair[1]).join(', ')).join('\n');
    legend.appendChild(chip);
  });
}

function fragmentSupplierControls(assembly, wrap, render) {
  const heading = document.createElement('b');
  heading.textContent = 'Overlapping fragment colors · display alternatives';
  wrap.appendChild(heading);
  let count = 0;
  assembly.models.forEach((model, index) => {
    const target = index === assembly.precursors.length;
    const id = target ? 'P' : 'R' + index;
    model.fragmentAlternatives.forEach(group => {
      count++;
      const row = document.createElement('div');
      const label = document.createElement('div');
      label.textContent = (target ? 'P atoms: ' : 'Merged R' + (index + 1) + ' atoms: ') +
        group.atoms.join(', ');
      row.appendChild(label);
      group.owners.forEach(owner => {
        const fragment = assembly.fragments[owner];
        const button = document.createElement('button');
        button.textContent = fragment.label;
        button.style.color = fragment.color;
        button.style.border = group.selected === owner ? '2px solid currentColor' : '1px solid #ccc';
        button.title = 'Show this fragment color only; the saved mapping is unchanged';
        button.onclick = () => {group.selected = owner; render(id, model);};
        row.appendChild(button);
      });
      wrap.appendChild(row);
    });
  });
  if (!count) wrap.textContent = 'Each colored atom has one fragment color. Identical source fragments in repeated copies share a color.';
}
