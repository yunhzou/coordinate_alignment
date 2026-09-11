"""Browser validation for any generated 3D search-trajectory bundle."""
import argparse
import json
import os
from pathlib import Path
import numpy as np
from playwright.sync_api import sync_playwright

ROOT=Path(__file__).resolve().parents[1]


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('bundle',type=Path)
    a=p.parse_args();bundle=a.bundle.resolve();data=json.loads((bundle/'trace.json').read_text())
    environment=dict(os.environ)
    libraries=ROOT/'manuscript/build/browser-libs/root/usr/lib/x86_64-linux-gnu'
    if libraries.exists():environment['LD_LIBRARY_PATH']=str(libraries)+':'+environment.get('LD_LIBRARY_PATH','')
    cached=Path('/h/399/yunhengzou/.cache/ms-playwright/chromium-1243/chrome-linux64/chrome')
    executable=os.environ.get('MANUSCRIPT_BROWSER',str(cached) if cached.exists() else None)
    errors,external,checks=[],[],[]
    with sync_playwright() as pw:
        browser=pw.chromium.launch(executable_path=executable,env=environment,args=['--no-sandbox','--enable-unsafe-swiftshader'])
        page=browser.new_page(viewport={'width':1280,'height':1050})
        page.on('pageerror',lambda e:errors.append(str(e)))
        page.on('request',lambda r:external.append(r.url) if r.url.startswith(('http:','https:')) else None)
        page.goto((bundle/'algorithm_trajectory.html').as_uri())
        page.wait_for_function('typeof searchTrajectory!=="undefined"')
        for ri,run in enumerate(data['runs']):
            page.evaluate('(i)=>searchTrajectory.setRun(i)',ri)
            assert page.locator('#path option').count()==len(run['paths'])
            for pi,path in enumerate(run['paths']):
                page.evaluate('(i)=>searchTrajectory.setPath(i)',pi)
                start=page.evaluate('searchTrajectory.getState()')['frame']
                frame=path['frames'][start]
                if run.get('focus'):
                    assert frame['kind']=='pop' and frame['trial']['edge']==run['focus']['source_edge']
                    assert page.locator('.selected-trial').count()==1
                    assert f'P{run["focus"]["target"]}' in page.locator('.selected-trial').inner_text()
                    assert len(frame['trial']['trials'][frame['preferred']])==page.locator('#trials tr').count()-1
                    compatible=next((r for r in frame['trial']['trials'][frame['preferred']] if r['status']=='compatible'),None)
                    if compatible:
                        page.get_by_role('button',name=f'P{compatible["target"]}',exact=True).click()
                        expected_preview={**frame['locked'],**compatible['witness']}
                        assert page.evaluate('searchTrajectory.mapping()')==expected_preview
                branching=next((i for i,f in enumerate(path['frames']) if len(f['candidates'])>1),None)
                if branching is not None:
                    page.evaluate('(i)=>searchTrajectory.setFrame(i)',branching)
                    page.locator('#candidates button').nth(1).click()
                    assert page.evaluate('searchTrajectory.getState().candidate')==1
                    expected_preview={**path['frames'][branching]['locked'],**path['frames'][branching]['candidates'][1]['witness']}
                    assert page.evaluate('searchTrajectory.mapping()')==expected_preview
                for fi in sorted({0,start,min(start+1,len(path['frames'])-1),len(path['frames'])-1}):
                    page.evaluate('(i)=>searchTrajectory.setFrame(i)',fi)
                    for side,key in [('r','reactant'),('p','product')]:
                        coords=page.evaluate('(side)=>viewers[side].selectedAtoms({}).map(a=>[a.x,a.y,a.z])',side)
                        np.testing.assert_allclose(coords,data['input'][key]['coordinates'],atol=1e-12)
                    assert page.locator('#candidates button').count()==len(path['frames'][fi]['candidates'])
                # Final overlays represent each actual saved event once, R losses / P gains.
                page.evaluate('''() => {
                    const original=bond;window.eventDraws=[];
                    bond=(side,pair,color,...args)=>{if(['red','green'].includes(color))window.eventDraws.push({side,pair,color});return original(side,pair,color,...args)};
                    paint();bond=original;
                }''')
                expected=sorted(('r' if e['kind'] in ('broken','weakened') else 'p',
                    'red' if e['kind'] in ('broken','weakened') else 'green',
                    tuple(sorted(e['r'] if e['kind'] in ('broken','weakened') else e['p']))) for e in path['events'])
                actual=sorted((e['side'],e['color'],tuple(sorted(e['pair']))) for e in page.evaluate('eventDraws'))
                assert actual==expected
                page.evaluate('searchTrajectory.setFrame(0)');page.locator('#next').click()
                assert page.evaluate('searchTrajectory.getState().frame')==1
                page.locator('#prev').click();assert page.evaluate('searchTrajectory.getState().frame')==0
                page.locator('#play').click();page.wait_for_function('searchTrajectory.getState().frame>=2')
                page.locator('#play').click();assert not page.evaluate('searchTrajectory.getState().playing')
                page.locator('#critical').click() if run.get('focus') else page.evaluate('searchTrajectory.setFrame(0)')
                if path['terminal']==run['default_terminal']:
                    page.screenshot(path=str(bundle/f'run_{ri}_decision.png'),full_page=True)
                checks.append(dict(context=run['context'],terminal=path['terminal'],fixed_geometry='passed',
                                   event_colors='passed',candidate_selection='passed',playback='passed'))
                print(f'Checked context {run["context"]}, terminal {path["terminal"]}.',flush=True)
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth')<=390
        assert not errors,errors
        assert not external,external
        browser.close()
    (bundle/'browser_validation.json').write_text(json.dumps(dict(status='passed',checks=checks,
        javascript_errors=errors,external_requests=external),indent=2)+'\n')


if __name__=='__main__':main()
