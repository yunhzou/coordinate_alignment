"""Browser checks for shared molecular/catalog styles and saved-map identity."""
import json
import os
from pathlib import Path

from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'docs/viewer-validation'


def check_endpoint_overlays(frame):
    return frame.evaluate('''() => {
      const mech=findMech(currentMechId),expected={R:[],P:[]},seen={R:[],P:[]};
      const key=(color,pair)=>color+':'+[...pair].sort((a,b)=>a-b).join(',');
      for(const event of mech.event_records) {
        const side=event.wbo[1]<event.wbo[0]?'R':'P',color=side==='R'?'red':'green';
        expected[side].push(key(color,side==='R'||rOrdered?event.r:event.p));
      }
      for(const call of window.bondOverlayCalls) seen[call.side].push(key(call.color,call.pair));
      for(const side of ['R','P']) {
        if(JSON.stringify(seen[side].sort())!==JSON.stringify(expected[side].sort()))
          throw Error('Incorrect or duplicated '+side+' bond overlays: '+JSON.stringify(seen));
      }
      return {case:DATA.index,mechanism:mech.id,product_frame:rOrdered?'aligned':'native',
              R_red:seen.R.length,P_green:seen.P.length,total:seen.R.length+seen.P.length,status:'passed'};
    }''')


def main():
    OUT.mkdir(exist_ok=True)
    environment = dict(os.environ)
    libraries = ROOT / 'manuscript/build/browser-libs/root/usr/lib/x86_64-linux-gnu'
    if libraries.exists():
        environment['LD_LIBRARY_PATH'] = str(libraries) + ':' + environment.get('LD_LIBRARY_PATH', '')
    cached = Path('/h/399/yunhengzou/.cache/ms-playwright/chromium-1243/chrome-linux64/chrome')
    executable = os.environ.get('MANUSCRIPT_BROWSER', str(cached) if cached.exists() else None)
    errors, external, checks, color_checks = [], [], [], []
    with sync_playwright() as pw:
        browser = pw.chromium.launch(executable_path=executable, env=environment, headless=True,
                                     args=['--no-sandbox', '--enable-unsafe-swiftshader'])
        page = browser.new_page(viewport={'width':1440, 'height':950})
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.on('request', lambda r: external.append(r.url) if r.url.startswith(('http:', 'https:')) else None)
        folder = ROOT / 'reports/holdout_missing_pattern_seeds_20260910'
        page.goto((folder / 'viewer.html').as_uri())
        for option, index, count in ((0, 11, 2), (1, 64, 3), (2, 101, 4)):
            page.locator('#case').select_option(str(option))
            frame = page.frames[1]
            frame.wait_for_function(f'typeof DATA !== "undefined" && DATA.index === {index}')
            frame.wait_for_selector('#vw_R canvas')
            frame.evaluate('''() => {
              window.bondOverlayCalls=[];
              const draw=drawBonds,paint=render;
              drawBonds=(v,xyz,pairs,color)=>{
                const m=findMech(currentMechId),side=xyz===DATA.reactant.coords?'R':
                  [DATA.product.coords,m.product_xyz_in_R,m.product_xyz_in_R_aligned].includes(xyz)?'P':null;
                if(side) for(const pair of pairs) window.bondOverlayCalls.push({side,pair:[...pair],color});
                return draw(v,xyz,pairs,color);
              };
              render=()=>{window.bondOverlayCalls=[];return paint()};
            }''')
            assert frame.locator('#mech-sel button').count() == count
            assert '__TITLE__' not in frame.locator('h2').inner_text()
            for mechanism in range(1, count + 1):
                frame.locator(f'#mech-sel button[data-id="{mechanism}"]').click()
                assert frame.evaluate('currentMechId') == mechanism
                assert frame.evaluate("findMech(currentMechId).event_records.length") == (5 if index == 64 else 4)
                color_checks.append(check_endpoint_overlays(frame))
                # Display fitting may rotate/translate P; the selected correspondence stays fixed.
                assert frame.evaluate('''() => {
                  const m=findMech(currentMechId),r=molecularViewers.vw_R.selectedAtoms({}),p=molecularViewers.vw_P.selectedAtoms({});
                  return r.length===DATA.n_atoms && p.length===DATA.n_atoms && r.every((a,i)=>
                    [a.x,a.y,a.z].every((x,j)=>Math.abs(x-DATA.reactant.coords[i][j])<1e-9)) && p.every((a,i)=>
                    [a.x,a.y,a.z].every((x,j)=>Math.abs(x-m.product_xyz_in_R_aligned[i][j])<1e-9));
                }''')
                frame.locator('#rOrdered').uncheck()
                color_checks.append(check_endpoint_overlays(frame))
                assert frame.evaluate('''() => molecularViewers.vw_P.selectedAtoms({}).every((a,i)=>
                    [a.x,a.y,a.z].every((x,j)=>Math.abs(x-DATA.product.coords[i][j])<1e-9))''')
                frame.locator('#rOrdered').check()
                checks.append({'case':index, 'mechanism':mechanism, 'native_and_aligned_geometry':'passed'})
            frame.locator('#mech-sel button[data-id="1"]').click()
            frame.locator('#showAtomIndices').check()
            frame.locator('#showAtomIndices').uncheck()
            frame.locator('#event-details').evaluate('(node)=>node.open=true')
            page.screenshot(path=str(OUT / f'reaction-{index}.png'), full_page=True)
            print(f'Checked original-style case {index}: {count} saved witnesses.', flush=True)
        page.set_viewport_size({'width':390, 'height':844})
        assert page.evaluate('document.documentElement.scrollWidth') <= 390
        assert page.frames[1].evaluate('document.documentElement.scrollWidth') <= 390
        page.set_viewport_size({'width':1440, 'height':950})
        for relative, layout in [
            ('reports/holdout_forward_cap1000_seed1_20260910/viewer.html', 'minimum'),
            ('reports/golden_remaining_21/viewer.html', 'golden_remaining'),
            ('reports/golden_alternatives_20260907/case9/viewer.html', 'golden_mapping'),
            ('reports/eight_case_studies/t02_chloro_quinoline_alcohol/viewer.html', 'catalog'),
            ('manuscript/animations/index.html', 'growth_animation'),
        ]:
            page.goto((ROOT / relative).as_uri())
            assert page.locator('body').get_attribute('data-viewer-layout') == layout
            if layout == 'minimum':
                page.locator('#case').select_option('64')
                assert page.locator('#cards .card').count() >= 3
            elif layout == 'golden_remaining':
                page.locator('#candidate').select_option(index=1)
                assert page.locator('.atom-hit').count() > 0
                page.screenshot(path=str(ROOT / 'reports/golden_remaining_21/browser_check.png'), full_page=True)
            elif layout == 'golden_mapping':
                page.locator('#mode').select_option(index=1)
                assert page.locator('.atom-hit').count() > 0
            elif layout == 'catalog':
                assert page.locator('canvas').count() >= 2
                page.screenshot(path=str(OUT / 'catalog.png'))
            else:
                page.locator('#next').click()
                assert page.evaluate('animationAPI.getState().step') == 1
            checks.append({'page':relative, 'layout':layout, 'status':'passed'})
        assert not errors, errors
        assert not external, external
        browser.close()
    (OUT / 'checks.json').write_text(json.dumps(dict(status='passed', checks=checks,
        javascript_errors=errors, external_requests=external,
        scope='Presentation and actual saved mapping checks; no benchmark search.'), indent=2) + '\n')
    (OUT / 'red-green-checks.json').write_text(json.dumps(dict(status='passed',
        convention='R: red broken/weakened; P: green formed/strengthened; each event shown once.',
        checks=color_checks, javascript_errors=errors), indent=2) + '\n')
    print('Shared viewer browser checks passed.')


if __name__ == '__main__':
    main()
