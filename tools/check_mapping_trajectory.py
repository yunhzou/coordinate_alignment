"""Check the original 3D trajectory against every saved diagnostic witness."""
import json
import os
from pathlib import Path
import tempfile

import numpy as np
from playwright.sync_api import sync_playwright

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / 'reports/holdout_missing_pattern_seeds_20260910'
OUT = ROOT / 'docs/viewer-validation'


def main():
    source = json.loads((REPORT / 'viewer_data.json').read_text())
    environment = dict(os.environ)
    libraries = ROOT / 'manuscript/build/browser-libs/root/usr/lib/x86_64-linux-gnu'
    if libraries.exists():
        environment['LD_LIBRARY_PATH'] = str(libraries) + ':' + environment.get('LD_LIBRARY_PATH', '')
    executable = os.environ.get('MANUSCRIPT_BROWSER',
        '/h/399/yunhengzou/.cache/ms-playwright/chromium-1243/chrome-linux64/chrome')
    errors, external, checks = [], [], []
    with sync_playwright() as pw, tempfile.TemporaryDirectory() as downloads:
        browser = pw.chromium.launch(executable_path=executable, env=environment,
            args=['--no-sandbox', '--enable-unsafe-swiftshader'])
        page = browser.new_page(viewport={'width':1440, 'height':1100}, accept_downloads=True)
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.on('request', lambda r: external.append(r.url) if r.url.startswith(('http:', 'https:')) else None)
        page.goto((REPORT / 'viewer.html').as_uri() + '#case=101')
        assert page.locator('#case').input_value() == '2'
        for ci, case in enumerate(source):
            page.locator('#case').select_option(str(ci))
            frame = page.frames[1]
            frame.wait_for_function(f'typeof DATA !== "undefined" && DATA.index === {case["index"]} && interpViewer !== null')
            frame.evaluate('''() => {
                const draw=drawBonds,paint=renderInterpolationFrame;
                window.trajectoryOverlays=[];
                drawBonds=(v,xyz,pairs,color)=>{
                    if(v===interpViewer) for(const pair of pairs) window.trajectoryOverlays.push({pair:[...pair],color});
                    return draw(v,xyz,pairs,color);
                };
                renderInterpolationFrame=()=>{window.trajectoryOverlays=[];return paint()};
            }''')
            for mi, record in enumerate(case['records'], 1):
                frame.locator(f'#mech-sel button[data-id="{mi}"]').click()
                assert frame.locator('#interp-panel').is_visible()
                actual = frame.evaluate('findMech(currentMechId)')
                mapping = dict(enumerate(record['mapping']))
                assert actual['mapping_RP'] == {str(a):b for a,b in mapping.items()}
                path = actual['endpoint_interpolation']['frames']
                assert len(path) == 101
                coords = np.array([f['coords'] for f in path])
                assert np.isfinite(coords).all()
                np.testing.assert_allclose(coords[0], case['endpoints'][0]['coordinates'], atol=1e-10)
                np.testing.assert_allclose(coords[-1], actual['product_xyz_in_R_aligned'], atol=1e-10)
                # Every final pair distance must match the actual mapped product.
                product = np.array(case['endpoints'][1]['coordinates'])[record['mapping']]
                np.testing.assert_allclose(np.linalg.norm(coords[-1,:,None]-coords[-1,None,:], axis=-1),
                    np.linalg.norm(product[:,None]-product[None,:], axis=-1), atol=1e-10)
                expected = sorted(('red' if e['kind'] in ('broken','weakened') else 'green', tuple(sorted(e['r'])))
                                  for e in record['events'])
                for index in (0, 25, 50, 75, 100):
                    frame.evaluate('(i)=>trajectoryAPI.setFrame(i)', index)
                    drawn = frame.evaluate('window.trajectoryOverlays')
                    assert sorted((e['color'],tuple(sorted(e['pair']))) for e in drawn) == expected
                    displayed = frame.evaluate('interpViewer.selectedAtoms({}).map(a=>[a.x,a.y,a.z])')
                    np.testing.assert_allclose(displayed, coords[index], atol=1e-10)
                frame.evaluate('window.trajectoryOverlays=[]')
                frame.locator('#interpEvent').select_option('0')
                assert len(frame.evaluate('window.trajectoryOverlays')) == 1
                frame.locator('#interpEvent').select_option('all')
                frame.evaluate('trajectoryAPI.setFrame(0)')
                frame.locator('#interpNext').click()
                assert frame.evaluate('trajectoryAPI.getState().frame') == 1
                frame.locator('#interpPrev').click()
                assert frame.evaluate('trajectoryAPI.getState().frame') == 0
                frame.locator('#interpPlayBtn').click()
                frame.wait_for_function('trajectoryAPI.getState().frame >= 2')
                frame.locator('#interpPlayBtn').click()
                assert not frame.evaluate('trajectoryAPI.getState().playing')
                frame.evaluate('trajectoryAPI.setFrame(50)')
                with page.expect_download() as info:
                    frame.locator('#interpDownload').click()
                download = info.value
                target = Path(downloads) / download.suggested_filename
                download.save_as(target)
                lines = target.read_text().splitlines()
                n = len(mapping)
                assert len(lines) == 101 * (n + 2)
                exported = np.array([[list(map(float,line.split()[1:])) for line in lines[i+2:i+n+2]]
                                     for i in range(0,len(lines),n+2)])
                np.testing.assert_allclose(exported, coords, atol=5.1e-7)
                if case['index'] == 101 and mi in (1, 2):
                    frame.locator('#interp-panel').screenshot(path=str(OUT / f'pr7-trajectory-{mi}.png'))
                checks.append(dict(case=case['index'], mechanism=mi, frames=101,
                    endpoint_identity='passed', event_overlays='passed', playback='passed', xyz_download='passed'))
                print(f'Checked case {case["index"]}, mechanism {mi}: 101 frames, endpoints, colors, playback, XYZ.', flush=True)
        page.set_viewport_size({'width':390,'height':844})
        assert page.evaluate('document.documentElement.scrollWidth') <= 390
        assert page.frames[1].evaluate('document.documentElement.scrollWidth') <= 390
        assert not errors, errors
        assert not external, external
        browser.close()
    (OUT / 'trajectory-checks.json').write_text(json.dumps(dict(status='passed',
        checks=checks, javascript_errors=errors, external_requests=external), indent=2)+'\n')


if __name__ == '__main__':
    main()
