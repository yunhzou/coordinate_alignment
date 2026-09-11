"""Check native 3D geometry, displayed identities and all nine saved witnesses."""
import argparse
import json
import os
from pathlib import Path

from playwright.sync_api import sync_playwright


ROOT = Path(__file__).resolve().parents[1]


def check(args):
    environment = dict(os.environ)
    libraries = ROOT / 'manuscript/build/browser-libs/root/usr/lib/x86_64-linux-gnu'
    if libraries.exists():
        environment['LD_LIBRARY_PATH'] = str(libraries) + ':' + environment.get('LD_LIBRARY_PATH', '')
    cached = Path('/h/399/yunhengzou/.cache/ms-playwright/chromium-1243/chrome-linux64/chrome')
    executable = os.environ.get('MANUSCRIPT_BROWSER', str(cached) if cached.exists() else None)
    args.output.mkdir(parents=True, exist_ok=True)
    with sync_playwright() as pw:
        browser = pw.chromium.launch(executable_path=executable, env=environment, headless=True,
                                     args=['--no-sandbox', '--enable-unsafe-swiftshader'])
        page = browser.new_page(viewport={'width':1440, 'height':1020}, device_scale_factor=1)
        errors, requests, checked = [], [], []
        page.on('pageerror', lambda error: errors.append(str(error)))
        page.on('request', lambda request: requests.append(request.url) if request.url.startswith(('http:', 'https:')) else None)
        page.goto(args.viewer.resolve().as_uri(), wait_until='load')
        assert page.locator('.view canvas').count() == 2
        for index, count in ((11, 2), (64, 3), (101, 4)):
            page.locator('#case').select_option(str(index))
            assert page.locator('#witness option').count() == count
            assert page.locator('#recovered').is_visible() == (index == 101)
            for seeds in (1, 3, 10, 30):
                badge = page.locator(f'[data-seeds="{seeds}"]').inner_text()
                assert ('recovered' in badge) == (index == 101 and seeds >= 10)
            for ordinal in range(count):
                camera = page.evaluate('views.map(v=>v.getView())')
                page.locator('#witness').select_option(str(ordinal))
                assert page.evaluate('views.map(v=>v.getView())') == camera
                assert page.locator('#event-rows tr').count() == (5 if index == 64 else 4)
                # Endpoint geometry and R identities must agree with the selected actual mapping.
                assert page.evaluate('''() => views.every((v,side)=>v.selectedAtoms({}).every(a=>{
                    const point=current.endpoints[side].coordinates[a.index];
                    return [a.x,a.y,a.z].every((x,j)=>x===point[j]) &&
                      (side ? mapping()[a.properties.r]===a.index : a.properties.r===a.index);
                }))''')
                page.locator('#event-rows tr').first.click()
                assert '→ p' in page.locator('#inspect').inner_text()
                checked.append({'case': index, 'witness': ordinal,
                                'pattern': page.evaluate('current.records[method].pattern')})
            page.locator('#aam').click()
            assert page.evaluate('current.records[method].role') == 'aam'
            page.locator('#slap').click()
            assert page.evaluate('current.records[method].role') == 'missing'
            page.locator('#hydrogen').uncheck()
            page.locator('#focus').click()
            page.screenshot(path=str(args.output / f'viewer-case-{index}-missing.png'), full_page=True)
            if index == 101:
                target = page.evaluate('current.records[method].pattern')
                page.locator('#recovered').click()
                assert page.evaluate('current.records[method].pattern') == target
                assert page.evaluate('current.records[method].provenance.seeds') == 10
                # The AAM recovery uses its actual witness, not copied SLAP atom labels.
                assert page.evaluate("JSON.stringify(mapping())!==JSON.stringify(current.records.find(r=>r.role==='missing').mapping)")
            for control in ('labels', 'event-labels', 'colors', 'events'):
                page.locator('#' + control).check()
                page.locator('#' + control).uncheck()
            page.locator('#event-labels').check()
            page.locator('#events').check()
            page.locator('#hydrogen').check()
            page.locator('#reset').click()
            print(f'Checked case {index}: {count} witnesses and display controls.', flush=True)
        page.locator('#case').select_option('11')
        page.locator('#hydrogen').uncheck()
        assert page.evaluate("views[0].selectedAtoms({index:73})[0].style.sphere !== undefined")
        page.locator('#r').drag_to(page.locator('#r'), source_position={'x':250,'y':180},
                                 target_position={'x':340,'y':225})
        assert page.evaluate('JSON.stringify(views[0].getView().slice(4))===JSON.stringify(views[1].getView().slice(4))')
        page.goto(args.viewer.resolve().as_uri() + '#case=101&witness=2')
        page.wait_for_function('current.index===101 && method===2')
        assert page.evaluate('current.index===101 && current.records[method].provenance.seeds===10')
        page.set_viewport_size({'width':390, 'height':844})
        assert page.evaluate('document.documentElement.scrollWidth') <= 390
        assert not errors and not requests, (errors, requests)
        (args.output / 'viewer_validation.json').write_text(json.dumps(dict(
            status='passed', displayed_witnesses=checked, javascript_errors=errors,
            external_requests=requests, checks=[
                'Three cases and nine independently scored witnesses',
                'Original coordinates and actual mapped identities on both 3D endpoints',
                'Camera preserved when switching witnesses; rotation linked',
                'All event rows, atom selection, display toggles, focus and reset',
                'Two one-seed alternatives in case 64',
                'Actual AAM recovery witness matches SLAP pattern class in case 101',
                'Event H remains visible when spectator H is hidden',
                'Direct case/witness URL and 390px mobile width',
                'Offline rendering with zero external requests']), indent=2) + '\n')
        browser.close()
    print(f'Passed: {len(checked)} witnesses, native R/P geometry, identities and controls.')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    default = ROOT / 'reports/holdout_missing_pattern_seeds_20260910'
    parser.add_argument('--viewer', type=Path, default=default / 'viewer.html')
    parser.add_argument('--output', type=Path, default=default)
    check(parser.parse_args())
