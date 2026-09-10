"""Exercise the offline viewer and save review screenshots."""
import json
import os
from pathlib import Path
from playwright.sync_api import sync_playwright

MAN=Path(__file__).resolve().parents[1]
OUT=MAN/'build';OUT.mkdir(exist_ok=True)
with sync_playwright() as pw:
    cached=Path('/h/399/yunhengzou/.cache/ms-playwright/chromium-1243/chrome-linux64/chrome')
    executable=os.environ.get('MANUSCRIPT_BROWSER',str(cached) if cached.exists() else None)
    local_libs=OUT/'browser-libs/root/usr/lib/x86_64-linux-gnu'
    environment=dict(os.environ)
    if local_libs.exists():
        environment['LD_LIBRARY_PATH']=str(local_libs)+':'+environment.get('LD_LIBRARY_PATH','')
    browser=pw.chromium.launch(
        executable_path=executable,env=environment,headless=True,args=['--no-sandbox'])
    page=browser.new_page(viewport={'width':1280,'height':1040},device_scale_factor=1)
    errors=[];page.on('pageerror',lambda e:errors.append(str(e)))
    page.goto((MAN/'animations/index.html').as_uri(),wait_until='load')
    assert page.locator('#matched').inner_text()=='1 / 18'
    assert page.locator('#candidates').inner_text()=='5'
    page.locator('#next').click()
    assert page.evaluate('animationAPI.getState().step')==1
    page.locator('#steps').select_option('6')
    assert page.locator('#candidates').inner_text()=='42'
    page.locator('#witness').select_option('2')
    page.screenshot(path=str(OUT/'viewer-growth.png'),full_page=True)
    page.locator('#hydrogens').uncheck();page.locator('#labels').uncheck()
    page.locator('#zoom').select_option('1.4')
    page.locator('#r').drag_to(page.locator('#r'),source_position={'x':120,'y':100},target_position={'x':190,'y':160})
    page.locator('#reset').click()
    page.locator('#hydrogens').check();page.locator('#labels').check()
    page.locator('#speed').select_option('300');page.locator('#play').click()
    page.wait_for_timeout(730)
    page.locator('#play').click()
    assert page.evaluate('animationAPI.getState().step')>=8
    # Click an actual terminal node, using the rendered graph's hit coordinates.
    pt=page.evaluate('graphHit[0]');page.locator('#graph').click(position={'x':pt['x'],'y':pt['y']})
    assert page.locator('#matched').inner_text()=='18 / 18'
    assert 'Retained terminal' in page.locator('#title').inner_text()
    page.screenshot(path=str(OUT/'viewer-terminal.png'),full_page=True)
    page.locator('#case').select_option('1')
    assert page.evaluate('animationAPI.getState().frames')==70
    deferred=page.evaluate("DATA[1].frames.findIndex(f=>f.event==='consumed')")
    page.locator('#steps').select_option(str(deferred))
    assert 'Defer boundary' in page.locator('#title').inner_text()
    page.screenshot(path=str(OUT/'viewer-deferred.png'),full_page=True)
    page.set_viewport_size({'width':390,'height':844})
    assert page.evaluate('document.documentElement.scrollWidth')<=390
    page.screenshot(path=str(OUT/'viewer-mobile.png'),full_page=True)
    assert not errors,errors
    (OUT/'viewer-validation.json').write_text(json.dumps({
        'status':'passed','javascript_errors':errors,
        'checks':['offline load','event stepping','candidate witness selection','hydrogen/label toggles',
                  'coordinate rotation/zoom/reset','timed playback','terminal graph click',
                  'both case traces','deferred-edge event','390px mobile width']},indent=2)+'\n')
    browser.close()
print('Offline viewer checks passed.')
