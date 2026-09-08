// Offline browser QA; PLAYWRIGHT_MODULE points to the installed test runner.
const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
const fs=require('fs');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox']});
 const page=await browser.newPage({viewport:{width:1800,height:1200}});
 const errors=[],requests=[];page.on('pageerror',e=>errors.push(e.message));
 page.on('request',r=>{if(/^https?:/.test(r.url()))requests.push(r.url())});
 await page.goto('file://'+process.argv[2],{waitUntil:'load'});
 const count=await page.locator('#case option').count();if(count!==21)throw new Error('Expected 21 cases');
 const checked=[];
 for(let i=0;i<count;i++){
  await page.locator('#case').selectOption(String(i));
  if(await page.locator('.panel svg').count()!==4)throw new Error('Missing R/P drawings at '+i);
  const choices=await page.locator('#candidate option').count();
  for(let j=1;j<=choices;j++){
   await page.locator('#candidate').selectOption(String(j));
   if(await page.locator('#detected .atom-hit').count()===0)throw new Error('Empty detected diagram');
  }
  await page.locator('#candidate').selectOption('1');
  const heavy=await page.locator('#detected .atom-hit').count();
  await page.locator('#hydrogens').check();
  const all=await page.locator('#detected .atom-hit').count();
  if(all<=heavy)throw new Error('H toggle failed '+i);
  await page.locator('#hydrogens').uncheck();
  await page.locator('#link').selectOption('p');
  await page.locator('#detected .atom-hit[data-side="P"]').first().click();
  if(!await page.locator('#reference .selected').count()||!await page.locator('#detected .selected').count())throw new Error('Linked atom highlighting failed');
  checked.push({position:i,choices,heavy,all});
 }
 await page.locator('#case').selectOption('10');
 await page.evaluate(()=>window.scrollTo(0,0));
 await page.screenshot({path:process.argv[3]+'.png',fullPage:false});
 await page.setViewportSize({width:900,height:1000});
 const dimensions=await page.locator('.comparison').boundingBox();
 if(dimensions.x+dimensions.width>901)throw new Error('Page overflow');
 if(errors.length||requests.length)throw new Error(JSON.stringify({errors,requests}));
 fs.writeFileSync(process.argv[3]+'.json',JSON.stringify({cases:count,checked,errors,externalRequests:requests},null,2));
 await browser.close();console.log('All 21 cases, candidates, H toggles and linked selections passed offline.');
})().catch(e=>{console.error(e);process.exit(1)});
