const {chromium}=require(process.env.PLAYWRIGHT_MODULE||'playwright');
(async()=>{
 const browser=await chromium.launch({headless:true,args:['--no-sandbox','--enable-unsafe-swiftshader']});
 const page=await browser.newPage({viewport:{width:1500,height:1100}});
 const errors=[],requests=[];page.on('pageerror',e=>errors.push(e.message));
 page.on('request',r=>{if(/^https?:/.test(r.url()))requests.push(r.url())});
 await page.goto('file://'+process.argv[2]);
 console.log('Loaded offline viewer');
 for(let i=0;i<3;i++){
  await page.locator('#case').selectOption(String(i));
  if(await page.locator('.view canvas').count()!==4)throw Error('Missing canvas');
  await page.locator('[data-event]').first().click();
  if(!(await page.locator('#inspect').textContent()).includes('AAM →'))throw Error('Selection failed');
  await page.locator('#hydrogen').uncheck();await page.locator('#hydrogen').check();
  await page.locator('#labels').check();await page.locator('#labels').uncheck();
  console.log('Checked case '+i);
 }
 await page.locator('#case').selectOption('0');
 await page.screenshot({path:process.argv[3],fullPage:true});
 if(errors.length||requests.length)throw Error(JSON.stringify({errors,requests}));
 console.log('PASS: three cases, four canvases, selection and toggles; zero JS errors or external requests.');
 await browser.close();
})().catch(e=>{console.error(e);process.exit(1)});
