// ══════════════════════════════════════════════════════════════════════════════
// Mineralogical Spectral Indices – Harmonized Landsat & Sentinel-2 (HLS v2)
// ──────────────────────────────────────────────────────────────────────────────
// Indices mapped:
//   1. Iron Oxide Index        (IOI)  – Red / Blue
//   2. Clay Alteration Index   (CAI)  – SWIR1 / SWIR2
//   3. Quartz Index            (QI)   – SWIR2 / SWIR1  |  TIR2/TIR1 (Landsat)
//   4. Hydrothermal Alteration (HTI)  – (SWIR1 + Red) / (NIR + Blue)
//
// Data: NASA/HLS/HLSL30/v002 (Landsat 30 m) + NASA/HLS/HLSS30/v002 (S2 30 m)
// Both collections are atmospherically corrected surface reflectance.
// ══════════════════════════════════════════════════════════════════════════════


// ── 1. AOI – uploaded shapefile asset ────────────────────────────────────────
var aoi = ee.FeatureCollection('projects/tavadaisheprojects/assets/MASVINGO');

Map.centerObject(aoi, 9);
Map.addLayer(aoi, {color: 'white'}, 'Masvingo AOI');


// ── 2. Cloud mask using HLS Fmask band ───────────────────────────────────────
function maskHLS(img) {
  var fmask = img.select('Fmask');
  var clear = fmask.bitwiseAnd(1 << 1).eq(0)   // bit 1: cloud
    .and(fmask.bitwiseAnd(1 << 2).eq(0))        // bit 2: adjacent cloud
    .and(fmask.bitwiseAnd(1 << 3).eq(0))        // bit 3: cloud shadow
    .and(fmask.bitwiseAnd(1 << 4).eq(0));       // bit 4: snow/ice
  return img.updateMask(clear)
            .divide(10000)                       // scale to 0–1 reflectance
            .copyProperties(img, img.propertyNames());
}


// ── 3. Load HLSL30 (Landsat) and harmonise band names ────────────────────────
//    HLSL30 original → common name
//    B02 Blue  | B03 Green | B04 Red | B05 NIR
//    B06 SWIR1 | B07 SWIR2 | B10 TIR1 | B11 TIR2

var hlsl = ee.ImageCollection('NASA/HLS/HLSL30/v002')
  .filterBounds(aoi)
  .filterDate('2020-01-01', '2022-12-31')
  .filter(ee.Filter.lt('CLOUD_COVERAGE', 20))
  .map(maskHLS)
  .select(
    ['B02', 'B03', 'B04', 'B05', 'B06', 'B07', 'B10',  'B11'],
    ['Blue','Green','Red', 'NIR','SWIR1','SWIR2','TIR1','TIR2']
  );

print(hlsl, 'HLSL30 collection');


// ── 4. Load HLSS30 (Sentinel-2) and harmonise band names ─────────────────────
//    HLSS30 original → common name
//    B02 Blue  | B03 Green | B04 Red | B8A NIR (narrow)
//    B11 SWIR1 | B12 SWIR2
//    (no thermal bands in S2)

var hlss = ee.ImageCollection('NASA/HLS/HLSS30/v002')
  .filterBounds(aoi)
  .filterDate('2020-01-01', '2022-12-31')
  .filter(ee.Filter.lt('CLOUD_COVERAGE', 20))
  .map(maskHLS)
  .select(
    ['B02', 'B03', 'B04', 'B8A', 'B11',  'B12'],
    ['Blue','Green','Red', 'NIR','SWIR1','SWIR2']
  );

print(hlss, 'HLSS30 collection');


// ── 5. Merged median composite (optical bands common to both sensors) ─────────
var composite = hlsl.select(['Blue','Green','Red','NIR','SWIR1','SWIR2'])
  .merge(hlss.select(['Blue','Green','Red','NIR','SWIR1','SWIR2']))
  .median()
  .clip(aoi);

print(composite, 'HLS Merged Composite');

// True colour preview
Map.addLayer(composite,
  {bands: ['Red','Green','Blue'], min: 0, max: 0.3},
  'True Colour – HLS Composite');

// SWIR False Colour (good for geology: SWIR2–SWIR1–Red)
Map.addLayer(composite,
  {bands: ['SWIR2','SWIR1','Red'], min: 0, max: 0.4},
  'SWIR False Colour – Geology');


// ── 6. Landsat-only composite (needed for thermal Quartz Index) ───────────────
var compositeL = hlsl.median().clip(aoi);


// ── 7. Compute spectral indices ───────────────────────────────────────────────

// 7a. Iron Oxide Index (IOI) ──────────────────────────────────────────────────
//   IOI = Red / Blue
//   High values → iron-bearing minerals (goethite, haematite, limonite)
//   Useful for mapping oxidised ore deposits and ferruginous crusts.
var IOI = composite.select('Red').divide(composite.select('Blue'))
  .rename('IOI');

Map.addLayer(IOI,
  {min: 0.5, max: 3.0, palette: ['#ffffcc','#fd8d3c','#800026']},
  'Iron Oxide Index (IOI)');


// 7b. Clay Alteration Index (CAI) ─────────────────────────────────────────────
//   CAI = SWIR1 / SWIR2
//   High values → clay / Al-OH / Mg-OH minerals (kaolinite, montmorillonite,
//   sericite, chlorite) associated with hydrothermal and weathering zones.
var CAI = composite.select('SWIR1').divide(composite.select('SWIR2'))
  .rename('CAI');

Map.addLayer(CAI,
  {min: 0.6, max: 1.6, palette: ['#f7f7f7','#9ecae1','#08306b']},
  'Clay Alteration Index (CAI)');


// 7c-i. Quartz Index – SWIR ratio (works for both Landsat & Sentinel-2) ───────
//   QI_SWIR = SWIR2 / SWIR1
//   Quartz-rich rocks absorb at ~8.6 µm; this SWIR proxy highlights silica-
//   enriched zones (quartzite, silicified veins, felsic intrusives).
var QI_SWIR = composite.select('SWIR2').divide(composite.select('SWIR1'))
  .rename('QI_SWIR');

Map.addLayer(QI_SWIR,
  {min: 0.6, max: 1.6, palette: ['#edf8e9','#74c476','#00441b']},
  'Quartz Index – SWIR (QI_SWIR)');


// 7c-ii. Quartz Index – Thermal (Landsat only) ────────────────────────────────
//   QI_TIR = TIR1 / TIR2  (Band 10 / Band 11)
//   Quartz emits strongly at TIR1 relative to TIR2; high ratio = quartz-rich.
//   Landsat TIRS bands are 100 m native, resampled to 30 m in HLS.
var QI_TIR = compositeL.select('TIR1').divide(compositeL.select('TIR2'))
  .rename('QI_TIR');

Map.addLayer(QI_TIR,
  {min: 0.92, max: 1.08, palette: ['#313695','#ffffbf','#a50026']},
  'Quartz Index – Thermal/TIR (QI_TIR) [Landsat only]');


// 7d. Hydrothermal Alteration Index (HTI) ─────────────────────────────────────
//   HTI = (SWIR1 + Red) / (NIR + Blue)
//   Combines iron sensitivity (Red) with OH-mineral sensitivity (SWIR1);
//   elevated values indicate hydrothermally altered rocks, gossans, and
//   epithermal/porphyry alteration halos.
var HTI = composite.select('SWIR1').add(composite.select('Red'))
  .divide(composite.select('NIR').add(composite.select('Blue')))
  .rename('HTI');

Map.addLayer(HTI,
  {min: 0.5, max: 2.5, palette: ['#fff5eb','#fd8d3c','#7f0000']},
  'Hydrothermal Alteration Index (HTI)');


// ── 8. Stack all indices into one image for inspection & export ───────────────
var indexStack = IOI
  .addBands(CAI)
  .addBands(QI_SWIR)
  .addBands(HTI);

print(indexStack, 'Index Stack (IOI | CAI | QI_SWIR | HTI)');

// Optional: stack thermal QI with the rest (Landsat footprint only)
var indexStackFull = IOI
  .addBands(CAI)
  .addBands(QI_SWIR)
  .addBands(QI_TIR)
  .addBands(HTI);

print(indexStackFull, 'Full Index Stack incl. Thermal QI');


// ── 9. Export to Google Drive ─────────────────────────────────────────────────
// Each index exported as a separate single-band GeoTIFF, plus a stacked file.
// CRS: EPSG:32736 – UTM Zone 36S (covers Masvingo, Zimbabwe)

var exportCRS    = 'EPSG:32736';   // UTM Zone 36S – Masvingo, Zimbabwe
var exportScale  = 30;             // HLS native resolution
var exportRegion = aoi;

Export.image.toDrive({
  image: IOI.toFloat(), description: 'IOI_Iron_Oxide_Index',
  folder: 'GEE_SpectralIndices', crs: exportCRS,
  region: exportRegion, scale: exportScale, maxPixels: 1e13
});

Export.image.toDrive({
  image: CAI.toFloat(), description: 'CAI_Clay_Alteration_Index',
  folder: 'GEE_SpectralIndices', crs: exportCRS,
  region: exportRegion, scale: exportScale, maxPixels: 1e13
});

Export.image.toDrive({
  image: QI_SWIR.toFloat(), description: 'QI_Quartz_Index_SWIR',
  folder: 'GEE_SpectralIndices', crs: exportCRS,
  region: exportRegion, scale: exportScale, maxPixels: 1e13
});

Export.image.toDrive({
  image: QI_TIR.toFloat(), description: 'QI_Quartz_Index_Thermal_Landsat',
  folder: 'GEE_SpectralIndices', crs: exportCRS,
  region: exportRegion, scale: exportScale, maxPixels: 1e13
});

Export.image.toDrive({
  image: HTI.toFloat(), description: 'HTI_Hydrothermal_Alteration_Index',
  folder: 'GEE_SpectralIndices', crs: exportCRS,
  region: exportRegion, scale: exportScale, maxPixels: 1e13
});

// Stacked export (all 4 optical indices in one multi-band file)
Export.image.toDrive({
  image: indexStack.toFloat(), description: 'All_Indices_Stacked_HLS',
  folder: 'GEE_SpectralIndices', crs: exportCRS,
  region: exportRegion, scale: exportScale, maxPixels: 1e13
});
