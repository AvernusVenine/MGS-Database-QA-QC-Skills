"""Read-only loading of the QDI texture view through its ArcGIS layer file.

Only arcpy.da.SearchCursor is used. The layer file on disk is never saved.

The published layer uses `relateid` as its object-ID field. That field is not unique (one borehole has many
samples), so ArcGIS quietly drops some rows (71,396 of 71,903 on 2026-10-15). The in-memory copy of the layer is
switched to `objectid`, which is unique in qdtx, so that every row is read.
"""
import os

import pandas as pd

DEFAULT_LYR = r'L:\gdrs\data\pub\edu_umn_mngs\geos_quaternary_data_index\Texture Data.lyr'


def open_layer(lyr_path):
    import arcpy
    layer = arcpy.mp.LayerFile(lyr_path).listLayers()[0]
    cim = layer.getDefinition('V3')
    conn = cim.featureTable.dataConnection
    published_oid = conn.oIDFields
    if published_oid != 'objectid':
        conn.oIDFields = 'objectid'
        layer.setDefinition(cim)          # in-memory only; LayerFile.save() is never called
    return layer, published_oid


def load(lyr_path=DEFAULT_LYR, cache=None, refresh=False):
    """DataFrame of every row of the texture view (+ x, y), and facts about the layer for the report."""
    if cache and not refresh and os.path.exists(cache):
        return pd.read_pickle(cache), pd.read_pickle(cache + '.meta')
    import arcpy
    raw_layer = arcpy.mp.LayerFile(lyr_path).listLayers()[0]
    published_count = sum(1 for _ in arcpy.da.SearchCursor(raw_layer, ['objectid']))
    layer, published_oid = open_layer(lyr_path)
    fields = [f.name for f in arcpy.ListFields(layer) if f.type not in ('Geometry', 'Blob')
              and f.name.upper() != 'ESRI_OID']
    rows = arcpy.da.SearchCursor(layer, fields + ['SHAPE@X', 'SHAPE@Y'])
    df = pd.DataFrame(list(rows), columns=fields + ['x', 'y'])
    meta = pd.Series({'lyr': lyr_path, 'published_oid_field': published_oid,
                      'rows_seen_by_published_layer': published_count, 'rows_in_view': len(df),
                      'max_date_modify': str(pd.to_datetime(df['date_modify'], errors='coerce').max())})
    if cache:
        os.makedirs(os.path.dirname(cache) or '.', exist_ok=True)
        df.to_pickle(cache)
        meta.to_pickle(cache + '.meta')
    return df, meta
