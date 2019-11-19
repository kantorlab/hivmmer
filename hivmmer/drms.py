"""
"""
import pandas as pd
from importlib import resources

def _load_drm_list():
    """
    Load a pre-packaged list of IAS and Stanford DRMs.
    """
    with resources.open_binary("hivmmer", "drms.csv") as f:
        return pd.read_csv(f)

def drms(aafile, outfile, coverage=10, frequency=0.01):
    """
    """
    # Load the AA table
    aa = pd.read_excel(aafile)
    del aa["hxb2"]

    # Unpivot and remove zero- or low-coverage calls
    aa = aa.melt(id_vars=["region", "position", "coverage"], var_name="variant", value_name="count")
    aa = aa[aa["count"] >= 10]

    # Retain calls above the frequency threshold
    aa["frequency"] = aa["count"] / aa["coverage"]
    aa = aa[aa["frequency"] >= frequency]
    del aa["coverage"]
    del aa["count"]

    # Inner join the DRM list and output
    aa = aa.merge(_load_drm_list(), on=["region", "position", "variant"], how="inner")
    aa.to_csv(outfile, index=False)

# vim: expandtab sw=4 ts=4
