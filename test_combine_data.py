"""Test module.

Run test with:

    pytest test_combine_data

from the root of the repo.

It assumes data is stored in ../Pakistan
"""

from pathlib import Path
from combine_data import (
    read_data_pakistan,
    read_flexpart_univie,
    read_lagranto_chc,
    read_flexpart_xu,
    read_flexpart_tatfancheng,
    read_wam2layers,
    read_utrack,
    read_tracmass,
    read_uvigo,
    read_btrims,
    read_ughent,
    read_2ldrm,
)


def test_read_data_pakistan():
    """Check if function completes with no error."""
    # Assume data is stored relative to this file like so:
    folder = "../Pakistan/"

    # Test data read function
    read_data_pakistan(folder, "Pakistan")


def test_individual_datasets():
    basedir = Path("../")

    for casename in ["Pakistan", "Australia", "Scotland"]:
        read_flexpart_univie(basedir, casename)
        read_lagranto_chc(basedir, casename)
        read_flexpart_xu(basedir, casename)
        read_flexpart_tatfancheng(basedir, casename)
        read_wam2layers(basedir, casename)
        read_utrack(basedir, casename)
        read_tracmass(basedir, casename)
        read_uvigo(basedir, casename)
        read_btrims(basedir, casename)
        read_ughent(basedir, casename)
        read_2ldrm(basedir, casename)
