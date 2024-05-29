"""Test module.

Run test with:

    pytest test_combine_data

from the root of the repo.

It assumes data is stored in ../Pakistan
"""

from combine_data import read_data_pakistan


def test_read_data_pakistan():
    """Check if function completes with no error."""
    # Assume data is stored relative to this file like so:
    folder = "../Pakistan/"

    # Test data read function
    read_data_pakistan(folder)
