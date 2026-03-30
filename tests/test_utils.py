from pathlib import Path

import pytest

from flukit.utils import utils


@pytest.fixture
def setup():
    test_dir = Path(__file__).parent
    example_genbank = test_dir / "test-data" / "h1n1.gb"
    return example_genbank


def test_load_features(setup):
    example_genbank = setup
    features = utils.load_features(example_genbank)

    assert isinstance(features, dict)
