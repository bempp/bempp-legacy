"""Test cases for parallel utils"""
from unittest import TestCase
import bempp
import numpy as np


class TestParallel(TestCase):
    """Test cases for parallel utils"""

    def test_calculateblocks(self):
        """Test that calulate blocks works in a simple case"""
        from bempp.api.utils.parallel import calculateblocks
        chunks, nrowseng, ncolseng = calculateblocks(100, 100, 4)

        assert nrowseng == 2
        assert nrowseng == 2
        low = (0, 50)
        high = (50, 100)

        expectedchunks = [(low, low),
                          (low, high),
                          (high, low),
                          (high, high)]
        np.testing.assert_array_equal(chunks, expectedchunks)


    def test_gatherresults(self):
        pass


if __name__ == "__main__":
    from unittest import main

    main()
