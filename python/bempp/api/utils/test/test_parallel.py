"""Test cases for parallel utils"""
from unittest import TestCase
import bempp
import numpy as np


class TestCalculateBlocks(TestCase):
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


class TestGatherResults(TestCase):

    def test_gatherresults(self):
        """Test that gatherresults works as expected for a simple example"""
        shape = (10, 10)
        engines = (2, 2)
        self.__gathertest(shape, engines)

    def test_gatherresults_asym(self):
        """Test that gatherresults works as expected for a simple example"""
        shape = (20, 10)
        engines = (4, 2)
        self.__gathertest(shape, engines)

    def __gathertest(self, shape, engines):
        from bempp.api.utils.parallel import gatherresults

        numperengine = (shape[0]//engines[0], shape[1]//engines[1])
        arraysize = shape[0]*shape[1]
        expectedarray = np.arange(arraysize).reshape(shape)
        splitarray = []
        for i in range(engines[0]):
            for j in range(engines[1]):
                splitarray.append(
                expectedarray[i*numperengine[0]:(i+1)*numperengine[0],
                              j*numperengine[1]:(j+1)*numperengine[1]])

        viewmock = dict()
        viewmock['arrayname'] = splitarray
        resultarray = gatherresults(viewmock, 'arrayname',
                                    engines[0], engines[1])

        np.testing.assert_array_equal(resultarray, expectedarray)

if __name__ == "__main__":
    from unittest import main

    main()
