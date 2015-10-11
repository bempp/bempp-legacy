"""Test cases for parallel utils"""
from unittest import TestCase
import bempp
import numpy as np


class TestCalculateBlocks(TestCase):
    """Test cases for parallel utils"""

    def test_calculateblocks(self):
        """Test that calulate blocks works in a simple case"""
        size = (100, 100)
        nengines = 4
        expectedengines = (2, 2)
        expectedchunks = self.__expected4enginesplit(size[0])
        self.__calculateblocks_test(size, nengines, expectedengines,
                                    expectedchunks)

    def test_calculateblocks_reminder(self):
        """Test that calulate blocks works when the matrix is not dividable with
           the number of engines"""
        size = (117, 117)
        nengines = 4
        expectedengines = (2, 2)
        expectedchunks = self.__expected4enginesplit(size[0])
        self.__calculateblocks_test(size, nengines, expectedengines,
                                    expectedchunks)

    def test_calculateblocks_nonsquare(self):
        """Test correct split in a non square matrix"""
        size = (150, 100)
        nengines = 4
        expectedengines = (4, 1)
        row1 = (0, 37)
        row2 = (37, 74)
        row3 = (74, 111)
        row4 = (111, 150)
        column = (0, 100)
        expectedchunks = [(row1, column),
                          (row2, column),
                          (row3, column),
                          (row4, column)]
        self.__calculateblocks_test(size, nengines, expectedengines,
                                    expectedchunks)

    def test_calculateblocks_nonsquare2(self):
        """Test correct split in a non square matrix"""
        size = (400, 117)
        nengines = 8
        expectedengines = (8, 1)
        row1 = (0, 50)
        row2 = (50, 100)
        row3 = (100, 150)
        row4 = (150, 200)
        row5 = (200, 250)
        row6 = (250, 300)
        row7 = (300, 350)
        row8 = (350, 400)
        column = (0, 117)
        expectedchunks = [(row1, column),
                          (row2, column),
                          (row3, column),
                          (row4, column),
                          (row5, column),
                          (row6, column),
                          (row7, column),
                          (row8, column)]
        self.__calculateblocks_test(size, nengines, expectedengines,
                                    expectedchunks)

    def __calculateblocks_test(self, size, engines, expectedengines,
                               expectedchunks):
        from bempp.api.utils.parallel import calculateblocks
        chunks, nrowseng, ncolseng = calculateblocks(size[0], size[1], engines)
        assert nrowseng == expectedengines[0]
        assert ncolseng == expectedengines[1]
        np.testing.assert_array_equal(chunks, expectedchunks)

    def __expected4enginesplit(self, size):
        lowrow = (0, size//2)
        highrow = (size//2, size)
        lowcolumn = lowrow
        highcolumn = highrow
        expectedchunks = [(lowrow, lowcolumn),
                          (lowrow, highcolumn),
                          (highrow, lowcolumn),
                          (highrow, highcolumn)]
        return expectedchunks


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
