"Some helper functions for parallel calculations"


def calculateblocks(nrows, ncols, numengines):
    """"Given the x and y size of an array split it in upto numengine chunks"""

    # we have ncoleng engines along the xaxis and nroweng engines
    # along the y axis with ncoleng *nroweng > numengines
    # nchunks = nroweng*ncoleng
    ncoleng, nroweng = _calculateenginesplit(nrows, ncols, numengines)
    chunksizerows = nrows // nroweng
    chunksizecols = ncols // ncoleng
    blocks = list()
    for i in range(nroweng):
        for j in range(ncoleng):
            colsstart = chunksizecols*i
            if i == nroweng-1:
                colsend = ncols
            else:
                colsend = chunksizecols*(i+1)
            rowsstart = chunksizerows*j
            if j == ncoleng-1:
                rowsend = nrows
            else:
                rowsend = chunksizerows*(j+1)
            blocks.append(((colsstart, colsend), (rowsstart, rowsend)))
    return blocks, nroweng, ncoleng


def _calculateenginesplit(nrows, ncols, numengines):
    "Private helper function for calculateblocks"
    import math
    aspectratio = nrows//ncols
    if aspectratio > 1.1:
        raise NotImplementedError
    elif aspectratio < 0.9:
        raise NotImplementedError
    else:
        nroweng = int(math.floor(numengines**(1./2)))
        ncoleng = nroweng
    return (ncoleng, nroweng)


def gatherresults(view, resultname, nrows, ncols):
    """Gather results from engines"""
    import numpy as np
    partial = [np.concatenate(view[resultname][i*ncols:(i+1)*ncols],
                              axis=1) for i in range(nrows)]
    result = np.concatenate(partial, axis=0)
    return result
