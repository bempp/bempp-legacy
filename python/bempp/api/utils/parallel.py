"Some helper functions for parallel calculations"


def calculateblocks(nrows, ncols, numengines):
    """"Given the x and y size of an array split it in upto numengine chunks"""

    # we have ncoleng engines along the xaxis and nroweng engines
    # along the y axis with ncoleng *nroweng <= numengines
    nroweng, ncoleng = _calculateenginesplit(nrows, ncols, numengines)
    chunksizerows = nrows // nroweng
    chunksizecols = ncols // ncoleng
    blocks = list()
    for i in range(nroweng):
        for j in range(ncoleng):
            rowsstart = chunksizerows*i
            if i == nroweng-1:
                rowsend = nrows
            else:
                rowsend = chunksizerows*(i+1)
            colsstart = chunksizecols*j
            if j == ncoleng-1:
                colsend = ncols
            else:
                colsend = chunksizecols*(j+1)
            blocks.append(((rowsstart, rowsend), (colsstart, colsend)))
    return blocks, nroweng, ncoleng


def _calculateenginesplit(nrows, ncols, numengines):
    "Private helper function for calculateblocks"
    import math
    aspectratio = nrows/ncols

    nsmalldim = int(math.floor(numengines**(1./2)/aspectratio)) or 1
    nlargedim = int(math.floor(numengines/nsmalldim))
    if aspectratio > 1:
        nroweng = nlargedim
        ncoleng = nsmalldim
    else:
        nroweng = nsmalldim
        ncoleng = nlargedim
    return (nroweng, ncoleng)


def gatherresults(view, resultname, nrows, ncols):
    """Gather results from engines"""
    import numpy as np
    partial = [np.concatenate(view[resultname][i*ncols:(i+1)*ncols],
                              axis=1) for i in range(nrows)]
    result = np.concatenate(partial, axis=0)
    return result
