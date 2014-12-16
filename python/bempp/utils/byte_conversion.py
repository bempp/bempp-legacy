

def convert_to_bytes(s):
    res = s
    try:
        if not isinstance(s,bytes):
            res = res.encode('UTF-8')
    except:
        raise ValueError('String type expected.')
    return res
