import hashlib


def calculate_md5(path, blocksize=65535):
    """
    Calculates the MD5 on a file at the given path
    """
    hasher = hashlib.md5()

    with open(path, 'rb') as fileobj:
        buf = fileobj.read(blocksize)
        while len(buf) > 0:
            hasher.update(buf)
            buf = fileobj.read(blocksize)

    return hasher.hexdigest()
