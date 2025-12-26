def get_psite(read, offsets):
    """
    Infer P-site position for a read using length-specific offsets.

    Returns:
        int (psite genomic coordinate) or None
    """
    read_len = read.query_length
    if read_len not in offsets:
        return None

    offset = offsets[read_len]

    if read.is_reverse:
        return read.reference_end - offset - 1
    else:
        return read.reference_start + offset
