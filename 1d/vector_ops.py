def local_average(v, offset, W):
    """Given a vector -- iterable -- v of length N and a filter W, compute the average of v
    with W as a filter, centered on v[offset].  So local_average(v, W, 3)"""
        

    start = offset - len(W)//2    
       
    l = list(W[i] * v[start + i] for i in range(len(W)))

    return sum(l)/sum(W)
