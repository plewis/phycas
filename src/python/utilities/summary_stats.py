import sys,math

def summarize(v):
    """
    Returns basic summary statistics (mean, variance, stddev) for the supplied vector v in a tuple.
    
    """
    assert v is not None
    n = len(v)
    assert n > 1, 'summarize function requires a vector of at least two values'
    sum = 0.0
    sumsq = 0.0
    for x in v:
        sum += x
        sumsq += x**2.0
    mean = sum/float(n)
    meansq = mean**2.0
    variance = (sumsq - float(n)*meansq)/(float(n) - 1.0)
    stddev = math.sqrt(variance)
    return (mean, variance, stddev)    

