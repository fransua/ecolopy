

def rec_mult1 (la,lb):
    new = []
    len_a = len (la)
    len_b = len (lb)
    max_len = len_a + len_b
    if len_a > len_b:
        lb = lb + [0] * (len_a - len_b)
        len_b = len_a
    else:
        la = la + [0] * (len_b - len_a)
        len_a = len_b
    for i in xrange ((len_a + len_b - 1) / 2):
        new.append (sum ([la[j]*lb[i-j] for j in xrange (i)]))
    for i in xrange (len_a - 1, max_len - 1):
        new.append (sum ([la[j]*lb[i-j] for j in xrange (1+i-len_a, len_a)]))
    return new

def rec_mult2 (la,lb):
    len_a = len (la)
    len_b = len (lb)
    new = [0] * (len_a + len_b - 1)
    for i in xrange (len_a):
        for j in xrange (len_b):
            new[i+j] += la[i]*lb[j]
    return new

def main():

    la = range (5000)
    lb = range (5000)

    c = rec_mult1(la,lb)
    d = rec_mult2(la,lb)

    print c
    print d
