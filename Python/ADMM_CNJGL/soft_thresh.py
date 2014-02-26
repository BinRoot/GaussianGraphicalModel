import math

def soft_thresh(y, l):
    x = y
    if y == 0 :
        return 0
    else :
        return math.copysign(1, y) * max(abs(y) - l, 0)
