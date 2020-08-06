import numpy as np

def eformat(f, prec, exp_digits):
    s = "%.*e"%(prec, f)
    mantissa, exp = s.split('e')
    # add 1 to digits as 1 is taken by sign +/-
    return "%se%+0*d"%(mantissa, exp_digits+1, int(exp))

def h_format(f):
    return eformat(f,6,3).rjust(14)

def set_gw_depth(depth,fn,incr=0.01):
    X_RANGE=slice(6,20)
    H_RANGE=slice(21,35)
    BEFORE_H=slice(0,21)
    AFTER_H=slice(35,None)
    with open(fn) as fp:
        existing = fp.readlines()

    header = existing[:5]
    body = existing[5:-2]
    footer = existing[-2:]
    x_s = np.array([float(ln[X_RANGE]) for ln in body])
    delta = np.abs(x_s-depth)
    ix = np.argmin(delta)
    # make Hs and substitute
    h_s = [ln[H_RANGE] for ln in body]
    h0 = x_s[ix]
    new_h_s = np.arange(h0,h0+incr*len(h_s),incr)
    new_body = [ln[BEFORE_H]+h_format(h)+ln[AFTER_H] for ln,h in zip(body,new_h_s)]

    with open(fn,'w') as fp:
        fp.write(''.join(header + new_body + footer))
