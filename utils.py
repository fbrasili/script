from numpy import *

def power_law(x, A, B):
    return A * x**B

def calc_vals_err(xs, ys):
    X = sort( list( set(xs) ) )
    Y,dY = array( [ [ ys[argwhere(xs == x)].mean() , ys[argwhere(xs == x)].std()  ] for x in X] ).transpose()
    return array( [X,Y,dY] )

def gen_histogram(x_vals, y_vals, bin_width, computeErrors=True):
    x = arange( 1 , int( x_vals.max() ) , bin_width )
    left_edges = x - bin_width/2.
    cond = array([ len(y_vals[ (x_vals >= edge) * (x_vals <= edge + bin_width) ]) > 0 for edge in left_edges ])
    x = x[cond]
    y = array( [ y_vals[ (x_vals >= edge) * (x_vals <= edge + bin_width) ].mean() for edge in left_edges[cond] ] )
    if computeErrors:
        yerr = array( [ y_vals[ (x_vals >= edge) * (x_vals <= edge + bin_width) ].std() for edge in left_edges[cond] ] )
        return x,y,yerr
    else:
        return x,y

def weighted_mean(vals, weights):
    assert len(vals) == len(weights), "the lengths of the values and of the weights are different"
    return (vals * weights).sum() / weights.sum()

def running_mean(x, window):
    dNi = int( floor(window/2) )
    dNf = window - dNi
    return array([ x[max(0,i-dNi):min(len(x),i+dNf)].mean() for i in range(len(x)) ])

def rewrap(data_configuration, box):
    '''
    shape(data) = (N,M)
        N is the number of atoms
        M>5 is the number of values for each atom
    shape(box)  = (3,2)
    '''
    rangeX,rangeY,rangeZ = box
    lX = rangeX[1] - rangeX[0]
    lY = rangeY[1] - rangeY[0]
    lZ = rangeZ[1] - rangeZ[0]
    for i in arange(len(data_configuration)):
        x,y,z = data_configuration[i,2:5]
        while x < rangeX[0]: x += lX
        while x > rangeX[1]: x -= lX
        while y < rangeY[0]: y += lY
        while y > rangeY[1]: y -= lY
        while z < rangeZ[0]: z += lZ
        while z > rangeZ[1]: z -= lZ
        data_configuration[i,2:5] = x,y,z
    return data_configuration

def recenter(data_configuration):
    cond_mgel = ( data_configuration[:, 1].astype(int) == 1 ) + ( data_configuration[:, 1].astype(int) == 3 )
    data_mgel = data_configuration[cond_mgel, 2:5]
    center = data_mgel.mean(0)
    data_configuration[:, 2:5] = data_configuration[:, 2:5] - center
    return data_configuration

def compute_mean(X, Y, W, dX=None, edges=None, mode="upto"):
    '''
    mode=="upto" --> the mean is computed on x <= x[i]          for x[i]        in x
    mode=="from" --> the mean is computed on x >= x[i]          for x[i]        in x
    mode=="bins" --> the mean is computed on x[i] <= x < x[x+1] for x[i],x[i+1] in x 
    '''
    assert mode in ["upto","from","bins"], "the allowed modes are \"upto\", \"from\" and \"bins\""
    X,Y = array(X),array(Y)
    if dX=="custom":
        assert edges is not None, "specify the width or the edges of bins"
    else:
        edges = arange( X.min()//dX , dX * ( X.max()//dX + 3 )  , dX )
    cond = ( edges >= max(edges[edges<=X.min()]) ) * ( edges <= min(edges[edges> X.max()]) )
    xs = edges[cond]
    if mode == "upto":
        xs  = array([ x for x in xs if any(X<=x) ])
        if W is None:
            ys  = array([ average( Y[X<=x] , weights=W ) for x in xs ])
            dys = array([ sqrt( average( Y[X<=x]**2 , weights=W ) - y**2 ) for x,y in zip(xs,ys) ])
        else:
            ys  = array([ average( Y[X<=x] , weights=WW[X<=x] ) for x in xs ])
            dys = array([ sqrt( average( Y[X<=x]**2 , weights=W[X<=x] ) - y**2 ) for x,y in zip(xs,ys) ])
    elif mode == "from":
        xs  = array([ x for x in xs if any(X>=x) ])
        if W is None:
            ys  = array([ average( Y[X>=x] , weights=W ) for x in xs ])
            dys = array([ sqrt( average( Y[X>=x]**2 , weights=W ) - y**2 ) for x,y in zip(xs,ys) ])
        else:
            ys  = array([ average( Y[X>=x] , weights=W[X>=x] ) for x in xs ])
            dys = array([ sqrt( average( Y[X>=x]**2 , weights=W[X>=x] ) - y**2 ) for x,y in zip(xs,ys) ])
    elif mode == "bins":
        cond_bins = append( array([ sum( (X>=l)*(X<u) ) > 0 for l,u in zip(xs[:-1],xs[1:]) ]) , True )
        xs  = xs[cond_bins]
        if W is None:
            ys  = array([ average(Y[(X>=l)*(X<u)] , weights=W ) for l,u in zip(xs[:-1],xs[1:]) ])
            dys = array([ sqrt( average(Y[(X>=l)*(X<u)]**2 , weights=W ) - y**2 ) for l,u,y in zip(xs[:-1],xs[1:],ys) ])
        else:
            ys  = array([ average(Y[(X>=l)*(X<l+dX)] , weights=W[(X>=l)*(X<u)] ) for l,u in zip(xs[:-1],xs[1:]) ])
            dys = array([ sqrt( average(Y[(X>=l)*(X<l+dX)]**2 , weights=W[(X>=l)*(X<u)] ) - y**2 ) for l,u,y in zip(xs[:-1],xs[1:],ys) ])
        xs = xs[:-1]
    return xs,ys,dys

