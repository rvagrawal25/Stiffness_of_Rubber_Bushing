import sympy as sym

r = sym.symbols('r', positive=True)
theta = sym.symbols('theta')

delr_delx = sym.cos(theta)

delr_dely = sym.sin(theta)

deltheta_delx = -sym.sin(theta)/r

deltheta_dely = sym.cos(theta)/r

def del_delx(f):
    return delr_delx*sym.diff(f,r) + deltheta_delx*sym.diff(f,theta)

def del_dely(f):
    return delr_dely*sym.diff(f,r) + deltheta_dely*sym.diff(f,theta)

def del2_delx2(f):
    return del_delx(del_delx(f))

def del2_dely2(f):
    return del_dely(del_dely(f))

def polarlaplacian(f):
    return (del2_delx2(f)+del2_dely2(f)).simplify()

def polarbiharm(f):
    return polarlaplacian(polarlaplacian(f))

def sigma_xx(f):
    return del2_dely2(f)
def sigma_yy(f):
    return del2_delx2(f)
def sigma_xy(f):
    return -del_delx(del_dely(f))

def sigma_rect(f):
    return sym.Matrix([[sigma_xx(f), sigma_xy(f)], [sigma_xy(f), sigma_yy(f)]])

Q = sym.Matrix([[sym.cos(theta), sym.sin(theta)],[-sym.sin(theta), sym.cos(theta)]])

def sigma_polar(f):
    return Q*sigma_rect(f)*(Q.T)

def sigma_rr(f):
    return sigma_polar(f)[0,0].simplify().expand()
def sigma_tt(f):
    return sigma_polar(f)[1,1].simplify().expand()
def sigma_rt(f):
    return sigma_polar(f)[0,1].simplify().expand()
def sigma_tr(f):
    return sigma_polar(f)[1,0].simplify().expand()