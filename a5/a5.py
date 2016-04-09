#!/home/ddong/anaconda2/bin/python
import scipy.interpolate as inter
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

plt.rc('text', usetex=True)
plt.rc('font', family='serif')


tempf = 'temp.pdf'

# Flow Conditions
Re = 2e6
Min = 0.1
pstat = 1.0
gam = 1.4

# Grid Indices
imax = 513  # Size of I-index 1-based
jmax = 257  # Size of J-index 1-based
tes = 64    # TE Lower I-index 0-based
tee = 448   # TE Upper I-index 0-based
le = 256    # LE Index I-index 0-based
fname = 'naca0012_result.dat'

# Load Solution File
#data = np.loadtxt(fname, dtype = np.float64, unpack=False)
#np.save('bindata', data)
data = np.load('bindata.npy')

# Get Variables
x = np.reshape(data[:, 0], (imax, jmax), order='F')
y = np.reshape(data[:, 1], (imax, jmax), order='F')
rho = np.reshape(data[:, 2], (imax, jmax), order='F')
rhou = np.reshape(data[:, 3], (imax, jmax), order='F')
rhov = np.reshape(data[:, 4], (imax, jmax), order='F')
rhoe = np.reshape(data[:, 5], (imax, jmax), order='F')
mu = np.reshape(data[:, 6], (imax, jmax), order='F')
edv = np.reshape(data[:, 7], (imax, jmax), order='F')
d2 = np.reshape(data[:, 8], (imax, jmax), order='F')

U = np.empty([3, imax, jmax], dtype=np.float64)
U[0,:,:] = rhou/rho
U[1,:,:] = rhov/rho

p = (gam - 1.0) * (rhoe - (rhou**2.0 + rhov**2.0) / (2.0 * rho))
pw = (p[:, 0] + p[:, 1]) / 2.0
cpw = (pw/pstat - 1.0) / (2.0 * gam * Min**2.0)

def rotate_vec(xComp, yComp, t):
    return (np.cos(t) * xComp - np.sin(t) * yComp,
            np.sin(t) * xComp + np.cos(t) * yComp)

# Velocity Gradient (i, j, (u,v), (x,y))
jtau = 1;
dUdX = np.empty([imax, jtau, 2, 2], dtype=np.float64)
# Shear Stress (i, j, (x,y), (x,y))
tau = np.empty([imax, jtau, 2, 2], dtype=np.float64)
# Shear Stress Normal and Axial
tauNA = np.empty([imax, jtau, 2], dtype=np.float64)
# Shear Stress at the Wall (tangent to the surface)
tauW = np.empty([imax, jtau], dtype=np.float64)
# Angle theta at the surface
theta = np.empty([imax, jtau], dtype=np.float64)
# Surface Normals
normals = np.empty([imax, jtau, 2], dtype=np.float64)
# Surface Area
area = np.empty([imax, jtau], dtype=np.float64)

for i in range(imax-1):
    for j in range(jtau):
        dxdxi  = x[i+1,j] - x[i,j]
        dydxi  = y[i+1,j] - y[i,j]
        dxdeta = x[i,j+1] - x[i,j]
        dydeta = y[i,j+1] - y[i,j]
        J = 1 / (dxdxi*dydeta - dydxi*dxdeta)

        area[i,j] = np.sqrt(dxdxi**2.0 + dydxi**2.0)

        # Angle theta at the surface
        theta[i,j] = np.arctan(dydxi/dxdxi)
        # Surface Normals (90 degree rotation)
        normals[i,j,:] = rotate_vec(dxdxi, dydxi, np.pi/2.0) / area[i,j]

        # Velocity Gradient (i, j, (u,v), (x,y))
        for ui in range(2):
            dudxi  = U[ui, i+1, j] - U[ui, i, j]
            dudeta = U[ui, i, j+1] - U[ui, i, j]

            dUdX[i, j, ui, 0] = J * (dydeta * dudxi - dydxi * dudeta)
            dUdX[i, j, ui, 1] = J * (dxdxi * dudeta - dxdeta * dudxi)

        # Shear Stress (i, j, (x,y), (x,y))
        muSum = mu[i,j] + edv[i,j]
        # tau XX
        tau[i,j,0,0] = muSum * 2.0 * dUdX[i,j,0,0] \
                       + 2.0/3.0 * muSum * (dUdX[i,j,0,0] + dUdX[i,j,1,1])
        # tau YY
        tau[i,j,1,1] = muSum * 2.0 * dUdX[i,j,1,1] \
                       + 2.0/3.0 * muSum * (dUdX[i,j,0,0] + dUdX[i,j,1,1])
        # tau XY
        tau[i,j,0,1] = muSum * (dUdX[i,j,0,1] + dUdX[i,j,1,0])
        # tau YX
        tau[i,j,1,0] = tau[i,j,0,1]

        # tau Normal
        tauNA[i,j,0] = normals[i,j,0] * tau[i,j,0,1] \
                       + normals[i,j,1] * tau[i,j,1,1]
        # tau Axial
        tauNA[i,j,1] = normals[i,j,0] * tau[i,j,0,0] \
                       + normals[i,j,1] * tau[i,j,0,1]

        # tau at the wall (tangent to the wall)

        tauW[i,j] = tauNA[i,j,0] * np.cos(theta[i,j]) + tauNA[i,j,1] * np.sin(theta[i,j])


Ca_visc = sum(tauNA[tes:tee,0,1] * area[tes:tee,0]) \
          / (0.5 * gam * pstat * Min**2)

Cn_visc = sum(tauNA[tes:tee,0,0] * area[tes:tee,0]) \
          / (0.5 * gam * pstat * Min**2)

#Ca_pres = sum(-cpw[tes:tee] * normals[tes:tee,0,0] * area[tes:tee,0])
#Cn_pres = sum(-cpw[tes:tee] * normals[tes:tee,0,1] * area[tes:tee,0])
Ca_pres = sum(cpw[tes:tee] * np.sin(theta[tes:tee,0]) * area[tes:tee,0])
Cn_pres = sum(-cpw[tes:tee] * np.cos(theta[tes:tee,0]) * area[tes:tee,0])

print 'Ca_visc: ', Ca_visc
print 'Cn_visc: ', Cn_visc

print 'Ca_pres: ', Ca_pres
print 'Cn_pres: ', Cn_pres

print '\n'
[Cd_visc, Cl_visc] = rotate_vec(Ca_visc, Cn_visc, 8.0/180.0 * np.pi)
[Cd_pres, Cl_pres] = rotate_vec(Ca_pres, Cn_pres, 8.0/180.0 * np.pi)

print 'Cd_visc: ', Cd_visc
print 'Cl_visc: ', Cl_visc

print 'Cd_pres: ', Cd_pres
print 'Cl_pres: ', Cl_pres

# Skin Friction Coefficient
Cv = tauW / (0.5 *  gam * pstat * Min**2)

dUdX[-1, :, :, :] = 0;

if 1 == 1:
    pp = PdfPages(tempf)
    plt.figure(figsize = (6,6))
    for j in range(100):
        plt.plot(x[0:imax,j], y[0:imax,j],'-k')
    for i in range(imax):
        plt.plot(x[i,0:jmax], y[i,0:jmax],'-k')
    plt.xlim([1,1.0001])
    plt.ylim([0,0.0001])

#   plt.plot(range(tes,tee),
#            -cpw[tes:tee] * normals[tes:tee,0,0] * area[tes:tee,0],'-sb',ms=2)

#   plt.plot(range(tes,tee),normals[tes:tee,:,0],'-sb',ms=2)

#   plt.plot(range(tes,tee),theta[tes:tee, 0],'-sb',ms=2)
    plt.grid()

    pp.savefig(bbx_inches='tight')
    pp.close()

def question1():
    fname = 'report/Figures/q1.pdf'
    pp = PdfPages(fname)
    plt.figure(figsize = (6,6))

    plt.plot(x[tes:le, 0], cpw[tes:le], '-sb',
             ms=2, label = 'Lower Surface')
    plt.plot(x[le:tee, 0], cpw[le:tee], '-sr',
             ms=2, label = 'Upper Surface')
    plt.gca().invert_yaxis()
    plt.grid()

    plt.legend(loc='upper right')

    plt.title(r'Pressure Distribution')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$C_p$')

    pp.savefig(bbx_inches='tight')
    pp.close()

def question2():
    fname = tempf
#    fname = 'report/Figures/q2.pdf'
    pp = PdfPages(fname)
    plt.figure(figsize = (6,6))

    plt.plot(x[tes:le, 0], Cv[tes:le], '-sb',
             ms=2, label = 'Lower Surface')
    plt.plot(x[le:tee, 0], Cv[le:tee], '-sr',
             ms=2, label = 'Upper Surface')

    plt.grid()
    plt.legend(loc='upper right')

    plt.title(r'Skin Friction Coefficient')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$C_f$')

    pp.savefig(bbx_inches='tight')
    pp.close()

def uplus_yplus(index_i, max_j, fname):
    ut = rotate_vec(U[0,:,:], U[1,:,:], theta[index_i,0])
    print ut[0]
    uplus = ut[0][index_i, 0:max_j] \
            / np.sqrt( abs(tauW[index_i, 0]) / rho[index_i, 0] )

    yplus = np.sqrt(d2[index_i, 0:max_j]) \
            * np.sqrt( abs(tauW[index_i, 0]) / rho[index_i, 0] ) \
            / ((mu[index_i,0] + edv[index_i,0])/rho[index_i, 0])

    pp = PdfPages(fname)
    plt.figure(figsize = (6,6))

    plt.semilogx(yplus, uplus, '-sb',
             ms=2, label = 'Lower Surface')

    plt.grid()

    plt.title(r'Velocity Profile at I = %d' %index_i)
    plt.xlabel(r'$y^+$')
    plt.ylabel(r'$u^+$')

    pp.savefig(bbx_inches='tight')
    pp.close()
    return


def question6():
    pp = PdfPages(tempf)

    f, axarr = plt.subplots(5,sharex=True, figsize=(8,12))

    w = rhou
    xlocs = [-10, 0.25, 0.5, 1.0, 2.0]
    dx_init = 5.0e-6
    nbp = 90
    yend = 1

    # Overplot x = -10
    xloc = -10.0
    ystart = 0.0
    xn = xloc * np.ones([nbp * 2])
    yn = np.empty([nbp * 2])
    base = ((yend - ystart) / dx_init) ** (1.0/(nbp-1))
    for j in range(nbp):
        yn[j] = -(ystart + dx_init * base**j)
        yn[nbp + j] = ystart + dx_init * base**j

    wslice = inter.griddata((x.flat,
                            y.flat),
                            w.flat,
                            (xn, yn),
                            method='cubic')
    for (xi, xloc) in enumerate(xlocs):
        ax = axarr[xi]
        ax.set_title(r'Plane $x$ = %3.2f' % xloc)
        ax.plot(yn[0:nbp-1], wslice[0:nbp-1], '-s',
                ms=2, label='Lower Farfield')
        ax.plot(yn[nbp:2*nbp-1], wslice[nbp:2*nbp-1], '-o',
                ms=2, label='Upper Farfield')
        ax.grid()

    for (xi, xloc) in enumerate(xlocs[1:]):
        ystart = 0.0
        if xloc > 0.0 and xloc < 1.0:
            ystart = 5.0*0.12*1.0*(
                        0.2969 * np.sqrt(xloc)
                        - 0.1260 * xloc
                        - 0.3516 * xloc**2
                        + 0.2843 * xloc**3
                        - 0.1015 * xloc**4)
        xn = xloc * np.ones([nbp * 2])
        yn = np.empty([nbp * 2])
        base = ((yend - ystart) / dx_init) ** (1.0/(nbp-1))
        for j in range(nbp):
            yn[j] = -(ystart + dx_init * base**j)
            yn[nbp + j] = ystart + dx_init * base**j

        wslice = inter.griddata((x.flat,
                                y.flat),
                                w.flat,
                                (xn, yn),
                                method='cubic')
        ax = axarr[xi+1]
        ax.plot(yn[0:nbp-1], wslice[0:nbp-1], '-s',
                ms=2, label = 'Lower Slice')
        ax.plot(yn[nbp:2*nbp-1], wslice[nbp:2*nbp-1], '-o',
                ms=2, label = 'Upper Slice')


    plt.legend(loc=4,prop={'size':6})

    plt.suptitle(r'Momentum Profiles')
    plt.xlabel(r'$y$')
    plt.ylabel(r'$u$')

    pp.savefig(bbx_inches='tight')
    pp.close()

    return



#question1()
#question2()
#uplus_yplus(350, 90, 'temp.pdf')
question6()

