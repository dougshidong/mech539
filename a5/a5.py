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
aoa = 8.0

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
pw = (p[0:-1, 0] + p[1:, 0]) / 2.0
cpw = (pw/pstat - 1.0) / (0.5 * gam * Min**2.0)

def rotate_vec(xComp, yComp, t):
    return (np.cos(t) * xComp - np.sin(t) * yComp,
            np.sin(t) * xComp + np.cos(t) * yComp)

# Velocity Gradient (i, (u,v), (x,y))
dUdX = np.empty([imax, 2, 2], dtype=np.float64)
# Shear Stress (i, (x,y), (x,y))
tau = np.empty([imax, 2, 2], dtype=np.float64)
# Shear Stress Normal and Axial
tauNA = np.empty([imax, 2], dtype=np.float64)
# Shear Stress at the Wall (tangent to the surface)
tauW = np.empty([imax], dtype=np.float64)
# Angle theta at the surface
theta = np.empty([imax], dtype=np.float64)
# Surface Normals
normals = np.empty([imax, 2], dtype=np.float64)
# Surface Area
area = np.empty([imax], dtype=np.float64)


stag_index = np.logical_and(x[:le,0] < 0.5, U[0,:le,1] < 0.0).argmax()
for i in range(imax-1):
    dxdxi  = x[i+1,0] - x[i,0]
    dydxi  = y[i+1,0] - y[i,0]
    dxdeta = x[i,0+1] - x[i,0]
    dydeta = y[i,0+1] - y[i,0]
    J = 1 / (dxdxi*dydeta - dydxi*dxdeta)

    area[i] = np.sqrt(dxdxi**2.0 + dydxi**2.0)

    # Angle theta at the surface
    theta[i] = np.arctan2(dydxi,dxdxi)
    if(i < stag_index):
        theta[i] += np.pi

    # Surface Normals (90 degree rotation)
    normals[i,:] = rotate_vec(dxdxi, dydxi, np.pi/2.0) / area[i]

    # Velocity Gradient (i, (u,v), (x,y))
    for ui in range(2):
        dudxi  = U[ui, i+1, 0] - U[ui, i, 0]
        dudeta = U[ui, i,   1] - U[ui, i, 0]

        dUdX[i, ui, 0] = J * (dydeta * dudxi - dydxi * dudeta)
        dUdX[i, ui, 1] = J * (dxdxi * dudeta - dxdeta * dudxi)

    # Shear Stress (i, (x,y), (x,y))
    muSum = mu[i,0] + edv[i,0]
    # tau XX
    tau[i,0,0] = muSum * 2.0 * dUdX[i,0,0] \
                   + 2.0/3.0 * muSum * (dUdX[i,0,0] + dUdX[i,1,1])
    # tau YY
    tau[i,1,1] = muSum * 2.0 * dUdX[i,1,1] \
                   + 2.0/3.0 * muSum * (dUdX[i,0,0] + dUdX[i,1,1])
    # tau XY
    tau[i,0,1] = muSum * (dUdX[i,0,1] + dUdX[i,1,0])
    # tau YX
    tau[i,1,0] = tau[i,0,1]

    # tau Normal
    tauNA[i,0] = normals[i,0] * tau[i,0,1] \
                   + normals[i,1] * tau[i,1,1]
    # tau Axial
    tauNA[i,1] = normals[i,0] * tau[i,0,0] \
                   + normals[i,1] * tau[i,0,1]

    # tau at the wall (tangent to the wall)

    tauW[i] = tauNA[i,1] * np.cos(theta[i]) + tauNA[i,0] * np.sin(theta[i])



Ca_visc = sum(tauNA[tes:tee,1] * area[tes:tee]) \
          / (0.5 * gam * pstat * Min**2)

Cn_visc = sum(tauNA[tes:tee,0] * area[tes:tee]) \
          / (0.5 * gam * pstat * Min**2)

Ca_pres = sum(-cpw[tes:tee] * normals[tes:tee,0] * area[tes:tee])
Cn_pres = sum(-cpw[tes:tee] * normals[tes:tee,1] * area[tes:tee])

print 'Ca_visc: ', Ca_visc
print 'Cn_visc: ', Cn_visc

print 'Ca_pres: ', Ca_pres
print 'Cn_pres: ', Cn_pres

print '\n'

alpha = aoa/180.0 * np.pi
Cd_visc = Cn_visc * np.sin(alpha) + Ca_visc * np.cos(alpha)
Cl_visc = Cn_visc * np.cos(alpha) - Ca_visc * np.sin(alpha)
Cd_pres = Cn_pres * np.sin(alpha) + Ca_pres * np.cos(alpha)
Cl_pres = Cn_pres * np.cos(alpha) - Ca_pres * np.sin(alpha)

print 'Cd_visc: ', Cd_visc
print 'Cl_visc: ', Cl_visc

print 'Cd_pres: ', Cd_pres
print 'Cl_pres: ', Cl_pres

# Skin Friction Coefficient
Cv = tauW / (0.5 *  gam * pstat * Min**2)

if 0 == 1:
    pp = PdfPages(tempf)
    plt.figure(figsize = (6,6))
#   for j in range(100):
#       plt.plot(x[0:imax,j], y[0:imax,j],'-k')
#   for i in range(imax):
#       plt.plot(x[i,0:jmax], y[i,0:jmax],'-k')
#   plt.xlim([1,1.0001])
#   plt.ylim([0,0.0001])
#    plot(theta

#   plt.plot(range(tes,tee),
#            -cpw[tes:tee] * normals[tes:tee,0,0] * area[tes:tee,0],'-sb',ms=2)

#   plt.plot(range(tes,tee),normals[tes:tee,:,0],'-sb',ms=2)

    plt.plot(range(tes,tee),theta[tes:tee, 0],'-sb',ms=2)
    plt.grid()

    pp.savefig(bbx_inches='tight')
    pp.close()

def question1():
    fname = tempf
#    fname = 'report/Figures/q1.pdf'
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

    plt.tight_layout()
    pp.savefig(bbx_inches='tight')
    pp.close()

def question4():
    fname = tempf
#    fname = 'report/Figures/q3airfoil.pdf'
    pp = PdfPages(fname)
    plt.figure(figsize = (5,2))

    upper_lamsep = np.logical_and(x[le:tee,0] < 0.5, Cv[le:tee] < 0.0).argmax() + le

    upper_turbstart = np.logical_and(x[upper_lamsep:tee,0] < 0.5,
                                    Cv[upper_lamsep:tee] > 0.0).argmax() + upper_lamsep

    upper_turbsep = np.logical_and(x[upper_turbstart:tee,0] > 0.5,
                                 Cv[upper_turbstart:tee] < 0.0).argmax() + upper_turbstart

    upper_sep = np.logical_and(x[upper_turbsep:tee,0] > 0.5,
                              Cv[upper_turbsep:tee] > 0.0).argmax() + upper_turbsep

    lower_lamsep = le - (Cv[le:tes:-1] < 0.0).argmax()
    lower_turbstart = lower_lamsep - (Cv[lower_lamsep:tes:-1] > 0.0).argmax()

    print Cv[le:tee]
    print tes
    print lower_turbstart
    l1 = plt.plot(x[lower_lamsep:upper_lamsep, 0],
             y[lower_lamsep:upper_lamsep, 0],
             '-b', label='Laminar Flow')

    l2 = plt.plot(x[tes:lower_turbstart, 0],
             y[tes:lower_turbstart, 0],
             '-r', label='Turbulent Flow')
    plt.plot(x[upper_turbstart:upper_sep, 0],
             y[upper_turbstart:upper_sep, 0],
             '-r', label='_Turbulent Flow')

    l3 = plt.plot(x[lower_lamsep:lower_turbstart:-1,0],
             y[lower_lamsep:lower_turbstart:-1,0],
             '-c', label='Recirculation')
    plt.plot(x[upper_lamsep:upper_turbstart,0],
             y[upper_lamsep:upper_turbstart,0],
             '-c', label='_Recirculation')
    plt.plot(x[upper_turbsep:upper_sep,0],
             y[upper_turbsep:upper_sep,0],
             '-c', label='_Recirculation')

#   d1 = plt.plot(x[lower_lamsep,0],y[lower_lamsep,0],
#                 'sk', ms=2, label = 'Laminar Separation')
#   plt.plot(x[upper_lamsep,0],y[upper_lamsep,0],
#            'sk', ms=2, label = 'Laminar Separation')
#   d2 = plt.plot(x[upper_turbsep,0],y[upper_turbsep,0],
#                 '*k', ms = 2,label = 'Turbulent Separation')
    plt.annotate('Laminar Separation', xy=(x[lower_lamsep,0],y[lower_lamsep,0]),
                fontsize = 7,
                xytext=(x[lower_lamsep,0]-0.5,y[lower_lamsep,0]-0.08),
                arrowprops=dict(arrowstyle="->")
               )
    plt.annotate('Laminar Separation', xy=(x[upper_lamsep,0]-0.01,y[upper_lamsep,0]),
                fontsize = 7,
                xytext=(x[upper_lamsep,0]-0.05,y[upper_lamsep,0]+0.07),
                arrowprops=dict(arrowstyle="->")
               )

    plt.annotate('Turbulent Separation', xy=(x[upper_turbsep,0],y[upper_turbsep,0]),
                fontsize = 7,
                xytext=(x[upper_turbsep,0]-0.4,y[upper_turbsep,0]+0.07),
                arrowprops=dict(arrowstyle="->")
               )


    plt.legend(loc=4,prop={'size':5})

    plt.title(r'Laminar and Turbulent Regions')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$y$')

    plt.axis('equal')
    plt.xlim([-0.1,1.1])

#   legend1 = plt.legend([l1,l2,l3],["test1","test2","test3"], loc=4,prop={'size':5})
#   plt.legend([d1, d2],["tes1","tes2"], loc=1,prop={'size':5})
#   plt.gca().add_artist(legend1)

    plt.tight_layout()

    pp.savefig(bbx_inches='tight')
    pp.close()
def uplus_yplus(index_i, max_j, fname):
    ut = np.abs(rotate_vec(U[0,:,:], U[1,:,:], theta[index_i]))

    uplus = ut[0][index_i, 0:max_j] \
            / np.sqrt( abs(tauW[index_i]) / rho[index_i, 0] )

    yplus = np.sqrt(d2[index_i, 0:max_j]) \
            * np.sqrt( abs(tauW[index_i]) / rho[index_i, 0] ) \
            / ((mu[index_i,0] + edv[index_i,0])/rho[index_i, 0])

    pp = PdfPages(fname)
    plt.figure(figsize = (6,6))

    plt.semilogx(yplus, yplus, '-k',
             ms=2, label = 'Lower Surface')
    plt.semilogx(yplus, np.log(yplus) / 0.41 + 5.15, '-k',
             ms=2, label = 'Lower Surface')
    plt.semilogx(yplus, uplus, '-sb',
             ms=2, label = 'Lower Surface')

    plt.grid()

    print uplus
    plt.ylim([0,max(uplus)])
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
                                method='linear')
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
question4()
uplus_yplus(151, 100, 'temp.pdf')
uplus_yplus(362, 100, 'temp.pdf')
#question6()

