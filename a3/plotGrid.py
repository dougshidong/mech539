#!/home/ddong/anaconda2/bin/python
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
#symb = ('o', 'v', '<', 's', 'p', '*', 'h', 'H', 'D', 'd')
col=('b', 'g', 'r', 'c', 'm', 'k')

def question2():
    # Coefficient of Pressure on Surface
    pp = PdfPages('./report/Figures/cpsurf.pdf')
    plt.figure(figsize = (9,5))
    i = 0
    for mach in np.linspace(80, 90, 6):
        fname = './q2results/solq2_%d.dat' % mach
        x,y,cp = np.loadtxt(fname, dtype=np.float64, unpack=True)
        imax = 80
        jmax = 40

        x = np.reshape(x, (imax, jmax))
        y = np.reshape(y, (imax, jmax))
        cp = np.reshape(cp, (imax, jmax))

        plt.xlim(19.5, 21.5)

        plt.plot(x[:, 0], cp[:, 0],
                 color = col[i],
                 marker = 'o', mec = col[i], mfc = 'None', ms = 4,
                 label = r'Mach = ' + str(mach/100.0))
        i = i + 1

    plt.gca().invert_yaxis()
    plt.legend(loc='upper right', fontsize = 10)
    plt.title(r'$C_p$ for Varying Freestream Mach')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$C_p$')

    plt.tight_layout()

    pp.savefig(bbox_inches='tight')
    pp.close()
    # Coefficient of Pressure Contour Plot
    for mach in np.linspace(80, 90, 6):
        fname = './q2results/solq2_%d.dat' % mach
        x,y,cp = np.loadtxt(fname, dtype=np.float64, unpack=True)
        imax = 159
        jmax = 79
        imax = 80
        jmax = 40

        x = np.reshape(x, (imax, jmax))
        y = np.reshape(y, (imax, jmax))
        cp = np.reshape(cp, (imax, jmax))

        pp = PdfPages('./report/Figures/cpcontour_%d.pdf' % mach)

        plt.figure(figsize = (10,4))
        plt.xlim(19.5, 21.5)
        plt.ylim(0, 1)
        plt.gca().set_aspect('equal', adjustable='box')
        CS = plt.contour(x, y, cp, 10)

        plt.clabel(CS, CS.levels, inline=True, fmt='%4.3f', fontsize=8)
        plt.title(r'Coefficient of Pressure Distribution, Mach = %.2f' % (mach/100.0))
        plt.xlabel(r'$x$')
        plt.ylabel(r'$y$')

        pp.savefig(bbox_inches='tight')
        pp.close()

    # Convergence Plot
    pp = PdfPages('./report/Figures/convergenceq2.pdf')
    plt.figure(figsize = (9,5))
    i = 0
    for mach in np.linspace(80, 90, 6):
        fname = './q2results/convsolq2_%d.dat' % mach
        resi, times = np.loadtxt(fname, dtype=np.float64, unpack=True)
        plt.semilogy(range(1, len(resi)+1), resi,
                     label = r'Mach = ' + str(mach/100.0))
        #plt.semilogy(range(1, 10000+1), resi[0:10000],
        #             label = r'Mach = ' + str(mach/100.0))
        i = i + 1

    plt.legend(loc='upper right', fontsize = 10)
    plt.title(r'Convergence for Varying Freestream Mach')
    plt.xlabel(r'Iterations')
    plt.ylabel(r'Residual')

    pp.savefig(bbox_inches='tight')
    pp.close()

def question3():
    # Coefficient of Pressure on Surface
    pp = PdfPages('./report/Figures/q3cpsurf.pdf')
    plt.figure(figsize = (9,5))
    i = 0
    mults = [1, 2, 4]
    imaxl = [80, 159, 317]
    jmaxl = [40, 79, 157]

    sname = (r'Coarse', r'Medium', r'Fine')
    for mult in mults:
        fname = './q3results/solq3_%d.dat' % mult
        x,y,cp = np.loadtxt(fname, dtype=np.float64, unpack=True)

        imax = imaxl[i]
        jmax = jmaxl[i]

        x = np.reshape(x, (imax, jmax))
        y = np.reshape(y, (imax, jmax))
        cp = np.reshape(cp, (imax, jmax))

        plt.xlim(19.5, 21.5)

        plt.plot(x[:, 0], cp[:, 0],
                 color = col[i],
                 marker = 'o', mec = col[i], mfc = 'None', ms = 4,
                 label = sname[i])
        i = i + 1

    plt.gca().invert_yaxis()
    plt.legend(loc='upper right', fontsize = 10)
    plt.title(r'$C_p$ for Varying Freestream Grid Sizes')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$C_p$')

    pp.savefig(bbox_inches='tight')
    pp.close()
def question4():
    # Time vs Iterations
    i = 0
    isolv = [2, 3]
    mults = [1, 2]
    imaxl = [80, 159, 317]
    jmaxl = [40, 79, 157]

    sname1 = (r'Gauss-Seidel', r'Line Implicit')
    sname2 = (r'Coarse', r'Medium', r'Fine')

    pp = PdfPages('./report/Figures/q4res.pdf')
    plt.figure(figsize = (9,5))
    for iisolv in range(2):
        i = 0
        for mult in mults:
            fname = ('./q4results/convsolq4_%d%d.dat' % (mult,iisolv+2))
            resi, times = np.loadtxt(fname, dtype=np.float64, unpack=True)
            plt.semilogy(range(1, len(resi)+1), resi,
                         label = (sname1[iisolv],sname2[i]))
            i = i + 1

    plt.legend(loc='upper right', fontsize = 10)
    plt.title(r'Residual Convergence')
    plt.xlabel(r'Iterations')
    plt.ylabel(r'max(\Delta \phi)')

    pp.savefig(bbox_inches='tight')
    pp.close()

    pp = PdfPages('./report/Figures/q4time.pdf')
    plt.figure(figsize = (9,5))
    for iisolv in range(2):
        i = 0
        for mult in mults:
            fname = ('./q4results/convsolq4_%d%d.dat' % (mult,iisolv+2))
            resi, times = np.loadtxt(fname, dtype=np.float64, unpack=True)
            plt.semilogy(times, resi,
                     label = (sname1[iisolv],sname2[i]))
            i = i + 1

    plt.legend(loc='upper right', fontsize = 10)
    plt.title(r'Residual over CPU Time')
    plt.xlabel(r'CPU Time (s)')
    plt.ylabel(r'max(\Delta \phi)')

    pp.savefig(bbox_inches='tight')
    pp.close()

    # Coefficient of Pressure on Surface
    pp = PdfPages('./report/Figures/q4cpsurf.pdf')
    plt.figure(figsize = (9,5))
    i = 0
    for iisolv in range(2):
        fname = './q4results/solq4_1%d.dat' % (iisolv+2)
        x,y,cp = np.loadtxt(fname, dtype=np.float64, unpack=True)
        imax = 80
        jmax = 40

        x = np.reshape(x, (imax, jmax))
        y = np.reshape(y, (imax, jmax))
        cp = np.reshape(cp, (imax, jmax))

        plt.xlim(19.5, 21.5)

        plt.plot(x[:, 0], cp[:, 0],
                 color = col[i],
                 marker = 'o', mec = col[i], mfc = 'None', ms = 4,
                 label = str(iisolv))
        i = i + 1

    plt.gca().invert_yaxis()
    plt.legend(loc='upper right', fontsize = 10)
    plt.title(r'$C_p$ for Varying Freestream Mach')
    plt.xlabel(r'$x$')
    plt.ylabel(r'$C_p$')

    plt.tight_layout()

    pp.savefig(bbox_inches='tight')
    pp.close()
    
    pp = PdfPages('./report/Figures/q4timeres.pdf')
    plt.figure(figsize = (9,5))
    for iisolv in range(2):
        i = 0
        for mult in mults:
            fname = ('./q4results/convsolq4_%d%d.dat' % (mult,iisolv+2))
            resi, times = np.loadtxt(fname, dtype=np.float64, unpack=True)
            plt.plot(times, range(1, len(times)+1),
                     label = (sname1[iisolv],sname2[i]))
            i = i + 1

    plt.legend(loc='upper right', fontsize = 10)
    plt.title(r'Residual over CPU Time')
    plt.xlabel(r'CPU Time (s)')
    plt.ylabel(r'max(\Delta \phi)')

    pp.savefig(bbox_inches='tight')
    pp.close()

#question2()
#question3()
question4()
