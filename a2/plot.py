#!/home/ddong/anaconda2/bin/python
'''Ploting results'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

q3 = 1
q4 = 0
q5 = 0
q6 = 0
q6b = 0
q7 = 0
# Latex Font
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


if(q3 == 1):
    rfname = './q45results/RESULTS1400'
    
    # Read results and assign to variables
    data = np.loadtxt(rfname)
    
    lbound = 0
    # Read results
    ni   = data[lbound]
    lbound += 1
    nj   = data[lbound]
    lbound += 1
    U = data[lbound:lbound+ni*nj]
    lbound += ni*nj
    
    X = np.linspace(0.0, 1.0, ni)
    Y = np.linspace(0.0, 1.0, nj)
    U = np.reshape(U, (ni, nj))

    pp=PdfPages(r'./report/Figures/U400.pdf')
    plt.title(r'U distribution')
    plt.xlabel(r'X')
    plt.ylabel(r'Y')
    plt.axis('equal')
    plt.axis([0, 1, 0, 1])
    plt.gca().set_aspect('equal', adjustable='box')
    CS = plt.contourf(X, Y, U, 10)
    plt.colorbar(CS)
    pp.savefig()
    pp.close()

if (q4 == 1):
    grids = [100, 200, 400]
    solvn = ['Jacobi', 'Gauss-Seidel', 'SOR']
    for inx in grids:
        plt.figure(figsize = (12,5))
        for isolv in range(3):
            rfname = './q45results/RESULTS%d%d' %(isolv + 1, inx)
            lbound = 0
            data = np.loadtxt(rfname)
            
            ni   = data[lbound]
            lbound += 1
            nj   = data[lbound]
            lbound += 1
            lbound += ni*nj
            it= data[lbound]
            lbound += 1
            res = data[lbound:lbound+it]
            lbound += it
            plt.semilogy(range(1,int(it) + 1), res, label = solvn[isolv])
            plt.xlabel(r'Iterations')
            plt.ylabel(r'Residual')
            plt.title(r'Convergence nx = %d' % inx)
        plt.legend()
        pp=PdfPages('./report/Figures/Conv%d.pdf' % inx)
        pp.savefig()
        pp.close()
if (q5 == 1):
    grids = [100, 200, 400]
    solvn = ['Jacobi', 'Gauss-Seidel', 'SOR']
    for inx in grids:
        plt.figure(figsize = (12,5))
        for isolv in range(3):
            rfname = './q45results/RESULTS%d%d' %(isolv + 1, inx)
            lbound = 0
            data = np.loadtxt(rfname)
            
            ni   = data[lbound]
            lbound += 1
            nj   = data[lbound]
            lbound += 1
            lbound += ni*nj
            it= data[lbound]
            lbound += 1
            res = data[lbound:lbound+it]
            lbound += it
            cputime = data[lbound:lbound+it]
            lbound += it
            print cputime
            plt.semilogy(cputime, res, label = solvn[isolv])
            plt.xlabel(r'CPU\_TIME')
            plt.ylabel(r'Residual')
            plt.title(r'Convergence nx = %d' % inx)
        plt.legend()
        pp=PdfPages('./report/Figures/time%d.pdf' % inx)
        pp.savefig()
        pp.close()
if(q6 == 1):
    grids = np.array([50, 100, 200, 300, 400])
    ncoeff = 3600
    calcExact = 0
    relError = np.empty(len(grids))
    solvn = ['Jacobi', 'Gauss-Seidel', 'SOR']
    pp=PdfPages(r'./report/Figures/Error.pdf')
    plt.figure()
    for isolv in range(3):
        print isolv
        for igrid, inx in enumerate(grids):
            print inx
            rfname = './q6results/RESULTS%d%d' %(isolv + 1, inx)
            # Read results and assign to variables
            data = np.loadtxt(rfname)
            
            lbound = 0
            # Read results
            ni   = data[lbound]
            lbound += 1
            nj   = data[lbound]
            lbound += 1
            U = data[lbound:lbound+ni*nj]
            lbound += ni*nj
            
            X = np.linspace(0.0, 1.0, ni, dtype = np.float128)
            Y = np.linspace(0.0, 1.0, nj, dtype = np.float128)
            U = np.reshape(U, (ni, nj))
            Uexact = np.zeros((ni, nj))
            if(calcExact == 1):
                XX, YY = np.meshgrid(X, Y)
                for nn in range(1, ncoeff + 1):
                    nnpi = np.float128(nn * np.pi)
                    cn = -2 * (np.cos(nnpi) - 1) / (nnpi * np.sinh(nnpi))
                    Uexact += cn * np.sinh(nnpi * YY) * np.sin(nnpi * XX)
                np.save('./exactSol/exactSol%d' % inx, Uexact)
            else:
                Uexact = np.load('./exactSol/exactSol%d.npy' % inx)
            Uexact_vector = Uexact[1:-2, 1:-2].flatten()
            U_vector = U[1:-2, 1:-2].flatten()
            relError[igrid] = np.linalg.norm(Uexact_vector - U_vector, 1) \
                              / np.linalg.norm(Uexact_vector, 1)
        plt.semilogy(-1 * np.log10(1.0/grids), relError, '-o', label=solvn[isolv])
        plt.legend()
    plt.title(r'Order of Accuracy')
    plt.xlabel(r'$-log_{10}(\Delta x)$')
    plt.ylabel(r'Error')
    pp.savefig()
    pp.close()
if(q6b == 1):
    grids = [200]
    ncoeff = 3600
    calcExact = 0
    relError = np.empty(len(grids))
    solvn = ['Jacobi', 'Gauss-Seidel', 'SOR']
    isolv = 2
    print isolv
    for igrid, inx in enumerate(grids):
        print inx
        rfname = './q6results/RESULTS%d%d' %(isolv + 1, inx)
        # Read results and assign to variables
        data = np.loadtxt(rfname)
        
        lbound = 0
        # Read results
        ni   = data[lbound]
        lbound += 1
        nj   = data[lbound]
        lbound += 1
        U = data[lbound:lbound+ni*nj]
        lbound += ni*nj
        
        X = np.linspace(0.0, 1.0, ni, dtype = np.float128)
        Y = np.linspace(0.0, 1.0, nj, dtype = np.float128)
        U = np.reshape(U, (ni, nj))
        Uexact = np.zeros((ni, nj))
        if(calcExact == 1):
            XX, YY = np.meshgrid(X, Y)
            for nn in range(1, ncoeff + 1):
                nnpi = np.float128(nn * np.pi)
                cn = -2 * (np.cos(nnpi) - 1) / (nnpi * np.sinh(nnpi))
                Uexact += cn * np.sinh(nnpi * YY) * np.sin(nnpi * XX)
            np.save('./exactSol/exactSol%d' % inx, Uexact)
        else:
            Uexact = np.load('./exactSol/exactSol%d.npy' % inx)
        Uexact_vector = Uexact[0:-1, 0:-1].flatten()
        U_vector = U[0:-1, 0:-1].flatten()
        relError[igrid] = np.linalg.norm(Uexact_vector - U_vector, 1) \
                          / np.linalg.norm(Uexact_vector, 1)
        pp=PdfPages(r'./report/Figures/ErrorContour.pdf')
        plt.title(r' Absolute Error Contour')
        plt.xlabel(r'X')
        plt.ylabel(r'Y')
        plt.gca().set_aspect('equal', adjustable='box')
        CS = plt.contourf(np.abs(X[0:-1]), np.abs(Y[0:-1]), 
                          np.abs(Uexact[0:-1, 0:-1] - U[0:-1,0:-1]),
                          norm = LogNorm())
        plt.colorbar(CS)
        pp.savefig()
        pp.close()
if(q7 == 1):
    grids = [100, 200, 400]
    weights = [0.75, 1.00, 1.25, 1.5, 1.75, 1.90, 1.95, 1.99]
    for inx in grids:
        plt.figure(figsize = (12,5))
        for iw, w in enumerate(weights):
            rfname = './q7results/RESULTS%d%d' %(iw + 1, inx)
            lbound = 0
            data = np.loadtxt(rfname)
            
            ni   = data[lbound]
            lbound += 1
            nj   = data[lbound]
            lbound += 1
            lbound += ni*nj
            it= data[lbound]
            lbound += 1
            res = data[lbound:lbound+it]
            lbound += it
            if(iw != 6):
                plt.semilogy(range(1,int(it) + 1), res, label = r'SOR = %.2f' %w)
        plt.xlabel(r'Iterations')
        plt.ylabel(r'Residual')
        plt.title(r'Effects of Relaxation on SOR Convergence $nx = %d$' % inx)
        plt.legend()
        pp=PdfPages('./report/Figures/SOR%d.pdf' % inx)
        pp.savefig()
        pp.close()
