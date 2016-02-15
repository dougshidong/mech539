#!/home/ddong/anaconda2/bin/python
'''Ploting results'''
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages

# Latex Font
plt.rc('text', usetex=True)
plt.rc('font', family='serif')


rfname = 'RESULTS1100'

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

#pp=PdfPages(r'./report/Figures/U400.pdf')
#plt.title(r'U distribution')
#plt.xlabel(r'X')
#plt.ylabel(r'Y')
#plt.axis('equal')
#plt.axis([0, 1, 0, 1])
#plt.gca().set_aspect('equal', adjustable='box')
#plt.draw()
#CS = plt.contourf(X, Y, U, 10)
#pp.savefig()
#pp.close()

it= data[lbound]
lbound += 1
res = data[lbound:lbound+it]
lbound += it

pp=PdfPages('./report/Figures/Conv100.pdf')
plt.title(r'Convergence N = 100')
plt.xlabel(r'Iterations')
plt.ylabel(r'Residual')
iterations = np.linspace(1, it, it)
print iterations
print it
plt.semilogy(iterations, res, '-')
pp.savefig()
pp.close()

grids = [100, 200, 400]
for inx in grids:
    plt.figure()
    for isolv in range(1, 4):
        rfname = 'RESULTS%d%d' %(isolv, inx)
        lbound = 0
        data = np.loadtxt(rfname)
        
        ni   = data[lbound]
        lbound += 1
        nj   = data[lbound]
        lbound += 1
        U = data[lbound:lbound+ni*nj]
        lbound += ni*nj
        it= data[lbound]
        lbound += 1
        res = data[lbound:lbound+it]
        lbound += it
        plt.semilogy(range(1,int(it) + 1), res)
    pp=PdfPages('./report/Figures/Conv%d.pdf' % inx)
    pp.savefig()
    pp.close()
