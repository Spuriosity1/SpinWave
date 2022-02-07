import numpy as np
from numpy import linalg
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

from rotation import Rx, Ry, Rz

def parameterised_path(path, basis=np.eye(2)):
    coords = np.matmul(np.array(list(path.values())), basis)
    
    x0y0 = coords[0]
    
    # stores the times at which special points appear
    t_list = [0]
    uvec_list = []
    
    for i in range(1,len(path)):
        xy = coords[i]
        arclen = linalg.norm(xy-x0y0)
        t_list.append(t_list[i-1]+arclen)
        uvec_list.append((xy-x0y0)/arclen)
        x0y0 = xy
    
    t_list = np.array(t_list)
    uvec_list = np.array(uvec_list)
    
    def p(t):
        retval = []
        for tau in t:
            idx = t_list.searchsorted(tau)
            if idx==0:
                retval.append(coords[0])
            elif idx < len(path):
                retval.append(coords[idx-1] + (tau-t_list[idx-1])*uvec_list[idx-1])
        return np.array(retval)
            
    return p, t_list


def axisEqual3D(ax):
    extents = np.array([getattr(ax, 'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents, axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers, 'xyz'):
        getattr(ax, 'set_{}lim'.format(dim))(ctr - r, ctr + r)
        
    ax.set_box_aspect([1,1,1])
        
        
        
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class Lattice(object):
    ''' 
    This class understands all information about
    relative locations of atoms, 
    and the rotors on each site
    a1,a2 - primitive LATTICE unit cell vectors
    a1m,a2m - primitive MAGNETIC unit cell vectors
    '''
    def __init__(self, a1,a2,a1m=None,a2m=None):
        
        if (a1m is None or a2m is None):
            a1m = a1
            a2m = a2
        
        self.sites = {}
        self.couplings = {}

        self.natoms = 0
        self.a1 = np.array(a1)
        self.a2 = np.array(a2)
        self.a1m = np.array(a1m)
        self.a2m = np.array(a2m)
        
        self.RS_units = (1., 'arb.')
        
        
    def set_dist_units(self, num, string):
        self.RS_units=(num, string)
    
    @property
    def area(self):
        return np.abs(self.a1[0]*self.a2[1] - self.a1[1]*self.a2[0])
    
    @property
    def mag_area(self):
        return np.abs(self.a1m[0]*self.a2m[1] - self.a1m[1]*self.a2m[0])
        
    @property
    def reciprocal_vectors(self):
        prefactor = 2*np.pi/self.area
        b1 = ( self.a2[1], -self.a2[0])
        b2 = (-self.a1[1],  self.a1[0])
        return prefactor*np.array((b1, b2))
    
    @property
    def mag_reciprocal_vectors(self):
        # returns the reciprocal vectors in the format
        # B[0,:];  B[1,:]
        prefactor = 2*np.pi/self.mag_area
        b1 = ( self.a2m[1], -self.a2m[0])
        b2 = (-self.a1m[1],  self.a1m[0])
        return prefactor*np.array((b1, b2))
        
    # adds an atom
    # works in REAL SPACE units
    def add_atom(self, label, xy, R = np.eye(3)):
#         if (xy[0] > 0.5 or xy[0] < -0.5 or xy[1] > 0.5 or xy[1] < -0.5):
#             print("[WARN] atom is outside unit cell")
        self.sites[label] = ({
            "index": self.natoms,
            "rotor":  R,
            "pos":   np.array(xy)
        })
    
        self.natoms += 1
        
        
    def define_coupling(self, name, matrix=np.zeros((3,3))):
        assert matrix.shape==(3,3)
        if name in self.couplings:
            self.couplings[name]['mat'] = matrix
        else:
            self.couplings[name] = {
                'mat': matrix,
                'bond_list': []
            }
            
            
        
    # Adds a connectivity tensor J between sites with labels "label1", "label2", 
    # r2 - r1 is made to be the closest matching site to the posistion vector label1[xy] + delta
    def add_coupling(self, label1, label2, coupling_name, delta=(0.01,0.009),color='k', check_duplicate=True):
#         assert J.shape == (3,3)
        s1 = self.sites[label1]
        s2 = self.sites[label2]
        
        # find the closest matching site
        delta = np.array(delta)
        
        r = s2['pos'] - s1['pos'] - delta
        possibilities = np.array([
            r, r+self.a1m, r-self.a1m, r+self.a2m, r-self.a2m,
            r+self.a1m+self.a2m, r-self.a1m-self.a2m, r+self.a1m-self.a2m, r-self.a1m+self.a2m])
        
        idx = np.argmin(linalg.norm(possibilities,axis=1))
        
        new_bond = {
            'i': label1,
            'j': label2,
            'rj - ri': possibilities[idx]+delta,
            'color': color
        }
        
        if check_duplicate:
            for c in self.couplings.values():
                for bond in c['bond_list']:
                    ij = (bond['i'], bond['j'])
                    if  ij == (new_bond['i'], new_bond['j']):
                        # Adding another bond between existing sites, make sure we aren't doubling up!
                        # HACK: just use 1e-10, this will cause problem if someone tries to use Angstrom!
                        if linalg.norm(bond['rj - ri'] - new_bond['rj - ri']) < 1e-10:
                            print('WARN: Adding duplicate bond between sites `{}` and `{}`'.format(*ij))

                    if  ij == (new_bond['j'], new_bond['i']):
                        # Adding another bond between existing sites, make sure we aren't doubling up!
                        # HACK: just use 1e-10, this will cause problems if someone tries to use Angstrom!
                        if linalg.norm(bond['rj - ri'] + new_bond['rj - ri']) < 1e-10:
                            print('WARN: Adding duplicate bond between sites `{}` and `{}`'.format(*ij))

        self.couplings[coupling_name]['bond_list'].append(new_bond)

        
    def clear_couplings(self):
        for c in self.couplings:
            self.couplings[c]['bond_list'] = []
        
    # sets the rotator for site 'label'
    # should be a 3x3 matrix such that Rz = G
    def set_rotor(self, label, R):
#         if linalg.norm(np.dot(R, [0,0,1]) - self.sites[label]['G']) > 1e-6
        assert R.shape == (3,3)
        assert np.abs(linalg.det(R) - 1) < 1e-9
        self.sites[label]["rotor"] = R
    
    # sets the CGS vector of site 'label'
    # substitute for set_rotor
    def set_CGS(self, label, G):
        G = np.array(G,dtype=np.float64)
        assert G.shape == (3,)
        
        G = (1./linalg.norm(G))*G
        theta = np.arctan2(np.sqrt(G[0]**2+G[1]**2), G[2])
        phi = np.arctan2(G[1],G[0])
        
        self.sites[label]["rotor"] = np.matmul(Rz(phi),Ry(theta))
    
    def get_CGS(self, label):
        return self.sites[label]["rotor"][:,2]
    
        
        
        

def plot_ellipsoid(ax, coefs=(1,1,1),basis=np.eye(3),origin=np.zeros(3),**kwargs):
    #plagiarised from https://stackoverflow.com/questions/7819498/plotting-ellipsoid-with-matplotlib
#     Plots an ellipsoid of lengths set by 'coefs', along axes sepcified by 'basis'

    # Set of all spherical angles:
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)

    # Cartesian coordinates that correspond to the spherical angles:
    # (this is the equation of an ellipsoid):
    xyz = np.array( [ coefs[0]*np.outer(np.cos(u), np.sin(v)),
                      coefs[1]*np.outer(np.sin(u), np.sin(v)),
                      coefs[2]*np.outer(np.ones_like(u), np.cos(v)) ] )
    
    x,y,z = np.einsum('ij,jnm->inm',basis, xyz)
    
    # Plot:
    ax.plot_surface(x+origin[0], y+origin[1], z+origin[2], rstride=4, cstride=4, **kwargs)
        
        
class PlottingLattice(Lattice):
    '''
    Includes all the functions responsible for plotting real-space and k-space representations of the lattice
    
    '''
    
    def __init__(self, a1,a2,a1m=None,a2m=None):
        super().__init__(a1,a2,a1m,a2m)
        self.fig_KS = None
        self.ax_KS = None
        self.fig_RS = None
        self.ax_RS = None
        self.cellvert_dict = {
            "primitive": np.array([(1,1), (1, -1), (-1, -1), (-1, 1)])/2,
            "hexagon": np.array([(1,1), (-1, 2), (-2, 1), (-1, -1), (1, -2), (2, -1)])/3,
            "hexagon2": np.array([(2,1), (1, -1), (-1, -2), (-2, -1), (-1, 1), (1, 2)])/3,
        }
    
    
    # Plots a representation of the lattice composed of n1 * n2 unit cells
    def plot_atoms(self,n1=1,n2=1, arrowlen=0.5, plot_GS=True, plot_labels=True):
        self.fig_RS = plt.figure(figsize=(7,4))
        self.ax_RS = self.fig_RS.add_subplot(111, projection='3d')
            
        # Draw the unit cell
        cnrs = np.array([[0,0],self.a1,self.a1+self.a2,self.a2,[0,0]]) - 0.5*(self.a1+self.a2)
        self.ax_RS.plot(cnrs[:,0],cnrs[:,1],np.zeros_like(cnrs[:,0]),'k:')
        
        cnrs = np.array([[0,0],self.a1m,self.a1m+self.a2m,self.a2m,[0,0]]) - 0.5*(self.a1m+self.a2m)
        self.ax_RS.plot(cnrs[:,0],cnrs[:,1],np.zeros_like(cnrs[:,0]),'g--')
        
        self.ax_RS.set_zlim([-1,1])
        
        for i in range(-n1//2+1,n1//2+1):
            for j in range(-n2//2+1,n2//2+1):
                XY = i*self.a1m + j*self.a2m
                for label in self.sites:
                    s = self.sites[label]
                    xy = s['pos'] + XY
                    
                    # plot the location
                    if i==0 and j==0:
                        self.ax_RS.plot(xy[0],xy[1],0,'o',color='r')
                    else:
                        self.ax_RS.plot(xy[0],xy[1],0,'o',color='k')
                    
                    # plot the GS vector
                    if plot_GS==True:
                        G = arrowlen*self.spinspace_basis@np.matmul(s['rotor'],[0,0,1])
                        self.ax_RS.quiver(xy[0],xy[1],0, G[0],G[1],G[2])
                    
                    # add labels
                    if (plot_labels==True and i==0 and j==0) or plot_labels=='all':
                        self.ax_RS.text(xy[0],xy[1],0,label)
        #equalise the axes
        axisEqual3D(self.ax_RS)
        
    def plot_vector_at_atom(self, label, G, **kwargs):
        if self.ax_RS is None:
            self.plot_atoms(**kwargs)
            
        xy = self.sites[label]['pos']
        self.ax_RS.quiver(xy[0],xy[1],0, G[0],G[1],G[2])
    
    
    # possible styles:
    #     'net'  only the bonds
    #     'eigs'  plot an ellipsoid along the principal axes 
    def plot_couplings(self, n1=1,n2=1, bond_width=1,style='decompose',scale=1):
        
        self.plot_atoms(n1,n2, plot_GS=False,plot_labels=True)
        
        for c in self.couplings.values():
            for bond in c['bond_list']:
                if style == 'eigs':
                    l,v = linalg.eig(c['mat'])
                    Lmax = np.max(np.abs(l))
                    if Lmax > 1e-9:
                        l = l/Lmax
                    
                for i in range(-n1//2+1,n1//2+1):
                    for j in range(-n2//2+1,n2//2+1):
                        XY = i*self.a1m + j*self.a2m

                        delta = bond['rj - ri']
                        xy = self.sites[bond['i']]['pos'] + XY
                        # plot the bond itself
                        self.ax_RS.plot([xy[0],xy[0]+delta[0]], [xy[1],xy[1]+delta[1]],0,color=bond['color'],lw=bond_width)

                        if style=='eigs':
                            plot_ellipsoid(self.ax_RS, coefs=l, basis=scale*v, origin=np.append(xy+0.5*delta,0), color=bond['color'])
    
        axisEqual3D(self.ax_RS)
                        
                        
                        
                    
            
    # plots the Brillouin zone
    def plot_BZ(self, unitcell="primitive", plot_MBZ=True, plot_LBZ=True, *args, **kwargs):
    
        if type(unitcell) is list:
            cellverts = np.array(unitcell)
        elif type(unitcell) is np.ndarray:
            cellverts = unitcell
        elif type(unitcell) is str:
            try:
                cellverts = self.cellvert_dict[unitcell]
            except KeyError:
                print("unit cell specification not recognised, `"+unitcell+"`")
                return
        else:
            print("Invalid unit cell specification, `"+str(unitcell)+"`")
            return


        cellverts = np.matmul(cellverts, self.reciprocal_vectors)
        cellverts = np.append(cellverts,[cellverts[0]],axis=0)
        
        cellverts_m = np.matmul(self.cellvert_dict['primitive'], self.mag_reciprocal_vectors)
        cellverts_m = np.append(cellverts_m,[cellverts_m[0]],axis=0)

        self.fig_KS = plt.figure(figsize=(4,4))
        self.ax_KS = self.fig_KS.add_axes((0.1,0.1,0.8,0.9))
        
        if plot_LBZ:
            self.ax_KS.plot(*cellverts.T, ':',label='Lattice BZ', *args, **kwargs)
        if plot_MBZ:
            self.ax_KS.plot(*cellverts_m.T, label='Magnetic BZ', *args, **kwargs)
        self.ax_KS.set_aspect('equal', adjustable='box')
        
    # plots a path throught the BZ
    # path: a dict of points where each key is a label, 
    #     each value is a pair of coordinates in reciprocal lattice units
    #     (based on a1,a2, not the magnetic unit cell!)
    def plot_BZ_path(self, path, *args, **kwargs):
        k, t = parameterised_path(path, self.reciprocal_vectors)

        T= np.linspace(0,t[-1],1000)

        self.ax_KS.plot(*k(T).T, *args, **kwargs)
        for k in path:
            self.ax_KS.annotate(k, np.dot(np.array(path[k]),self.reciprocal_vectors),xytext=(5,-10),textcoords='offset pixels')
            self.ax_KS.plot(*np.dot(np.array(path[k]).T,self.reciprocal_vectors),'D',color='k',markersize=2)
            
        self.ax_KS.set_aspect('equal', adjustable='box')

        
# returns N random unit vectors with the shape
# (N, 3)
def rand_unitvec(N):
    s = np.random.uniform(low=-1,high=1,size=(N, 3))
    # normalise...
    return (s.T/linalg.norm(s,axis=1)).T
        
        
class SpinWave(PlottingLattice):
    def __init__(self, a1,a2,a1m=None,a2m=None):
        super().__init__(a1,a2,a1m,a2m)
        self.spinspace_basis = np.eye(3)
        self.quiet = False
        
    # provide the spin basis in the form
    # [[v1], [v2], [v3]]
    def set_spin_basis(self, O):
        O = np.array(O).T
        if O.shape != (3,3):
            print("Provide spins in the format [e1, e2, e3]")
        elif not np.allclose(O@O.T, np.eye(3)):
            print("Vectors must be orthonormal")
        elif np.abs(linalg.det(O) - 1) > 1e-10:
            print("Coordinates must be right handed")
        else:
            self.spinspace_basis = O
        
    # returns the Z object defined above in the format
    # Z[i, j]
    def get_Zmat(self, alpha, beta, U=np.eye(3)):
        # build it blockwise
        
        # top left
        z_zb = np.zeros((self.natoms, self.natoms),dtype=np.complex128)
        # top right
        z_z = np.zeros_like(z_zb)
        # bottom left
        zb_zb = np.zeros_like(z_zb)
        # bottom right
        zb_z = np.zeros_like(z_zb)
        
        for l1 in self.sites:
            idx1 = self.sites[l1]['index']
            z1 = (self.spinspace_basis@self.get_z(l1))[alpha]
            for l2 in self.sites:
                idx2 = self.sites[l2]['index']
                z2 = (self.spinspace_basis@self.get_z(l2))[beta]
                
                z_zb[idx1, idx2] = z1*np.conj(z2)
                z_z[idx1, idx2] = z1*z2
                zb_zb[idx1, idx2] = np.conj(z1*z2)
                zb_z[idx1, idx2] = np.conj(z1)*z2

        return np.block([[z_zb, z_z], [zb_zb, zb_z]])
        
        
    def get_z(self, label):
        R = self.sites[label]["rotor"]
        return R[:,0] + 1j*R[:,1]
        
        
    def calc_Hamiltonian(self, K):
        # accepts K in the shape K[N, 2]
        # N is arbitrary
        # calculates the two independents bits, A and B corresponding to the left 2
        # natoms*natoms blocks
        num_k = K.shape[0] if K.ndim == 2 else 1
        
        A = np.zeros((num_k, self.natoms, self.natoms),dtype=np.complex128)
        B = np.zeros_like(A)
        
        for c in self.couplings.values():
            J = c['mat']
            for bond in c['bond_list']:
                zi = self.get_z(bond['i'])
                zj = self.get_z(bond['j'])

                xi = self.get_CGS(bond['i'])
                xj = self.get_CGS(bond['j'])

                si = self.sites[bond['i']]
                sj = self.sites[bond['j']]


                delta = bond['rj - ri']

                phase = np.exp(-1j*np.dot(K, delta))

                x = 0.5*np.einsum('a,ab,b',np.conj(zi),J,zj)*phase

                A[:, si['index'],sj['index']] += x
                A[:, sj['index'],si['index']] += np.conj(x)

                B[:, si['index'],sj['index']] += 0.5*np.einsum('a,ab,b',zi,J,zj)*phase
                B[:, sj['index'],si['index']] += 0.5*np.einsum('a,ab,b',zj,J,zi)*np.conj(phase)

                A[:, si['index'],si['index']] -= np.einsum('a,ab,b',xi,J,xj)
                A[:, sj['index'],sj['index']] -= np.einsum('a,ab,b',xj,J,xi)
            
        return np.concatenate((np.concatenate((A,B),axis=1),np.concatenate((B.conj().transpose(0,2,1), A),axis=1)),axis=2)
    
    
    @property
    def first_order_terms(self):
        residuals = np.zeros((self.natoms),dtype=complex)
        for c in self.couplings.values():
            J = c['mat']
            for bond in c['bond_list']:
                residuals[self.sites[bond['i']]['index']] += np.dot(self.get_z(bond['i']),np.matmul(J,self.get_CGS(bond['j'])))
                residuals[self.sites[bond['j']]['index']] += np.dot(self.get_z(bond['j']),np.matmul(J,self.get_CGS(bond['i'])))
        return residuals
    
    def energy(self, K, tol=1e-9):
        # accepts K in the shape K[N, 2]
        # K is in standard dimensionless Cartesian units
        G = np.diag(np.repeat([1,-1],self.natoms))
        E = linalg.eig(np.einsum('ir,nrj -> nij', G, self.calc_Hamiltonian(K)))[0]
        
        failed = linalg.norm(np.imag(E),axis=1)>tol
        
        if True in failed:
            print("Warning! Diagonalisation was not successful at {}% of K values".format(100*len(E[failed,0])/len(failed)))
            print(K[failed, :])
        return np.sort(np.real(E),axis=1)
    
    def isgapped(self, K, mincond=100):
        G = np.diag(np.repeat([1,-1],self.natoms))
        return linalg.cond(np.einsum('ir,nrj -> nij', G, self.calc_Hamiltonian(K))) < mincond
    
    def isposdef(self, K):
        return np.all(linalg.eigh(self.calc_Hamiltonian(K))[0]>0, axis=1)
        
    

    def Colpa_diagonalise(self, K):
        # Reference: Appendix A of Smit et al., 'Magnon daming in the zigzag phase of the Kitaev-Heisenberg-Γ model on a honeycomb lattice' 2020
        # Algorithm originally due to Colpa, 1978
        G = np.array([1,-1],dtype=np.complex128)
        G = np.diag(np.repeat(G,self.natoms))
#         E, Bp = linalg.eig(np.einsum('ir,nrj -> nij', G, self.calc_Hamiltonian(K)))

        H = self.calc_Hamiltonian(K) + 1e-8*np.eye(self.natoms*2)
        L = np.zeros_like(H)
#         for n in range(H.shape[0]):
        try:
            L = linalg.cholesky(H)
        except linalg.LinAlgError:
            if not self.quiet:
                print("Matrix not positive definite somewhere !")
            E = np.zeros((K.shape[0], self.natoms*2),np.float64)
            E[:] = np.nan
            B = np.zeros((K.shape[0], self.natoms*2,self.natoms*2),np.complex128)
            return (E,B)
    
        LH = L.transpose(0,2,1).conj()
    
        LhGL = np.einsum('nip, pq, nqj -> nij', LH, G, L)

        E, U = linalg.eigh(LhGL)
        
        # check for positive definiteness
        # KLUDGE: just use 1e-9 arbitrarily
        failed = linalg.norm(np.imag(E),axis=1)>1e-9
        if True in failed:
            print("Warning! Diagonalisation was not successful at {}% of K values".format(100*len(E[failed,0])/len(failed)))
            print(K[failed, :])

            
        # sort E and U in the order 
        # 1, 3.2, 5.8, ... -1, -2.2, -5.8, ...
        
        idx = np.argsort(np.real(E),axis=1)
        idx = np.hstack((idx[:, self.natoms:], np.flip(idx[:, :self.natoms],axis=1)))
        
        for n in range(E.shape[0]):
            E[n,:] = E[n, idx[n]]
            U[n, :] = U[n][:, idx[n]]
        
#         assert np.max(linalg.norm(np.einsum('nir,nrj -> nij',U.conj().transpose(0,2,1), U) - np.eye(self.natoms*2), axis=(1,2))) < 1e-7
        
        E = np.real(E)
        
        # Solve LH B = U (GE)^1/2
        # noting that GE is diagonal
        # XXX Note: needs to check that E has the claimed positivity structure
        B = linalg.solve(LH, np.einsum('nij, nj->nij',U, np.sqrt(np.abs(E))))
        
        
        return (E, B)
    
    def Bogo_diagonalise(self, K):
        G = np.array([1,-1],dtype=np.complex128)
        G = np.diag(np.repeat(G,self.natoms))
        E, B = linalg.eig(np.einsum('ij,njk->nik',G, self.calc_Hamiltonian(K)))
        
        # reorder the eigenvectors
        idx = np.argsort(np.real(E),axis=1)
        idx = np.hstack((idx[:, self.natoms:], np.flip(idx[:, :self.natoms],axis=1)))
        
        norms = np.zeros(self.natoms*2)
        
        for n in range(E.shape[0]):
            E[n,:] = E[n, idx[n]]
            b = B[n][:, idx[n]]
            norms = np.power(np.abs(np.diag(b.conj().T @ G @ b )),-0.5)
            B[n, :, :] = norms*b
        
        
        return (E, B)
    
    
    def get_ffactor(self, K, alpha, beta):
        Z = self.get_Zmat(alpha, beta)
        E, B = self.Colpa_diagonalise(K)
        
        return ( E, np.einsum('nip,pq,nqi->ni' , B.conj().transpose(0,2,1), Z, B) )
    

    def get_perp_ffactor(self, K):
        if K.shape[1] == 2:
            print(K.shape)
            K = np.c_[K, np.ones(K.shape[0])]
            
        
        tensorQ =  - (np.einsum('na,nb->nab',K, K).T/np.einsum('na,na->n',K, K)).T
        tensorQ += np.full_like(tensorQ,np.eye(3))
        Z = np.multiply.outer(tensorQ[:,0,0], self.get_Zmat(0,0)) + np.multiply.outer(tensorQ[:,1,1], self.get_Zmat(1,1)) + np.multiply.outer(tensorQ[:,2,2], self.get_Zmat(2,2))
        Z += np.multiply.outer(tensorQ[:,0,1], self.get_Zmat(0,1)) + np.multiply.outer(tensorQ[:,1,0], self.get_Zmat(1,0))
        Z += np.multiply.outer(tensorQ[:,0,2], self.get_Zmat(0,2)) + np.multiply.outer(tensorQ[:,2,0], self.get_Zmat(2,0))
        Z += np.multiply.outer(tensorQ[:,2,1], self.get_Zmat(2,1)) + np.multiply.outer(tensorQ[:,1,2], self.get_Zmat(1,2))

#         Z = self.get_Zmat(0,0) + self.get_Zmat(1,1)

        E, B = self.Colpa_diagonalise(K[:, :2])
        
        return ( E, np.einsum('nip,npq,nqi->ni' , B.conj().transpose(0,2,1), Z, B) )
    
    @property
    def spinham(self):
        H = np.zeros((self.natoms*3,self.natoms*3),dtype=np.float64)
        for c in self.couplings.values():
            J = c['mat']
            for bond in c['bond_list']:
                si = self.sites[bond['i']]
                sj = self.sites[bond['j']]
                i = 3*si['index']
                j = 3*sj['index']
                H[i:i+3,j:j+3] += J
                H[j:j+3,i:i+3] += J
        return H
        
    
    
    ################################################################
    ###### Plotting Functions    
    
    def plot_bandgraph(self, path, *args, **kwargs):
        # Plots the dispersion relation omega along the given path (reciprocal lattice units)
        # Based on the self-computed/ self-constructed Hamiltonian
        
        if np.max(np.abs(self.first_order_terms)) > 1e-10:
            print("WARNING: first order terms do not vanish")
            print("Suspect that the classical ground state is incorrect")
            print(self.first_order_terms)
        
        self.plot_bandgraph_custom(path, self.energy, *args, **kwargs)
        
    def plot_bandgraph_custom(self, path, omega, *args, **kwargs):
        # Plots the supplied dispersion relation omega along the given path
        B = self.reciprocal_vectors
        
        k, t = parameterised_path(path, B)
        T= np.linspace(0,t[-1],1000)

        K = k(T)

        fig, ax = plt.subplots(figsize=(7,5))

        E = omega(K)
        
        ax.plot(T, E, *args, **kwargs)
        ax.set_ylabel('Energy')
        ax.set_ylim([0,np.max(E)*1.15])
        ax.set_xticks(t)
        ax.set_xticklabels(labels = list(path.keys()))
        ax.xaxis.grid()
        
    def plot_band_2D(self, idx='min', npoints = (50,50), color_kwargs={}):
        
        if type(idx) is int:
            assert idx >= 0
            assert idx < self.natoms*2
        elif type(idx) is str:
            assert idx in ['min', 'max']
        else:
            raise Exception('unrecognised idx directive {} (must be number between 0 and {}, or else "min" or "max")'.format(idx, self.natoms*2))
                
        
        self.plot_BZ()
        xmin, xmax = self.ax_KS.get_xlim()
        ymin, ymax = self.ax_KS.get_ylim()
        
        X = np.linspace(xmin, xmax, npoints[0])
        Y = np.linspace(ymin, ymax, npoints[1])
        Z = np.zeros(npoints)
        
        one = np.ones_like(X)
        
        for iy in range(npoints[1]):
            K = np.array([X, one*Y[iy]]).T
            
            E = self.energy(K)
            if idx == 'min':
                Z[iy] = np.min(E,axis=1)
            elif idx == 'max':
                Z[iy] = np.max(E,axis=1)
            else:
                Z[iy] = E[:,idx]
                
        # Make the units right....
        X /= self.RS_units[0]
        Y /= self.RS_units[0]
    
        c = self.ax_KS.pcolormesh(X,Y,Z,shading='nearest',**color_kwargs)
        self.fig_KS.colorbar(c)
        
from scipy.special import voigt_profile
from scipy.interpolate import interp1d
       
    
# look, this is a bot of a mess ngl
class PowderAverage(object):
    def __init__(self, spinw, gauss_broaden =0.135, intrinsic_broaden=0):
        self.spinw = spinw
        self.ax_cuts = None
        self.fig_cuts = None
        self.S=None
        self.Egrid=None
        self.Qgrid=None
        self.gauss_broaden = gauss_broaden
        self.intrinsic_broaden = intrinsic_broaden
        self.form_factor = None
        
    def set_form_factor(self, func):
        self.form_factor = func
        
    def set_form_factor_arr(self, Q, f):
        self.form_factor = interp1d(Q, f, kind='cubic')
        
    # Don't ask why this function is in a class called PowderAverage, okay?
    def plot_single_crystal(self, path, Egrid, ab_index, use_ffactor=False, plot_dispersion=True, lw=0.3):
        
        ab_index = ab_index.lower()
        translate = {
            'x': 0,
            'y': 1,
            'z': 2
        }
        
        # Convert the pathspec into a numerical path 
        B = self.spinw.reciprocal_vectors
        k, t = parameterised_path(path, B)
        T= np.linspace(0,t[-1],1000)
        K = k(T)


        # Figure out what S actually looks like
        if ab_index == 'perp':
            E,S = self.spinw.get_perp_ffactor(K)
            ab_index=r'\perp'
        elif ab_index in 'xx xy xz yx yy yz zx zy zz'.split():
            alpha = translate[list(ab_index)[0]]
            beta = translate[list(ab_index)[1]]
            E,S = self.spinw.get_ffactor(K, alpha, beta)
        else:
            raise Exception('Unrecognised index specification'+ab_index)
        
        E = np.real(E)
        S = np.real(S)
    
        # take only the positive energies (rest are redundant)
        permute = np.argsort(E,axis=1)
        E = np.take_along_axis(E, permute, axis=1)[:, self.natoms:]
        S = np.take_along_axis(S, permute, axis=1)[:, self.natoms:]
        
        onevec = np.ones_like(T)
        
        Emesh = np.outer(onevec, Egrid)
        Spec = np.zeros_like(Emesh)
        
        for i in range(self.natoms):
            Spec += voigt_profile(Emesh-np.outer(E[:,i],np.ones(num_E)),self.gauss_broaden, self.intrinsic_broaden) *np.outer(S[:,i],np.ones(num_E))

            
        fig, ax = plt.subplots(figsize=(7,5))
        
        # Deals with the magnetic form factors
        if use_ffactor:
            if self.form_factor is not None:
                F2 = self.form_factor(self.Qgrid)**2
            else:
                raise RuntimeError('No form factor data is known!')
        else:
            F2 = np.ones_like(self.Qgrid)
        
        CB = ax.pcolormesh(T,Emesh,F2*Spec.T,shading='nearest')
        
        if plot_dispersion:
            ax.plot(T, E, color='white',lw=lw)
        
        ax.set_ylabel('Energy')
        ax.set_xticks(t)
        ax.set_xticklabels(labels = list(path.keys()))
        ax.xaxis.grid()
        ax.set_title(r'$ S^{'+str(ab_index)+r'}(\mathbf{Q},\omega)$')
        
        fig.colorbar(CB)
        
    def calc_paverage(self, Egrid, Qgrid=None, num_samples=500):
        
        if Qgrid is None:
            Qgrid = np.linspace(0, np.max(linalg.norm(self.spinw.reciprocal_vectors,axis=1)),50)
        else:
            Qgrid = np.sort(Qgrid)
        
        
        
        # Uniformly choose random points on the unit sphere
        

        S = np.zeros((Qgrid.shape[0], Egrid.shape[0]))
        
        for j in range(Qgrid.shape[0]):
            S[j] = self.sphere_average(Qgrid[j], Egrid, num_samples)
        
        self.S = S
        self.Qgrid = Qgrid 
        self.Egrid = Egrid
        
        
        
    def plot_paverage(self, logscale=False, vmin=None, vmax=None, cmap='jet', use_ffactor=False, ax=None, **kwargs):
        if self.S is None:
            raise RuntimeError('No powder average data found, run PowderAverage.calc_paverage first!')
        
        if use_ffactor:
            if self.form_factor is not None:
                F2 = self.form_factor(self.Qgrid)**2
            else:
                raise RuntimeError('No form factor data is known!')
        else:
            F2 = np.ones_like(self.Qgrid)

            
        make_fig = False # flag for whether we're handling the figure creation or the user
        if ax is None:
            make_fig = True
            fig, ax = plt.subplots()
            
        if logscale:
            c = ax.pcolormesh(self.Qgrid, self.Egrid, F2*self.S.T+1e-10, shading='nearest', cmap=cmap, vmin=vmin, vmax=vmax, norm=matplotlib.colors.LogNorm())
        else:
            c = ax.pcolormesh(self.Qgrid, self.Egrid, F2*self.S.T, shading='nearest', cmap=cmap, vmin=vmin, vmax=vmax)
        
        if make_fig:
            fig.colorbar(c)
            ax.set_xlabel("Momentum transfer $Q$, %s$^{-1}$" % self.spinw.RS_units[1])
            ax.set_ylabel("Energy transfer, meV")
        return c
        

    
    def sphere_average(self, Q, Egrid, num_samples, dQ=0):
            
        if dQ > 0:
            Q = np.random.uniform(low=(Q-dQ)**3, high=(Q+dQ)**3, size=num_samples)**(1/3)
            
        K = (self.spinw.RS_units[0]*rand_unitvec(num_samples).T*Q).T
        ee, s = self.spinw.get_perp_ffactor(K)
        one_ee = np.ones_like(ee[:,0])
        S = np.zeros_like(Egrid)

        
        for i in range(self.spinw.natoms):
            S += np.real(np.dot(s[:, i],voigt_profile(np.outer(one_ee,Egrid) - np.outer(ee[:,i],np.ones_like(Egrid)), self.gauss_broaden, self.intrinsic_broaden)))
       
        return S/num_samples
                
            
    def plot_cut(self, Q, Egrid, npoints=1000, dQ=0, use_ffactor=True, ax=None, **kwargs):
        if ax is None:
            if self.ax_cuts is None:
                self.fig_cuts, self.ax_cuts = plt.subplots()
            ax = self.ax_cuts
            
        if use_ffactor:
            if self.form_factor is not None:
                F2 = self.form_factor(Q)**2
            else:
                raise RuntimeError('No form factor data is known!')
        else:
            F2 = 1
            
        S = self.sphere_average(Q, Egrid, npoints, dQ)

        ax.plot(Egrid, F2*S, label='Q=[%.2f - %.2f] %s$^{-1}$' % (Q-dQ,Q+dQ, self.spinw.RS_units[1]), **kwargs)
        ax.set_xlabel('Energy')
        ax.set_ylabel('Intensity')
        
        
        
        
        
        
######## Classical ground state finding
        
def gradient_optimise(spinw, num_guesses = 100, cut=20, gamma = 1e-2, rtol = 1e-6, vtol=1e-5, niter = 500, invert_atom=None, invert_vec=[0,0,1]):
    M = spinw.spinham

    ######## Use MC scattershot to get some convincing guesses

    # initialise to randomness
    S = []
    for i in range(spinw.natoms):
        S.append(rand_unitvec(num_guesses))
    S = np.block(S)

    # calculate the energy of the current guesses
    energy = np.einsum('nj,jk,nk->n',S,M,S)

    ##### Alternate between gradient descent and renormalisation to get spins back to the unit spehere. #####

    result = {'success': False,
             'iterations': niter,
              'message': ''
             }

    counter = 0
    for n in range(niter):
        counter += 1
        # descend based on the Jacobian
        newS =  S - gamma*2*np.dot(S, M)

        # renormalise
        for i in range(0,len(S),3):
            newS[:, i:i+3] = (newS[:, i:i+3].T/np.linalg.norm(newS[:,i:i+3],axis=1)).T

        # exit if we have converged
        min_idx = np.argmin(np.einsum('nj,jk,nk->n',newS,M,newS))

        if (np.linalg.norm(newS[min_idx]-S[min_idx]) < rtol*np.linalg.norm(newS[min_idx])):
            S = newS
            result['success']=True
            result['iterations']=n
            break
        S = newS



    final_energies = np.einsum('nj,jk,nk->n',S,M,S)
    final_spins = S.reshape(num_guesses, spinw.natoms, 3)

    idx = np.argpartition(final_energies, cut)

    final_spins = final_spins[idx[:cut]]
    final_energies = final_energies[idx[:cut]]

    if invert_vec is not None and invert_atom is not None:
        # rotate optimised spins
        if type(invert_vec) is not int:
            invert_vec = [0,0,1]
        for n in range(num_guesses):
            if final_spins[n, invert_atom].dot(invert_vec) < 0:
                final_spins[n] = -final_spins[n]

    result['energy_mean'] = np.var(final_energies)
    result['energy_var'] = np.var(final_energies)
    result['spin_mean'] = np.mean(final_spins,axis=0)
    result['spin_var'] = np.var(final_spins,axis=0)


    good_energy = result['energy_var'] < vtol**2
    good_spins =  np.max(result['spin_var'])< vtol

    if good_energy:
        result['success'] = True
        if not good_spins:
            result['message'] = "Energy converged, but spins did not\npossible continuous family of degenerate solutions?"
    else:
        result['success'] = False
        if not good_spins:
            result['message'] = "Energy did not converge to within tolerance %f, possible phase competition?" % vtol**2
        else:
            result['message'] = "Method did not converge to within tolerance %f." % vtol

    result['values'] = (final_spins, final_energies)
    idx = np.argmin(final_energies)
    return (final_spins[idx], final_energies[idx], result)