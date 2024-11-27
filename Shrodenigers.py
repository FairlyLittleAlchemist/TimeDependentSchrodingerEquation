    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation
    from matplotlib.patches import Rectangle
    from scipy.sparse import csr_matrix
    from scipy.sparse.linalg import norm



    import time

    def conjugate_gradient_sparse(A, b, x0=None, tol=1e-1, max_iter=100):     
        if not isinstance(A, csr_matrix):
            A = csr_matrix(A)
        
        
        n = A.shape[0]
        if x0 is None:
            x0 = np.zeros(n)
        
        
        r = b - A.dot(x0)
        p = r.copy()
        
        
        r_norm = np.linalg.norm(r)
        if r_norm < tol:
            return x0
        
        
        x = x0
        for k in range(max_iter):
            
            Ap = A.dot(p)
            
            
            rTr = r.dot(r)
            alpha = rTr / p.dot(Ap)
            
            
            x = x + alpha * p
            
            
            r_new = r - alpha * Ap
            r_new_norm = np.linalg.norm(r_new)  
            
            
            if r_new_norm < tol:
                return x
            
            
            beta = r_new.dot(r_new) / rTr
            
            
            p = r_new + beta * p
            
            
            r = r_new
        
        return x
        


    def is_positive_definite(A):
        try:
            np.linalg.cholesky(A)
            return True
        except np.linalg.LinAlgError:
            return False



    def psi0(x, y, x0, y0, sigma=0.5, k=15*np.pi):
        
        """
        Proposed wave function for the initial time t=0.
        Initial position: (x0, y0)
        Default parameters:
            - sigma = 0.5 -> Gaussian dispersion.
            - k = 15*np.pi -> Proportional to the momentum.
            
        Note: if Dy=0.1 use np.exp(-1j*k*(x-x0)), if Dy=0.05 use 
            np.exp(1j*k*(x-x0)) so that the particle will move 
            to the right.
        """
        
        return np.exp(-1/2*((x-x0)**2 + (y-y0)**2)/sigma**2)*np.exp(1j*k*(x-x0))






    L = 7
    Dy = 0.05 
    Dt = Dy**2/4 
    Nx = int(L/Dy) + 1 
    Ny = int(L/Dy) + 1 
    Nt = 500 
    rx = -Dt/(2j*Dy**2) 
    ry = -Dt/(2j*Dy**2) 


    x0 = L/5
    y0 = L/2


    w = 0.6 
    s = 0.8 
    a = 0.2 



    j0 = int(1/(2*Dy)*(L-w)) 
    j1 = int(1/(2*Dy)*(L+w)) 

    i0 = int(1/(2*Dy)*(L+s) + a/Dy) 
    i1 = int(1/(2*Dy)*(L+s))        
    i2 = int(1/(2*Dy)*(L-s))        
    i3 = int(1/(2*Dy)*(L-s) - a/Dy) 


    v0 = 200
    v = np.zeros((Ny,Ny), complex) 
    v[0:i3, j0:j1] = v0
    v[i2:i1,j0:j1] = v0
    v[i0:,  j0:j1] = v0
        
    Ni = (Nx-2)*(Ny-2)  

    A = np.zeros((Ni,Ni), complex)
    M = np.zeros((Ni,Ni), complex)


    for k in range(Ni):     
        
        
        i = 1 + k//(Ny-2)
        j = 1 + k%(Ny-2)
        
        
        A[k,k] = 1 + 2*rx + 2*ry + 1j*Dt/2*v[i,j]
        M[k,k] = 1 - 2*rx - 2*ry - 1j*Dt/2*v[i,j]
        
        if i != 1: 
            A[k,(i-2)*(Ny-2)+j-1] = -ry 
            M[k,(i-2)*(Ny-2)+j-1] = ry
            
        if i != Nx-2: 
            A[k,i*(Ny-2)+j-1] = -ry
            M[k,i*(Ny-2)+j-1] = ry
        
        if j != 1: 
            A[k,k-1] = -rx 
            M[k,k-1] = rx 

        if j != Ny-2: 
            A[k,k+1] = -rx
            M[k,k+1] = rx


    print(A.shape)

    Asp = csr_matrix(A)  
    x = np.linspace(0, L, Ny-2)  
    y = np.linspace(0, L, Ny-2)
    x, y = np.meshgrid(x, y)
    psis = []  


    psi = psi0(x, y, x0, y0)
    psi[0,:] = psi[-1,:] = psi[:,0] = psi[:,-1] = 0  
    psis.append(np.copy(psi))  


    psi = psi0(x, y, x0, y0)
    psi_vect = psi.reshape((Ni))
    b = np.matmul(M, psi_vect)  


    for i in range(1,Nt):
        psi_vect = psi.reshape((Ni)) 
        b = np.matmul(M,psi_vect) 
        psi_vect =conjugate_gradient_sparse(Asp, b, x0=psi_vect) 
        psi = psi_vect.reshape((Nx-2,Ny-2)) 
        psis.append(np.copy(psi)) 



    mod_psis = [] 
    for wavefunc in psis:
        re = np.real(wavefunc) 
        im = np.imag(wavefunc) 
        mod = np.sqrt(re**2 + im**2) 
        mod_psis.append(mod) 























    fig = plt.figure() 
    ax = fig.add_subplot(111, xlim=(0,L), ylim=(0,L)) 

    img = ax.imshow(mod_psis[0], extent=[0,L,0,L], cmap=plt.get_cmap("hot"), vmin=0, vmax=np.max(mod_psis), zorder=1) 


    slitcolor = "w" 
    slitalpha = 0.08 
    wall_bottom = Rectangle((j0*Dy,0),     w, i3*Dy,      color=slitcolor, zorder=50, alpha=slitalpha) 
    wall_middle = Rectangle((j0*Dy,i2*Dy), w, (i1-i2)*Dy, color=slitcolor, zorder=50, alpha=slitalpha)
    wall_top    = Rectangle((j0*Dy,i0*Dy), w, i3*Dy,      color=slitcolor, zorder=50, alpha=slitalpha)


    ax.add_patch(wall_bottom)
    ax.add_patch(wall_middle)
    ax.add_patch(wall_top)



    def animate(i):
        
        """
        Animation function. Paints each frame. Function for Matplotlib's 
        FuncAnimation.
        """
        
        img.set_data(mod_psis[i]) 
        img.set_zorder(1)
        
        return img, 


    anim = FuncAnimation(fig, animate, interval=1, frames =np.arange(0,Nt,2), repeat=False, blit=0) 
    plt.show() 