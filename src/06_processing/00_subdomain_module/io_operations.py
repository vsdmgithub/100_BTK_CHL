import numpy as np
import concurrent.futures
from tqdm import tqdm

def write_ucvcwc_file(U,V,W,num=0,addr='./'):
    """
    Writes ucvcwc restart file at a given timestep.

    :param addr:   name of the folder to write to
    :param [U,V,W]: velocity components
    :param num:     number of the file to read (default 0)
    """
    # Note: main code expects z,x,y ordering

    fname = addr+'01uc00.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(U)))
    fout.close()

    fname = addr + '02vc00.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(V)))
    fout.close()

    fname = addr + '03wc00.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(W)))
    fout.close()

def write_boundary_file(U,V,W,num=0,addr='./',flag=-1):
    """
    Writes boundary files for subdomain calculation. 

    :param [U,V,W]: velocity components
    :param num:     timestep
    :param addr:    name of the folder to write to
    :param flag:    flag indicating which boundary to write:
                    -1, write all boundaries
                    1, bot boundary
                    2, top boundary
                    3, left boundary
                    4, right boundary
                    5, back boundary
                    6, front boundary
    :return: None
    """
    if flag == -1 or flag == 1:
        #bot boundary data
        Ubc = U[:,:,0] 
        Vbc = V[:,:,1] 
        Wbc = W[:,:,0]
        fname = addr+'bc_bot.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc)))
        fout.write(bytes(np.transpose(Vbc)))
        fout.write(bytes(np.transpose(Wbc)))
        fout.close()

    if flag == -1 or flag == 2:
        #top boundary data
        Ubc = U[:,:,-1] 
        Vbc = V[:,:,-2] 
        Wbc = W[:,:,-1] 
        fname = addr+'bc_top.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc)))
        fout.write(bytes(np.transpose(Vbc)))
        fout.write(bytes(np.transpose(Wbc)))
        fout.close()

    if flag == -1 or flag == 3:
        #left boundary data
        Ubc = U[:,1,:]
        Vbc = V[:,0,:] 
        Wbc = W[:,0,:] 
        fname = addr+'bc_lft.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc)))
        fout.write(bytes(np.transpose(Vbc)))
        fout.write(bytes(np.transpose(Wbc)))
        fout.close()

    if flag == -1 or flag == 4:
        #rht boundary data
        Ubc = U[:,-2,:]
        Vbc = V[:,-1,:]
        Wbc = W[:,-1,:]
        fname = addr+'bc_rht.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc)))
        fout.write(bytes(np.transpose(Vbc)))
        fout.write(bytes(np.transpose(Wbc)))
        fout.close()

    if flag == -1 or flag == 5:
        #back boundary data
        Ubc = U[0,:,:] 
        Vbc = V[0,:,:]
        Wbc = W[1,:,:] 
        fname = addr+'bc_bck.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc)))
        fout.write(bytes(np.transpose(Vbc)))
        fout.write(bytes(np.transpose(Wbc)))
        fout.close()

    if flag == -1 or flag == 6:
        #front boundary data
        Ubc = U[-1,:,:] 
        Vbc = V[-1,:,:]
        Wbc = W[-2,:,:] 
        fname = addr+'bc_frt.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc)))
        fout.write(bytes(np.transpose(Vbc)))
        fout.write(bytes(np.transpose(Wbc)))
        fout.close()

def write_boundary_file_from_six_boundaries(Ulft,Vlft,Wlft,\
                               Urht,Vrht,Wrht,\
                               Utop,Vtop,Wtop,\
                               Ubot,Vbot,Wbot,\
                               Ubck,Vbck,Wbck,\
                               Ufrt,Vfrt,Wfrt,\
                               num=0,bc_addr='./',):
    # * left boundary data
    """
    Writes boundary files for the given boundary conditions.

    Parameters:
    - Ulft, Vlft, Wlft: Left boundary velocity components
    - Urht, Vrht, Wrht: Right boundary velocity components
    - Ubot, Vbot, Wbot: Bottom boundary velocity components
    - Utop, Vtop, Wtop: Top boundary velocity components
    - Ubck, Vbck, Wbck: Back boundary velocity components
    - Ufrt, Vfrt, Wfrt: Front boundary velocity components
    - num: The timestep number (default is 0)
    - bc_addr: The address of the boundary file directory (default is './')
    """
    fname = bc_addr+'bc_lft.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Ulft[:,0,:])))
    fout.write(bytes(np.transpose(Vlft[:,0,:])))
    fout.write(bytes(np.transpose(Wlft[:,0,:])))
    fout.close()

    # * rht boundary data
    fname = bc_addr+'bc_rht.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Urht[:,0,:])))
    fout.write(bytes(np.transpose(Vrht[:,0,:])))
    fout.write(bytes(np.transpose(Wrht[:,0,:])))
    fout.close()

    # * bot boundary data
    fname = bc_addr+'bc_bot.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Ubot[:,:,0])))
    fout.write(bytes(np.transpose(Vbot[:,:,0])))
    fout.write(bytes(np.transpose(Wbot[:,:,0])))
    fout.close()

    # * top boundary data
    fname = bc_addr+'bc_top.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Utop[:,:,0])))
    fout.write(bytes(np.transpose(Vtop[:,:,0])))
    fout.write(bytes(np.transpose(Wtop[:,:,0])))
    fout.close()

    # * back boundary data
    fname = bc_addr+'bc_bck.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Ubck[0,:,:])))
    fout.write(bytes(np.transpose(Vbck[0,:,:])))
    fout.write(bytes(np.transpose(Wbck[0,:,:])))
    fout.close()

    # * front boundary data
    fname = bc_addr+'bc_frt.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Ufrt[0,:,:])))
    fout.write(bytes(np.transpose(Vfrt[0,:,:])))
    fout.write(bytes(np.transpose(Wfrt[0,:,:])))
    fout.close()

def write_boundary_file_from_boundary(Ubc,Vbc,Wbc,num=0,bc_addr='./',flag=-1, adjFlag=False):
    """
    Writes boundary files for the given boundary conditions.

    Parameters:
    - Ubc, Vbc, Wbc: velocity components of the boundary
    - num: The timestep number (default is 0)
    - bc_addr: The address of the boundary file directory (default is './')
    - flag: The flag indicating which boundary to write:
            1, bot boundary
            2, top boundary
            3, left boundary
            4, right boundary
            5, back boundary
            6, front boundary
    """
    if flag == 1:
        # * bot boundary data
        fname = bc_addr+'bc_bot.dble.'+"{:07d}".format(num)
        if adjFlag:
            fname = bc_addr+'adjbot.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc[:,:])))
        fout.write(bytes(np.transpose(Vbc[:,:])))
        fout.write(bytes(np.transpose(Wbc[:,:])))
        fout.close()

    elif flag == 2:
        # * top boundary data
        fname = bc_addr+'bc_top.dble.'+"{:07d}".format(num)
        if adjFlag:
            fname = bc_addr+'adjtop.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc[:,:])))
        fout.write(bytes(np.transpose(Vbc[:,:])))
        fout.write(bytes(np.transpose(Wbc[:,:])))
        fout.close()
    
    elif flag == 3:
        # * left boundary data
        fname = bc_addr+'bc_lft.dble.'+"{:07d}".format(num)
        if adjFlag:
            fname = bc_addr+'adjlft.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc[:,:])))
        fout.write(bytes(np.transpose(Vbc[:,:])))
        fout.write(bytes(np.transpose(Wbc[:,:])))
        fout.close()

    elif flag == 4:
        # * rht boundary data
        fname = bc_addr+'bc_rht.dble.'+"{:07d}".format(num)
        if adjFlag:
            fname = bc_addr+'adjrht.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc[:,:])))
        fout.write(bytes(np.transpose(Vbc[:,:])))
        fout.write(bytes(np.transpose(Wbc[:,:])))
        fout.close()
        fout.close()

    elif flag == 5:
        # * back boundary data
        fname = bc_addr+'bc_bck.dble.'+"{:07d}".format(num)
        if adjFlag:
            fname = bc_addr+'adjbck.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc[:,:])))
        fout.write(bytes(np.transpose(Vbc[:,:])))
        fout.write(bytes(np.transpose(Wbc[:,:])))
        fout.close()

    elif flag == 6:
        # * front boundary data
        fname = bc_addr+'bc_frt.dble.'+"{:07d}".format(num)
        if adjFlag:
            fname = bc_addr+'adjfrt.dble.'+"{:07d}".format(num)
        fout = open(fname,'wb')
        fout.write(bytes(np.transpose(Ubc[:,:])))
        fout.write(bytes(np.transpose(Vbc[:,:])))
        fout.write(bytes(np.transpose(Wbc[:,:])))
        fout.close()
    
    else:
        print("Invalid flag value",flush=True)

def read_ucvcwc_file(fileroot,ngx,ngy,ngz,num,adjointFlag=0,component=-1):
    """
    Reads U,V,W velocity data from binary fortran files 

    :param fileroot:  Location of data files
    :param ng[x,y,z]: number of grid points in each dimension
    :param num:       Time step to extract data
    :param adjointFlag: [1] read adjoint field [else] read forward field
    :param component:   [-1] all, [0] U, [1] V, [2] W


    :return: The velocity field components as numpy arrays 
        
            - U: The streamwise velocity field components 
            - V: The vertical velocity field components 
            - W: The spanwise velocity field components 
    """

    fname =  '01uc00.dble.'
    if (adjointFlag == 1):
        fname =  '01ucad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Ufield =  np.fromfile(filename,dtype=float,count=(ngx+2)*(ngz+1)*(ngy+1),offset=0)
    Ufield =  np.transpose(np.reshape(Ufield,(ngy+1,ngx+2,ngz+1))) 
    if component == 0:
        return num,Ufield

    fname = '02vc00.dble.'
    if (adjointFlag == 1):
        fname = '02vcad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Vfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(ngy+2),offset=0)
    Vfield =  np.transpose(np.reshape(Vfield,(ngy+2,ngx+1,ngz+1))) 
    if component == 1:
        return num,Vfield

    fname = '03wc00.dble.'
    if (adjointFlag == 1):
        fname = '03wcad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Wfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+2)*(ngy+1),offset=0)
    Wfield =  np.transpose(np.reshape(Wfield,(ngy+1,ngx+1,ngz+2))) 
    if component == 2:
        return num,Wfield

    return (num,Ufield,Vfield,Wfield)

def read_icvec_file(fileroot, ngx, ngy, ngz, num):
    """
    Reads U,V,W velocity data from new Fortran ICvec binary file.
    
    The ICvec file stores u, v, w consecutively in memory as flattened arrays.
    
    Parameters
    ----------
    fileroot : str
        Path to directory containing ICvec files.
    ngx, ngy, ngz : int
        Number of grid points in each direction.
    num : int
        Time step / iteration number.
    
    Returns
    -------
    tuple
        (U, V, W) numpy arrays with shapes:
        - U: (ngz-1, ngx-2, ngy-1)
        - V: (ngz-1, ngx-1, ngy-2)
        - W: (ngz-2, ngx-1, ngy-1)
    """
    
    # File name pattern from Fortran code
    filename = f"{fileroot}ICvec.{num:03d}"
    
    # Sizes from Fortran reshapes
    size_u = (ngz-1) * (ngx-2) * (ngy-1)
    size_v = (ngz-1) * (ngx-1) * (ngy-2)
    size_w = (ngz-2) * (ngx-1) * (ngy-1)
    
    # Read entire binary file
    data = np.fromfile(filename, dtype=np.float64)
    
    if data.size != size_u + size_v + size_w:
        raise ValueError(
            f"File {filename} has {data.size} entries, "
            f"expected {size_u + size_v + size_w}"
        )
    
    # Slice out the three fields
    u_flat = data[0:size_u]
    v_flat = data[size_u:size_u+size_v]
    w_flat = data[size_u+size_v:]
    
    # Reshape using Fortran ordering (matches Fortran's RESHAPE)
    U = np.reshape(u_flat, (ngz-1, ngx-2, ngy-1), order='F')
    V = np.reshape(v_flat, (ngz-1, ngx-1, ngy-2), order='F')
    W = np.reshape(w_flat, (ngz-2, ngx-1, ngy-1), order='F')
    
    return num, U, V, W
        
def read_ufvfwf_file(fileroot,ngx,ngy,ngz,num):
    """
    Reads U,V,W velocity flux data from binary fortran files 

    :param fileroot:  Location of data files
    :param ng[x,y,z]: number of grid points in each dimension
    :param num:       Time step to extract data

    :return: The velocity field components as numpy arrays 
        
            - U: The streamwise velocity field components 
            - V: The vertical velocity field components 
            - W: The spanwise velocity field components 
    """

    fname =  '04uu00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Ufield =  np.fromfile(filename,dtype=float,count=(ngx+2)*(ngz+1)*(ngy+1),offset=0)
    Ufield =  np.transpose(np.reshape(Ufield,(ngy+1,ngx+2,ngz+1))) 

    fname = '05vv00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Vfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(ngy+2),offset=0)
    Vfield =  np.transpose(np.reshape(Vfield,(ngy+2,ngx+1,ngz+1))) 

    fname = '06ww00.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Wfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+2)*(ngy+1),offset=0)
    Wfield =  np.transpose(np.reshape(Wfield,(ngy+1,ngx+1,ngz+2))) 

    return (num,Ufield,Vfield,Wfield)

def read_ucvcwcPres_file(fileroot,ngx,ngy,ngz,num,adjointFlag=0,component=-1):
    """
    Reads U,V,W velocity data from binary fortran files 

    :param fileroot:  Location of data files
    :param ng[x,y,z]: number of grid points in each dimension
    :param num:       Time step to extract data
    :param adjointFlag: [1] read adjoint field [else] read forward field
    :param component:   [-1] all, [0] U, [1] V, [2] W, [3] Pressure

    :return: The velocity field components as numpy arrays 
        
            - U: The streamwise velocity field components 
            - V: The vertical velocity field components 
            - W: The spanwise velocity field components 
            - P: The pressure field
    """

    fname    = 'pres_p.dble.'
    if (adjointFlag == 1):
        fname = 'phiadj.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Pfield   = np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(ngy+1),offset=0)
    Pfield   = np.transpose(np.reshape(Pfield,(ngy+1,ngx+1,ngz+1))) 
    if component == 3:
        return num, Pfield

    fname =  '01uc00.dble.'
    if (adjointFlag == 1):
        fname =  '01ucad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Ufield =  np.fromfile(filename,dtype=float,count=(ngx+2)*(ngz+1)*(ngy+1),offset=0)
    Ufield =  np.transpose(np.reshape(Ufield,(ngy+1,ngx+2,ngz+1))) 
    if component == 0:
        return num, Ufield

    fname = '02vc00.dble.'
    if (adjointFlag == 1):
        fname = '02vcad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Vfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(ngy+2),offset=0)
    Vfield =  np.transpose(np.reshape(Vfield,(ngy+2,ngx+1,ngz+1))) 
    if component == 1:
        return num, Vfield

    fname = '03wc00.dble.'
    if (adjointFlag == 1):
        fname = '03wcad.dble.'

    filename = fileroot + fname + "{:07d}".format(num)
    Wfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+2)*(ngy+1),offset=0)
    Wfield =  np.transpose(np.reshape(Wfield,(ngy+1,ngx+1,ngz+2))) 
    if component == 2:
        return num, Wfield

    return (num, Ufield,Vfield,Wfield,Pfield)

def read_obs_file(fileroot,ngx,ngy,ngz,num):
    """
    Read the observed velocity fields from binary files for a given time step.

    Args:
        fileroot (str): The root location of the data files.
        ngx (int): The number of grid points in the x-direction.
        ngy (int): The number of grid points in the y-direction.
        ngz (int): The number of grid points in the z-direction.
        num (int): The time step to extract data for.

    Returns:
        tuple: A tuple containing the observed velocity fields in the x, y, and z directions.
               Each velocity field is a 3D numpy array with dimensions (ngy, ngx, ngz).
    """
    fname = '04uobs.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Ufield =  np.fromfile(filename,dtype=float,count=ngx*ngy*ngz,offset=0)
    Ufield =  np.transpose(np.reshape(Ufield,(ngy,ngx,ngz))) 

    fname = '05vobs.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Vfield =  np.fromfile(filename,dtype=float,count=ngx*ngy*ngz,offset=0)
    Vfield =  np.transpose(np.reshape(Vfield,(ngy,ngx,ngz))) 

    fname = '06wobs.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Wfield =  np.fromfile(filename,dtype=float,count=ngx*ngy*ngz,offset=0)
    Wfield =  np.transpose(np.reshape(Wfield,(ngy,ngx,ngz)))

    return (Ufield,Vfield,Wfield)

def read_boundary_file(fileroot,ngx,ngy,ngz,num,flag=-1,adjFlag=False, adjOpt = 0):
    """
    read boundary data files for/from subdomain Fortran code

    :param fileroot:  Location of data files
    :param ng[x,y,z]: number of grid points in each dimension
    :param num:       Time step to extract data

    :return: A tuple of velocity data on the six subdomain faces 

        - Ulft,Vlft,Wlft:   Velocity Field on Left boundary
        - Urht,Vrht,Wrht:   Velocity Field on Right boundary
        - Utop,Vtop,Wtop:   Velocity Field on Top boundary
        - Ubot,Vbot,Wbot:   Velocity Field on Bottom boundary
        - Ubck,Vbck,Wbck:   Velocity Field on Back boundary
        - Ufrt,Vfrt,Wfrt:   Velocity Field on Right boundary
    """


    if flag == -1 or flag == 1:
        #read bot boundary 
        fname =  'bc_bot.dble.'
        if adjFlag:
            fname = return_gradientBC_prefix(adjOpt) + 'bot.dble.'
            # fname = 'adjbot.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Ubot =  np.fromfile(filename,dtype=float,count=(ngx+2)*(ngz+1)*(1),offset=0)
        Ubot =  np.transpose(np.reshape(Ubot,(1,ngx+2,ngz+1))) 

        Vbot =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(1),offset=(ngz+1)*(ngx+2)*8)
        Vbot =  np.transpose(np.reshape(Vbot,(1,ngx+1,ngz+1))) 

        Wbot =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+2)*(1),offset=(ngz+1)*(ngx+2)*8 + (ngz+1)*(ngx+1)*8)
        Wbot =  np.transpose(np.reshape(Wbot,(1,ngx+1,ngz+2))) 

        if flag ==1:
            return (num,Ubot,Vbot,Wbot)

    if flag == -1 or flag == 2:
        #read top boundary 
        fname =  'bc_top.dble.'
        if adjFlag:
            fname = return_gradientBC_prefix(adjOpt) + 'top.dble.'
            # fname = 'adjtop.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Utop   =  np.fromfile(filename,dtype=float,count=(ngx+2)*(ngz+1)*(1),offset=0)
        Utop =  np.transpose(np.reshape(Utop,(1,ngx+2,ngz+1))) 

        Vtop   =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(1),offset=(ngz+1)*(ngx+2)*8)
        Vtop =  np.transpose(np.reshape(Vtop,(1,ngx+1,ngz+1))) 

        Wtop   =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+2)*(1),offset=(ngz+1)*(ngx+2)*8 + (ngz+1)*(ngx+1)*8)
        Wtop =  np.transpose(np.reshape(Wtop,(1,ngx+1,ngz+2))) 

        if flag==2:
            return (num,Utop,Vtop,Wtop)

    if flag == -1 or flag == 3:
        #read left boundary 
        fname =  'bc_lft.dble.'
        if adjFlag:
            fname = return_gradientBC_prefix(adjOpt) + 'lft.dble.'
            # fname = 'adjlft.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Ulft =  np.fromfile(filename,dtype=float,count=(ngz+1)*(ngy+1),offset=0)
        Ulft =  np.transpose(np.reshape(Ulft,(ngy+1,1,ngz+1))) 

        Vlft =  np.fromfile(filename,dtype=float,count=(ngz+1)*(ngy+2),offset=(ngz+1)*(ngy+1)*8)
        Vlft =  np.transpose(np.reshape(Vlft,(ngy+2,1,ngz+1))) 

        Wlft =  np.fromfile(filename,dtype=float,count=(ngz+2)*(ngy+1),offset=((ngz+1)*(ngy+1) + (ngz+1)*(ngy+2))*8)
        Wlft =  np.transpose(np.reshape(Wlft,(ngy+1,1,ngz+2))) 

        if flag==3:
            return (num,Ulft,Vlft,Wlft)

    if flag == -1 or flag == 4:
        #read right boundary 
        fname =  'bc_rht.dble.'
        if adjFlag:
            fname = return_gradientBC_prefix(adjOpt) + 'rht.dble.'
            # fname = 'adjrht.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Urht =  np.fromfile(filename,dtype=float,count=(1)*(ngz+1)*(ngy+1),offset=0)
        Urht =  np.transpose(np.reshape(Urht,(ngy+1,1,ngz+1))) 

        Vrht =  np.fromfile(filename,dtype=float,count=(1)*(ngz+1)*(ngy+2),offset=(ngz+1)*(ngy+1)*8)
        Vrht =  np.transpose(np.reshape(Vrht,(ngy+2,1,ngz+1))) 

        Wrht =  np.fromfile(filename,dtype=float,count=(1)*(ngz+2)*(ngy+1),offset=(ngz+1)*(ngy+1)*8 + (ngz+1)*(ngy+2)*8)
        Wrht =  np.transpose(np.reshape(Wrht,(ngy+1,1,ngz+2))) 

        if flag == 4:
            return (num,Urht,Vrht,Wrht)

    if flag == -1 or flag == 5:
        #read back boundary 
        fname =  'bc_bck.dble.'
        if adjFlag:
            fname = return_gradientBC_prefix(adjOpt) + 'bck.dble.'
            # fname = 'adjbck.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Ubck =  np.fromfile(filename,dtype=float,count=(ngx+2)*(1)*(ngy+1),offset=0)
        Ubck =  np.transpose(np.reshape(Ubck,(ngy+1,ngx+2,1))) 

        Vbck =  np.fromfile(filename,dtype=float,count=(ngx+1)*(1)*(ngy+2),offset=(ngy+1)*(ngx+2)*8)
        Vbck =  np.transpose(np.reshape(Vbck,(ngy+2,ngx+1,1))) 

        Wbck =  np.fromfile(filename,dtype=float,count=(ngx+1)*(1)*(ngy+1),offset=(ngy+1)*(ngx+2)*8 + (ngx+1)*(ngy+2)*8)
        Wbck =  np.transpose(np.reshape(Wbck,(ngy+1,ngx+1,1))) 

        if flag == 5:
            return (num,Ubck,Vbck,Wbck)

    if flag == -1 or flag == 6:
        #read front boundary
        fname =  'bc_frt.dble.'
        if adjFlag:
            fname = return_gradientBC_prefix(adjOpt) + 'frt.dble.'
            # fname = 'adjfrt.dble.'
        filename = fileroot + fname + "{:07d}".format(num)
        Ufrt =  np.fromfile(filename,dtype=float,count=(ngx+2)*(1)*(ngy+1),offset=0)
        Ufrt =  np.transpose(np.reshape(Ufrt,(ngy+1,ngx+2,1))) 

        Vfrt =  np.fromfile(filename,dtype=float,count=(ngx+1)*(1)*(ngy+2),offset=(ngy+1)*(ngx+2)*8)
        Vfrt =  np.transpose(np.reshape(Vfrt,(ngy+2,ngx+1,1))) 

        Wfrt =  np.fromfile(filename,dtype=float,count=(ngx+1)*(1)*(ngy+1),offset=(ngy+1)*(ngx+2)*8 + (ngx+1)*(ngy+2)*8)
        Wfrt =  np.transpose(np.reshape(Wfrt,(ngy+1,ngx+1,1))) 

        if flag == 6:
            return (num,Ufrt,Vfrt,Wfrt)

    return (num,Ulft,Vlft,Wlft,Urht,Vrht,Wrht,Utop,Vtop,Wtop,Ubot,Vbot,Wbot,Ubck,Vbck,Wbck,Ufrt,Vfrt,Wfrt)
                                
def get_observation_at_boundary(loc, nox, noy, noz, frames_obs, sf_u, sf_v, sf_w, flag=-1):
    """
    Reads the observation data for a given boundary and frame number.

    Parameters
    ----------
    loc : str
        The root location of the observation data files.
    nox, noy, noz : int
        The number of grid points in the x, y, and z directions.
    frames_obs : list
        The list of time frames for which to read the observations.
    sf_u, sf_v, sf_w : array_like
        The flux-to-cartesian conversion factors for the u, v, and w components.
    flag : int, optional
        The flag specifying which boundary to read from. Valid values are:
        - 1: Back wall (x = 0)
        - 2: Front wall (x = Lx)
        - 3: Left wall (y = 0)
        - 4: Right wall (y = Ly)
        - 5: Bottom wall (z = 0)
        - 6: Top wall (z = Lz)

    Returns
    -------
    u_obs, v_obs, w_obs : array_like
        The observed velocity fields in the x, y, and z directions.
    """
    if flag == 1 or flag == 2:
        u_obs = np.zeros((len(frames_obs), noz, nox))
        v_obs = np.zeros_like(u_obs)
        w_obs = np.zeros_like(u_obs)
        for num, frame in enumerate(frames_obs):
                # * Read observation data
            [u_obs_f,v_obs_f,w_obs_f] = read_obs_file(loc,nox,noy,noz,frame)
                # * Convert flux to cartesian observations
            if flag == 1:
                u_obs[num,:,:] = u_obs_f[:,:,0] / sf_u[:,:,0]
                v_obs[num,:,:] = v_obs_f[:,:,0] / sf_v[:,:,0]
                w_obs[num,:,:] = w_obs_f[:,:,0] / sf_w[:,:,0]
            else:
                u_obs[num,:,:] = u_obs_f[:,:,-1] / sf_u[:,:,-1]
                v_obs[num,:,:] = v_obs_f[:,:,-1] / sf_v[:,:,-1]
                w_obs[num,:,:] = w_obs_f[:,:,-1] / sf_w[:,:,-1]
    elif flag == 3 or flag == 4:
        u_obs = np.zeros((len(frames_obs), noz, noy))
        v_obs = np.zeros_like(u_obs)
        w_obs = np.zeros_like(u_obs)
        for num, frame in enumerate(frames_obs):
                # * Read observation data
            [u_obs_f,v_obs_f,w_obs_f] = read_obs_file(loc,nox,noy,noz,frame)
                # * Convert flux to cartesian observations
            if flag == 3:
                u_obs[num,:,:] = u_obs_f[:,0,:] / sf_u[:,0,:]
                v_obs[num,:,:] = v_obs_f[:,0,:] / sf_v[:,0,:]
                w_obs[num,:,:] = w_obs_f[:,0,:] / sf_w[:,0,:]
            else:
                u_obs[num,:,:] = u_obs_f[:,-1,:] / sf_u[:,-1,:]
                v_obs[num,:,:] = v_obs_f[:,-1,:] / sf_v[:,-1,:]
                w_obs[num,:,:] = w_obs_f[:,-1,:] / sf_w[:,-1,:]
    elif flag == 5 or flag == 6:
        u_obs = np.zeros((len(frames_obs), nox, noy))
        v_obs = np.zeros_like(u_obs)
        w_obs = np.zeros_like(u_obs)
        for num, frame in enumerate(frames_obs):
                # * Read observation data
            [u_obs_f,v_obs_f,w_obs_f] = read_obs_file(loc,nox,noy,noz,frame)
                # * Convert flux to cartesian observations
            if flag == 5:
                u_obs[num,:,:] = u_obs_f[0,:,:] / sf_u[0,:,:]
                v_obs[num,:,:] = v_obs_f[0,:,:] / sf_v[0,:,:]
                w_obs[num,:,:] = w_obs_f[0,:,:] / sf_w[0,:,:]
            else:
                u_obs[num,:,:] = u_obs_f[-1,:,:] / sf_u[-1,:,:]
                v_obs[num,:,:] = v_obs_f[-1,:,:] / sf_v[-1,:,:]
                w_obs[num,:,:] = w_obs_f[-1,:,:] / sf_w[-1,:,:]
    else:
        raise Exception("Invalid observation boundary flag value")
    
    return u_obs, v_obs, w_obs

def get_boundary_data(loc, ngx, ngy, ngz, frames, flag, adjFlag = False, adjOpt = 0):
    total_frames = len(frames)
    skip = frames[1] - frames[0]

    if flag == 1 or flag == 2:
        u = np.zeros((total_frames, ngz+1, ngx+2))
        v = np.zeros((total_frames, ngz+1, ngx+1))
        w = np.zeros((total_frames, ngz+2, ngx+1))

        with concurrent.futures.ProcessPoolExecutor() as executor:
            # Submit all jobs
            future_to_frame = {
                executor.submit(read_boundary_file, loc, ngx, ngy, ngz, frame, flag, adjFlag, adjOpt): frame
                for frame in frames
            }

            # Wrap `as_completed` with `tqdm`
            for future in tqdm(concurrent.futures.as_completed(future_to_frame),
                               total=total_frames,
                               desc="Reading boundary data (flag 1/2)"):
                num, U, V, W = future.result()
                u[(num-frames[0])//skip, :, :] = U[:, :, 0]
                v[(num-frames[0])//skip, :, :] = V[:, :, 0]
                w[(num-frames[0])//skip, :, :] = W[:, :, 0]

        return u, v, w
    
    if flag == 3 or flag == 4:
        u = np.zeros((total_frames, ngz+1, ngy+1))
        v = np.zeros((total_frames, ngz+1, ngy+2))
        w = np.zeros((total_frames, ngz+2, ngy+1))

        with concurrent.futures.ProcessPoolExecutor() as executor:
            future_to_frame = {
                executor.submit(read_boundary_file, loc, ngx, ngy, ngz, frame, flag, adjFlag, adjOpt): frame
                for frame in frames
            }

            for future in tqdm(concurrent.futures.as_completed(future_to_frame),
                               total=total_frames,
                               desc="Reading boundary data (flag 3/4)"):
                num, U, V, W = future.result()
                u[(num-frames[0])//skip, :, :] = U[:, 0, :]
                v[(num-frames[0])//skip, :, :] = V[:, 0, :]
                w[(num-frames[0])//skip, :, :] = W[:, 0, :]

        return u, v, w

    if flag == 5 or flag == 6:
        u = np.zeros((total_frames, ngx+2, ngy+1))
        v = np.zeros((total_frames, ngx+1, ngy+2))
        w = np.zeros((total_frames, ngx+1, ngy+1))

        with concurrent.futures.ProcessPoolExecutor() as executor:
            future_to_frame = {
                executor.submit(read_boundary_file, loc, ngx, ngy, ngz, frame, flag, adjFlag, adjOpt): frame
                for frame in frames
            }

            for future in tqdm(concurrent.futures.as_completed(future_to_frame),
                               total=total_frames,
                               desc="Reading boundary data (flag 5/6)"):
                num, U, V, W = future.result()
                u[(num-frames[0])//skip, :, :] = U[0, :, :]
                v[(num-frames[0])//skip, :, :] = V[0, :, :]
                w[(num-frames[0])//skip, :, :] = W[0, :, :]

        return u, v, w

def return_gradientBC_prefix(case):
    if case == 0: 
      fname =  'adj'
    elif case == 1:
      fname =  'ap1'
    elif case == 2:
      fname =  'ap2'
    elif case == 3:
      fname =  'ap3'
    elif case == 4:
      fname =  'ap4'
    elif case == 5:
      fname =  'ap5'

    return fname

def write_I_plane_data_ASCII_tecplot(file_loc, u, v, w, p, z, y, i_ind=0, t_ind=0, time=0.0):
    """
    *------------------------------------------------------------*
    | Write a Z-Y plane slice (I-plane) of 3D flow data to an    |
    | ASCII Tecplot-readable file, with solution time and strand |
    *------------------------------------------------------------*
    
    Parameters:
    ----------
    file_loc : str
        # --> Directory to save the output file (end with '/')
    u, v, w : 2D np.ndarray
        # --> Velocity components in the Z-Y plane (shape: [nz, ny])
    p : 2D np.ndarray
        # --> Pressure values in the Z-Y plane (shape: [nz, ny])
    z, y : 1D np.ndarray
        # --> Grid coordinates along z and y directions
    i_ind : int
        # --> Index of the I-plane (used in filename and zone title)
    t_ind : int
        # --> Time step index (used in filename and zone title)
    time : float
        # --> Physical solution time (used as Tecplot SOLUTIONTIME)

    Returns:
    -------
    None (writes an ASCII .dat file for Tecplot)
    """

    # * Set the output filename
    filename = file_loc + f"I_{i_ind:03d}_T_{t_ind:03d}.dat"
    
    # * Create meshgrid for Z (horizontal) and Y (vertical) coordinates
    # ^ Tecplot expects Z-Y data for a vertical slice
    Z, Y = np.meshgrid(z, y, indexing='ij')  # --> Shape: (nz, ny)
    
    # * Open file and start writing Tecplot-formatted ASCII data
    with open(filename, "w") as f:
        # * Write the title with I-plane info
        f.write(f"TITLE = \"Data at i = {i_ind} \"\n")
        
        # * List the variables in the correct order for Tecplot
        f.write("VARIABLES = \"Z\", \"Y\", \"U\", \"V\", \"W\", \"P\"\n")
        
        # * Write the zone info including dimensions, format, strand ID, and solution time
        f.write(
            f"ZONE I={len(y)}, J={len(z)}, F=POINT, "
            f"STRANDID={i_ind}, SOLUTIONTIME={time:.2f}, T=\"t_{t_ind}\"\n"
        )
        
        # * Write one line per data point (row-major order)
        for iz in range(len(z)):
            for iy in range(len(y)):
                f.write(
                    f"{Z[iz, iy]:.6f} {Y[iz, iy]:.6f} "
                    f"{u[iz, iy]:.6e} {v[iz, iy]:.6e} "
                    f"{w[iz, iy]:.6e} {p[iz, iy]:.6e}\n"
                )
    
    # ^ Optional: uncomment for debug confirmation
    # print(f"--> Data written to {filename}")

def write_J_plane_data_ASCII_tecplot(file_loc, u, v, w, p, x, z, j_ind=0, t_ind=0, time=0.0):
    """
    *------------------------------------------------------------*
    | Write an X-Z plane slice (J-plane) of 3D flow data to an   |
    | ASCII Tecplot-readable file, with solution time and strand |
    *------------------------------------------------------------*

    Parameters:
    ----------
    file_loc : str
        # --> Directory to save the output file (should end with '/')
    u, v, w : 2D np.ndarray
        # --> Velocity components in the X-Z plane (shape: [nx, nz])
    p : 2D np.ndarray
        # --> Pressure values in the X-Z plane (shape: [nx, nz])
    x, z : 1D np.ndarray
        # --> Grid coordinates along x and z directions
    j_ind : int
        # --> Index of the J-plane (used in filename and zone title)
    t_ind : int
        # --> Time step index (used in filename and zone title)
    time : float
        # --> Physical solution time (used as Tecplot SOLUTIONTIME)

    Returns:
    -------
    None (writes an ASCII .dat file for Tecplot)
    """

    # * Set the output filename
    filename = file_loc + f"J_{j_ind:03d}_T_{t_ind:03d}.dat"
    
    # * Create meshgrid for X-Z coordinates
    # ^ Shape will be (nx, nz) matching the plane data
    X, Z = np.meshgrid(x, z, indexing='ij')  # --> Shape: (nx, nz)
    
    # * Open the output file and begin writing Tecplot format
    with open(filename, "w") as f:
        # * Write the plot title with plane index
        f.write(f"TITLE = \"Data at j = {j_ind} \"\n")
        
        # * Declare variable names
        f.write("VARIABLES = \"X\", \"Z\", \"U\", \"V\", \"W\", \"P\"\n")
        
        # * Write zone information with I, J sizes, strand, and time
        f.write(
            f"ZONE I={len(z)}, J={len(x)}, F=POINT, "
            f"STRANDID={j_ind}, SOLUTIONTIME={time:.6f}, T=\"t_{t_ind}\"\n"
        )
        
        # * Loop through and write data for each (x, z) point
        for ix in range(len(x)):
            for iz in range(len(z)):
                f.write(
                    f"{X[ix, iz]:.6f} {Z[ix, iz]:.6f} "
                    f"{u[ix, iz]:.6e} {v[ix, iz]:.6e} "
                    f"{w[ix, iz]:.6e} {p[ix, iz]:.6e}\n"
                )

    # ^ Optional: Uncomment for debugging
    # print(f"--> Data written to {filename}")

def write_K_plane_data_ASCII_tecplot(file_loc, u, v, w, p, x, y, k_ind=0, t_ind=0, time=0.0):
    """
    *-------------------------------------------------------------*
    | Write an X-Y plane slice (K-plane) of 3D flow data to an    |
    | ASCII Tecplot-readable file, with solution time and strand  |
    *-------------------------------------------------------------*

    Parameters:
    ----------
    file_loc : str
        # --> Directory to save the output file (should end with '/')
    u, v, w : 2D np.ndarray
        # --> Velocity components in the X-Y plane (shape: [nx, ny])
    p : 2D np.ndarray
        # --> Pressure values in the X-Y plane (shape: [nx, ny])
    x, y : 1D np.ndarray
        # --> Grid coordinates along x and y directions
    k_ind : int
        # --> Index of the K-plane (used in filename and zone title)
    t_ind : int
        # --> Time step index (used in filename and zone title)
    time : float
        # --> Physical solution time (used as Tecplot SOLUTIONTIME)

    Returns:
    -------
    None (writes an ASCII .dat file for Tecplot)
    """

    # * Set the output filename
    filename = file_loc + f"K_{k_ind:03d}_T_{t_ind:03d}.dat"
    
    # * Create meshgrid for X-Y coordinates
    # ^ Shape will be (nx, ny) matching the plane data
    X, Y = np.meshgrid(x, y, indexing='ij')  # --> Shape: (nx, ny)
    
    # * Open the output file and begin writing Tecplot format
    with open(filename, "w") as f:
        # * Write the plot title with plane index
        f.write(f"TITLE = \"Data at k = {k_ind} \"\n")
        
        # * Declare variable names
        f.write("VARIABLES = \"X\", \"Y\", \"U\", \"V\", \"W\", \"P\"\n")
        
        # * Write zone information with I, J sizes, strand, and time
        f.write(
            f"ZONE I={len(y)}, J={len(x)}, F=POINT, "
            f"STRANDID={k_ind}, SOLUTIONTIME={time:.6f}, T=\"t_{t_ind}\"\n"
        )
        
        # * Loop through and write data for each (x, y) point
        for ix in range(len(x)):
            for iy in range(len(y)):
                f.write(
                    f"{X[ix, iy]:.6f} {Y[ix, iy]:.6f} "
                    f"{u[ix, iy]:.6e} {v[ix, iy]:.6e} "
                    f"{w[ix, iy]:.6e} {p[ix, iy]:.6e}\n"
                )

    # ^ Optional: Uncomment for debugging
    # print(f"--> Data written to {filename}")

def write_wall_data_ASCII_tecplot(file_loc, Sxy, Szy, Pw, x, z, t_ind=0, time=0.0):

    # * Set the output filename
    filename = file_loc + f"wall_T_{t_ind:03d}.dat"
    
    # * Create meshgrid for X-Z coordinates
    # ^ Shape will be (nx, nz) matching the plane data
    X, Z = np.meshgrid(x, z, indexing='ij')  # --> Shape: (nx, nz)
    
    # * Open the output file and begin writing Tecplot format
    with open(filename, "w") as f:
        # * Write the plot title with plane index
        f.write(f"TITLE = \"Data at wall \"\n")
        
        # * Declare variable names
        f.write("VARIABLES = \"X\", \"Z\", \"Sxy\", \"Szy\",\"Pw\"\n")
        
        # * Write zone information with I, J sizes, strand, and time
        f.write(
            f"ZONE I={len(z)}, J={len(x)}, F=POINT, "
            f"STRANDID={1}, SOLUTIONTIME={time:.6f}, T=\"t_{t_ind}\"\n"
        )
        
        # * Loop through and write data for each (x, z) point
        for ix in range(len(x)):
            for iz in range(len(z)):
                f.write(
                    f"{X[ix, iz]:.6f} {Z[ix, iz]:.6f} "
                    f"{Sxy[ix, iz]:.6e} {Szy[ix, iz]:.6e} {Pw[ix, iz]:.6e}\n"
                )

    # ^ Optional: Uncomment for debugging
    # print(f"--> Data written to {filename}")

def write_obs_file(U,V,W,num=0,dir='./'):
    """
    Writes observation file for subdomain
 
    :param [U,V,W]: velocity components (should be in flux form)
    :param num:     timestep number to read in 

    :return: None

    """

    fname = '04uobs.dble.'+"{:07d}".format(num)
    fout = open(dir + fname,'wb')
    fout.write(bytes(np.transpose(U)))
    fout.close()

    fname = '05vobs.dble.'+"{:07d}".format(num)
    fout = open(dir + fname,'wb')
    fout.write(bytes(np.transpose(V)))
    fout.close()

    fname = '06wobs.dble.'+"{:07d}".format(num)
    fout = open(dir + fname,'wb')
    fout.write(bytes(np.transpose(W)))
    fout.close()

def read_DA_statistics(file_addr):
    u_mean     = np.loadtxt(file_addr + 'u_mean.txt')
    v_mean     = np.loadtxt(file_addr + 'v_mean.txt')
    w_mean     = np.loadtxt(file_addr + 'w_mean.txt')
    p_mean     = np.loadtxt(file_addr + 'w_mean.txt')
    uu         = np.loadtxt(file_addr + 'uu.txt')
    vv         = np.loadtxt(file_addr + 'vv.txt')
    ww         = np.loadtxt(file_addr + 'ww.txt')
    uv         = np.loadtxt(file_addr + 'uv.txt')
    pp         = np.loadtxt(file_addr + 'pp.txt')
    return u_mean, v_mean, w_mean, p_mean, uu, vv, ww, uv, pp

def read_exp_statistics(file_addr):
    """
    Reads experimental statistics from files and returns them.

    Parameters
    ----------
    file_addr: str
        Path to the directory containing the experimental statistics files.

    Returns
    -------
    u_mean_exp, v_mean_exp, w_mean_exp : 1D numpy arrays
        Mean velocities in x, y, and z directions.
    uu_exp, vv_exp, ww_exp : 1D numpy arrays
        Reynolds stresses in x, y, and z directions.
    uv_exp : 1D numpy array
        Reynolds shear stress in x and y directions.
    p_var_exp : 1D numpy array
        Turbulent kinetic energy.
    """
    u_mean_exp      = np.loadtxt(file_addr + 'u_mean.txt')
    v_mean_exp      = np.loadtxt(file_addr + 'v_mean.txt')
    w_mean_exp      = np.loadtxt(file_addr + 'w_mean.txt')
    uu_exp          = np.loadtxt(file_addr + 'uu.txt')
    vv_exp          = np.loadtxt(file_addr + 'vv.txt')
    ww_exp          = np.loadtxt(file_addr + 'ww.txt')
    uv_exp          = np.loadtxt(file_addr + 'uv.txt')
    # pp_exp          = np.loadtxt(file_addr + 'pp.txt')

    # return u_mean_exp, v_mean_exp, w_mean_exp, uu_exp, vv_exp, ww_exp, uv_exp, pp_exp
    return u_mean_exp, v_mean_exp, w_mean_exp, uu_exp, vv_exp, ww_exp, uv_exp

def write_obs_fullsized_file(U,V,W,num=0,addr='./'):
    """
    Writes ucvcwc restart file at a given timestep.

    :param addr:   name of the folder to write to
    :param [U,V,W]: velocity components
    :param num:     number of the file to read (default 0)
    """
    # Note: main code expects z,x,y ordering

    fname = addr+'01ucmn.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(U)))
    fout.close()

    fname = addr + '02vcmn.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(V)))
    fout.close()

    fname = addr + '03wcmn.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(W)))
    fout.close()

def read_obs_fullsized_file(fileroot, ngx, ngy, ngz, num=0, flag=1):
    if flag == 1:
        fname =  '01ucob.dble.'
    else: 
        fname =  '01ucmn.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Ufield =  np.fromfile(filename,dtype=float,count=(ngx+2)*(ngz+1)*(ngy+1),offset=0)
    Ufield =  np.transpose(np.reshape(Ufield,(ngy+1,ngx+2,ngz+1))) 
    if flag == 1:
        fname =  '02vcob.dble.'
    else: 
        fname =  '02vcmn.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Vfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+1)*(ngy+2),offset=0)
    Vfield =  np.transpose(np.reshape(Vfield,(ngy+2,ngx+1,ngz+1)))
    if flag == 1:
        fname =  '03wcob.dble.'
    else: 
        fname =  '03wcmn.dble.'
    filename = fileroot + fname + "{:07d}".format(num)
    Wfield =  np.fromfile(filename,dtype=float,count=(ngx+1)*(ngz+2)*(ngy+1),offset=0)
    Wfield =  np.transpose(np.reshape(Wfield,(ngy+1,ngx+1,ngz+2)))
    return Ufield, Vfield, Wfield

def write_obsKern(Ku,Kv,Kw,addr='./'):
    num = 0
    # Note: main code expects z,x,y ordering

    fname = addr+'obKeru.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Ku)))
    fout.close()

    fname = addr + 'obKerv.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Kv)))
    fout.close()

    fname = addr + 'obKerw.dble.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Kw)))
    fout.close()

    print(f'Observation kernel written to {addr}obKer*.dble.{num}')

def write_flags(Fu,Fv,Fw,addr='./'):
    num = 0
    # Note: main code expects z,x,y ordering

    fname = addr + 'obFlgu.intg.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Fu)))
    fout.close()
    fname = addr + 'obFlgv.intg.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Fv)))
    fout.close()
    fname = addr + 'obFlgw.intg.'+"{:07d}".format(num)
    fout = open(fname,'wb')
    fout.write(bytes(np.transpose(Fw)))
    fout.close()

    print(f'Observation flags written to {addr}obFlg*.intg.{num}')

def read_obsKern(addr, ngx, ngy, ngz, num=0):
    # Note: main code expects z,x,y ordering

    fname = addr + 'obKeru.dble.' + "{:07d}".format(num)
    Ku =  np.fromfile(fname,dtype=float,count=(ngx+2)*(ngz+1)*(ngy+1),offset=0)
    Ku =  np.transpose(np.reshape(Ku,(ngy+1,ngx+2,ngz+1))) 

    fname = addr + 'obKerv.dble.' + "{:07d}".format(num)
    Kv =  np.fromfile(fname,dtype=float,count=(ngx+1)*(ngz+1)*(ngy+2),offset=0)
    Kv =  np.transpose(np.reshape(Kv,(ngy+2,ngx+1,ngz+1)))

    fname = addr + 'obKerw.dble.' + "{:07d}".format(num)
    Kw =  np.fromfile(fname,dtype=float,count=(ngx+1)*(ngz+2)*(ngy+1),offset=0)
    Kw =  np.transpose(np.reshape(Kw,(ngy+1,ngx+1,ngz+2)))

    print(f'Observation kernel read from {addr}obKer*.dble.{num}')
    return Ku, Kv, Kw

def read_flags(addr, ngx, ngy, ngz, num=0):

    fname = addr + 'obFlgu.intg.' + "{:07d}".format(num)
    uflag = np.fromfile(fname,dtype=np.int32,count=(ngx+2)*(ngz+1)*(ngy+1),offset=0)
    uflag = np.transpose(np.reshape(uflag,(ngy+1,ngx+2,ngz+1)))

    fname = addr + 'obFlgv.intg.' + "{:07d}".format(num)
    vflag = np.fromfile(fname,dtype=np.int32,count=(ngx+1)*(ngz+1)*(ngy+2),offset=0)
    vflag = np.transpose(np.reshape(vflag,(ngy+2,ngx+1,ngz+1)))

    fname = addr + 'obFlgw.intg.' + "{:07d}".format(num)
    wflag = np.fromfile(fname,dtype=np.int32,count=(ngx+1)*(ngz+2)*(ngy+1),offset=0)
    wflag = np.transpose(np.reshape(wflag,(ngy+1,ngx+1,ngz+2)))

    print(f'Observation flags read from {addr}obKer*.flag.{num}')
    return uflag, vflag, wflag
