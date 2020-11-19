import numpy as np

from findiff import FinDiff

def update_velocity(v_x, v_y, v_z, t_xx, t_yy, t_zz, t_xy, t_xz, t_yz, 
	dt, d_dx, d_dy, d_dz, f_x, f_y, f_z, rho):
	"""Update the velocity field in the finite diferences method at the time n+1/2

	Parameters
    ----------
    v_x : ndarray 
    	3D ndarray with the velocity field in x direction
    v_y : ndarray 
    	3D ndarray with the velocity field in y direction
    v_z : ndarray 
    	3D ndarray with the velocity field in z direction
    t_xx : ndarray 
    	3D ndarray with the stress field xx component
    t_yy : ndarray 
    	3D ndarray with the stress field yy component
    t_zz : ndarray 
    	3D ndarray with the stress field zz component
    t_xy : ndarray 
    	3D ndarray with the stress field xy component
    t_xz : ndarray 
    	3D ndarray with the stress field xz component
    t_yz : ndarray 
    	3D ndarray with the stress field yz component
	dt : number
        Time step
	d_dx : number
        Differential operator in x direction
	d_dy : FinDiff
		Differential operator in y direction
	d_dz : number
        Differential operator in z direction
	f_x : ndarray 
    	3D ndarray with body force field in x directon
	f_y : ndarray 
    	3D ndarray with body force field in y directon
	f_z : ndarray 
    	3D ndarray with body force field in z directon
    rho : ndarray 
    	3D ndarray with values of density effective

    Returns
    -------
    list
        A list with the ndarrays of the velocity fields updated []
    """

    #Effective parameters see Randall et al. (1991),

	v_x = v_x + dt * (1 / rho) * (d_dx(t_xx) + d_dy(t_xy) + d_dz(t_xz) + f_x[it])

    v_y = v_y + dt * (1 / rho) * (d_dx(t_xy) + d_dy(t_yy) + d_dz(t_yz) + f_y[it])

    v_z = v_z + dt * (1 / rho) * (d_dx(t_xz) + d_dy(t_yz) + d_dz(t_zz) + f_z[it])

	return [v_x, v_y, v_z]

def update_stress(v_x, v_y, v_z, dt, d_dx, d_dy, d_dz, rho):
	"""Update the stress tensor field in the finite diferences method at the time n+1

	Parameters
    ----------
    v_x : ndarray 
    	3D ndarray with the velocity field in x direction
    v_y : ndarray 
    	3D ndarray with the velocity field in y direction
    v_z : ndarray 
    	3D ndarray with the velocity field in z direction
    t_xx : ndarray 
    	3D ndarray with the stress field xx component
    t_yy : ndarray 
    	3D ndarray with the stress field yy component
    t_zz : ndarray 
    	3D ndarray with the stress field zz component
    t_xy : ndarray 
    	3D ndarray with the stress field xy component
    t_xz : ndarray 
    	3D ndarray with the stress field xz component
    t_yz : ndarray 
    	3D ndarray with the stress field yz component
	dt : number
        Time step
	d_dx : number
        Differential operator in x direction
	d_dy : FinDiff
		Differential operator in y direction
	d_dz : number
        Differential operator in z direction
	rho : ndarray 
    	3D ndarray with values of density effective
    lambda_1 : ndarray 
    	3D ndarray with values of lambda effective (Lamé parameter)
    mu : ndarray 
    	3D ndarray with values of mu effective (Lamé parameter)
    lambda_2mu : ndarray 
    	3D ndarray with values of lambda + 2 mu effective (Lamé parameters)

    Returns
    -------
    list
        A list with the ndarrays of the stress tensor fields updated []
    """

    #Effective parameters see Randall et al. (1991),
    
	t_xx = t_xx + dt * ( lambda_2mu * d_dx(v_x) + lambda_l * (d_dy(v_y) + d_dz(v_z)))
    
    t_yy = t_yy + dt * ( lambda_2mu * d_dy(v_y) + lambda_l * (d_dx(v_x) + d_dz(v_z)))
    
    t_zz = t_zz + dt * ( lambda_2mu * d_dz(v_z) + lambda_l * (d_dx(v_x) + d_dy(v_y)))
    
    t_xy = t_xy + dt * ( mu * (d_dy(v_z) + d_dx(v_y)))
    
    t_xz = t_xz + dt * ( mu * (d_dz(v_x) + d_dx(v_x)))
    
    t_yz = t_yz + dt * ( mu * (d_dx(v_y) + d_dy(v_z)))

	return [v_x, v_y, v_z]

