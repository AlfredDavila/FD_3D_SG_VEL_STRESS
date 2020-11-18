import numpy as np
import matplotlib.pyplot as plt
#https://numpydoc.readthedocs.io/en/latest/format.html


def diff_2(field, N, dx, order):
    """Made the derivate of a plane in one direction

    1. First derivate all the field 
    2. Then corrects the boundary values
    3. Return the derivate field

    Function used in the finite differences formulation
    for the seismic elastic wave equation

    Parameters
    ----------
    field : ndarray
        The field (Stress, Velocity ....)
    N : int
        Number of elements per row
    dx : number
        Delta x in the derivation
    order : int 
        Order of the differential operator

    Returns
    -------
    ndarray
        A ndarray with the second derivate of the field
    """
	
    #Validate order
    if order not in [2,4]:

        raise ValueError('Order not valid: {0}'.format(order)) 
        

    #Firts derivate all the field
    indexes = np.arange(0, len(field))

    prev_idxs = indexes - 1

    post_idxs = indexes + 1

    post_idxs[-1] = 0

    #For order 4
    #prev_prev_idxs = indexes - 2

    post_post_idxs = indexes + 2

    post_post_idxs[-2:] = 0

    c0 = 9.0/8.0

    c1 = 1.0/24.0

    dx_inv = 1/(dx)

    if order == 2:

        field_diff = field[post_idxs] - field

        print(post_idxs[10], '-',prev_idxs[10])

    elif order == 4:
        # @todo no da buenos resultados, revisar

        field_diff = (c0 * (field[post_idxs] - field) )
        - ( c1 *(field[post_post_idxs] - field[prev_idxs]) )

        print(post_idxs[10], '-',10, post_post_idxs[10], '-', prev_idxs[10])

    ##Then corrects the boundary values
    boundary_idx = np.arange(0, len(field), N)

    # In order 2 the boundary index are iN and iN-1 (i integer)

    boundary_idx_prev = boundary_idx - 1

    field_diff[boundary_idx] = 0.0
    field_diff[boundary_idx_prev] = 0.0

    # In order 4 the boundary index are the same of order 2 plus iN+1 and iN-2 (i integer)
    if order == 4:
        boundary_idx_prev_prev = boundary_idx - 2
        boundary_idx_post = boundary_idx + 1
        boundary_idx_post[-1] = 0

        field_diff[boundary_idx_prev_prev] = 0.0
        field_diff[boundary_idx_post] = 0.0

    #Return the derivate field
    field_diff = dx_inv * field_diff

    return field_diff

dx = 0.05

x = np.arange(0,2 * np.pi,dx)

y = np.sin(x)

y_diff = np.cos(x)

field_diff_a = y_diff

field = y

N_rows = 3

for i in range(N_rows-1):
    field = np.concatenate((field,y))

    field_diff_a = np.concatenate((field_diff_a,y_diff))

field_diff_4 = diff_2(field, len(x), dx, 4)

field_diff_2 = diff_2(field, len(x), dx, 2)

x_plot = np.arange(0, dx * len(x) * N_rows, dx)

plt.plot(x_plot, field_diff_2)
plt.plot(x_plot, field_diff_4)
plt.plot(x_plot, field_diff_a)
plt.show()
