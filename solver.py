from methods import *
from scipy.integrate import solve_ivp
from file_manage import *


def single_sol(initial_conds, t_range, args, max_step_size):
    """
    Solves an IVP given a set of initial values and time step.

    Arguments:
        initial_conds: A vector of the initial conditions:
                            initial_conds = [x0,phi0,v0,w0]

        t_range: A tuple containing the time range:
                            t_range = (t_0, t_end)

        args: A vector of the constants to be used in the integration:
                            args = [m1, m2, l, mu, g, r]

        max_step_size: The maximum step between calculations, the smaller the more precise.

    Returns:
        sol: The solution of the differential equation, type OdeSolution.

    """
    print("Integrating...")
    sol = solve_ivp(diff_func_P, t_range, initial_conds, method='RK45', args=([args]), dense_output=True,
                    max_step=max_step_size)
    return sol


def multi_sol(initial_conds_arr, t_range, args_arr, max_step_size, save_to_file, save_file_name):
    """
        Solves an IVP for an array of initial values and constants.

        Arguments:
            initial_conds_arr: An array of vectors of the initial conditions:
                                initial_conds_arr = [[x0,phi0,v0,w0],[x1,phi1,v1,w1],...]

            t_range: A tuple containing the time range:
                                t_range = (t_0, t_end)

            args_arr: An array of vectors of the constants to be used in the integration:
                                args = [[m1, m2, l, mu, g, r],[m1, m2, l, mu, g, r],...]

            max_step_size: The maximum step between calculations, the smaller the more precise.

            save_to_file: A bool (True or False) whether the solution array is saved to a file.

            save_file_name: The filename that the solution is saved to.

        Returns:
            sol_arr: The array of solutions of the differential equation, type OdeSolution.

        """

    if save_to_file is True:
        if type(save_file_name) is not str:
            raise Exception("You must give a filename (string) to save the data to")

    if len(initial_conds_arr) != len(args_arr):
        raise Exception("The initials conditions are not the same size as the arguments!")

    sol_arr = []

    for idx, conds in enumerate(initial_conds_arr):
        print("set ", idx+1)

        sol_arr.append(single_sol(conds, t_range, args_arr[idx], max_step_size))

    if save_to_file is True:
        save_file(sol_arr, save_file_name)

    return sol_arr



