#!/usr/bin/env python3
import numpy as np
from solver import *
from plot import *
from methods import *

"""

    You can find descriptions of what each function does and how to use it under it's declaration.
    
    Examples of plots:
    
        Example 1: Will plot graphs with 3 solutions of differing m2 values:
    
            args_1 = [m1, 40, l, mu, g, r]
            args_2 = [m1, 60, l, mu, g, r]
            args_3 = [m1, 80, l, mu, g, r]
            
            initial_conds_1 = [x0, phi0, v0, w0]
            initial_conds_2 = [x0, phi0, v0, w0]
            initial_conds_3 = [x0, phi0, v0, w0]
            
            args_arr = [args_1, args_2, args_3]
            initial_conds_arr = [initial_conds_1, initial_conds_2, initial_conds_3]
            
            ########### The labels of the legend are in this list ###########
            
            label_legend_list = ["$m2 = 40$", "$m2 = 60$", "$m2 = 80$"]
    
            #################################################################        
                            
            ### The file will be saved to name 'different m2 values.data' ###
            
            save_file_name = 'different m2 values'
            
            #################################################################
            
            t_range = (t_0, t_end)
            max_step_size = 0.001
            
            # If you don't want to save the solutions to a file, set the True to False and leave the save_file_name 
            empty. 
            
            sol = multi_sol(initial_conds_arr, t_range, args_arr, max_step_size, True, save_file_name)
            
            
            plot_multi_sol_panels(sol, label_legend_list)
            plot_multi_spiral(sol, label_legend_list, args_arr)           

        Example 2: Will make a 3 plots of a single solution:
        
            args = [m1, m2, l, mu, g, r]
            initial_conds = [x0, phi0, v0, w0]
            t_range = (t_0, t_end)
            max_step_size = 0.1 # Decrease this for more accuracy/precision if the spiral plot looks disjointed
            
            sol = single_sol(initial_conds_1, t_range, args_1, max_step_size)
            
            plot_solution_panels(sol)
            plot_colour_spiral(sol, l, r)
            plot_spiral(sol, l, r)


"""



def main():
    args_1 = [m1, 40, l, mu, g, r]
    args_2 = [m1, 60, l, mu, g, r]
    args_3 = [m1, 80, l, mu, g, r]

    initial_conds_1 = [x0, phi0, v0, w0]
    initial_conds_2 = [x0, phi0, v0, w0]
    initial_conds_3 = [x0, phi0, v0, w0]

    args_arr = [args_1, args_2, args_3]
    initial_conds_arr = [initial_conds_1, initial_conds_2, initial_conds_3]

    ########### The labels of the legend are in this list ###########

    label_legend_list = ["$m2 = 40$", "$m2 = 60$", "$m2 = 80$"]

    #################################################################

    ### The file will be saved to name 'different m2 values.data' ###

    save_file_name = 'different m2 values'

    #################################################################

    t_range = (t_0, t_end)
    max_step_size = 0.001

    # If you don't want to save the solutions to a file, set the True to False and leave the save_file_name empty

    sol = multi_sol(initial_conds_arr, t_range, args_arr, max_step_size, True, save_file_name)

    plot_multi_sol_panels(sol, label_legend_list)
    plot_multi_spiral(sol, label_legend_list, args_arr)


if __name__ == '__main__':
    # Initial Conditions
    x0 = 0
    phi0 = np.pi / 2
    v0 = 0
    w0 = 0
    # Constant values
    m1 = 600  # Heavier mass
    m2 = 30  # Lighter mass
    l = 20
    mu = 0.5
    g = 9.81
    r = 1  # Radius of the "pulley"
    # Integration region
    t_0 = 0
    t_end = 10

    main()
