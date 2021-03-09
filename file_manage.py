import pickle


def save_file(sol_arr, name):
    """
    Saves a solution to a file, to make it faster to plot graphs by using data from the file rather than calculating
    a solution each time.
    Arguments:
        sol_arr: The array of solutions calculated

        name: The name of the file to be saved to

    """

    with open(name + '.data', 'wb') as filehandle:
        pickle.dump(sol_arr, filehandle)


def open_file(name):
    """
    Opens a .data file containing a set of solutions and stores it in sol_arr.
    Arguments:
        sol_arr: The array of solutions calculated previously.

        name: The name of the file to be opened, as a string, (not including the .data).
    Returns:
        sol_arr: The array of the solutions saved in the file.

    """
    with open(name + '.data', 'rb') as filehandle:
        sol_arr = pickle.load(filehandle)
    return sol_arr
