import matplotlib.pyplot as plt


def plot_solution_panels(sol):
    fig, axs = plt.subplots(2, 2)

    axs[0, 0].plot(sol.t, sol.y[0])
    axs[0, 0].set_title('x against t')
    axs[0, 0].set_ylabel('x (m)')
    axs[0, 0].grid()

    axs[0, 1].plot(sol.t, sol.y[1], 'tab:orange')
    axs[0, 1].set_title('phi against t')
    axs[0, 1].set_ylabel('phi (rads)')
    axs[0, 1].grid()

    axs[1, 0].plot(sol.t, sol.y[2], 'tab:green')
    axs[1, 0].set_title('v against t')
    axs[1, 0].set_ylabel('v (m/s)')
    axs[1, 0].grid()

    axs[1, 1].plot(sol.t, sol.y[3], 'tab:red')
    axs[1, 1].set_title('w against t')
    axs[1, 1].set_ylabel('w (rads/s)')
    axs[1, 1].grid()

    fig.tight_layout()
    plt.show()
