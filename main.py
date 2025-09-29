segments = [
    {"type": "pass", "solver": solver_pass1, "z_eval": np.linspace(0, pass1.L, 100)},
    {"type": "nozzle", "nozzle": rev1},
    {"type": "pass", "solver": solver_pass2, "z_eval": np.linspace(0, pass2.L, 100)},
    {"type": "nozzle", "nozzle": rev2},
    {"type": "pass", "solver": solver_pass3, "z_eval": np.linspace(0, pass3.L, 100)}
]
