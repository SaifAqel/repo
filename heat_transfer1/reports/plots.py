# reports/plots.py
import matplotlib.pyplot as plt

def plot_profiles(stage_result):
    xs = [c.x.to("meter").magnitude for c in stage_result.cells]
    Ts = [c.Tg.to("kelvin").magnitude for c in stage_result.cells]
    Ps = [c.Pg.to("pascal").magnitude for c in stage_result.cells]
    qs = [c.qpp.to("watt/meter**2").magnitude for c in stage_result.cells]

    plt.figure()
    plt.plot(xs, Ts)
    plt.xlabel("x [m]")
    plt.ylabel("T [K]")

    plt.figure()
    plt.plot(xs, Ps)
    plt.xlabel("x [m]")
    plt.ylabel("P [Pa]")

    plt.figure()
    plt.plot(xs, qs)
    plt.xlabel("x [m]")
    plt.ylabel("q'' [W/m²]")
