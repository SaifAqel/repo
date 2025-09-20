# reports/plots.py
import matplotlib.pyplot as plt
def plot_profiles(stage_result):
    xs = [c.x for c in stage_result.cells]
    Ts = [c.Tg for c in stage_result.cells]
    Ps = [c.Pg for c in stage_result.cells]
    qs = [c.qpp for c in stage_result.cells]
    plt.figure(); plt.plot(xs,Ts); plt.xlabel("x"); plt.ylabel("T")
    plt.figure(); plt.plot(xs,Ps); plt.xlabel("x"); plt.ylabel("P")
    plt.figure(); plt.plot(xs,qs); plt.xlabel("x"); plt.ylabel("q''")
