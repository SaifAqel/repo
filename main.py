from heat_transfer.functions.runner import run

if __name__ == "__main__":
    stages_path = "heat_transfer/config/stages.yaml"
    streams_path = "heat_transfer/config/streams.yaml"
    gas_hist, water_hist = run(stages_path, streams_path)
    g_out = gas_hist[-1]
    w_out = water_hist[-1]
    print("Gas T_out:", g_out.temperature)
    print("Gas P_out:", g_out.pressure)
    print("Water h_out:", w_out.enthalpy)