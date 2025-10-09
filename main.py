from heat_transfer.functions.runner import run

if __name__ == "__main__":
    stages_path = "heat_transfer/config/stages.yaml"
    streams_path = "heat_transfer/config/streams.yaml"

    result = run(stages_path, streams_path)

    print("Converged:", result["converged"])
    print("Iterations:", result["iterations"])
    print("hw0:", result["hw0"])
    print("Residual:", result["residual"])
