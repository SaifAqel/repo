# io/profile_writer.py
import pandas as pd
def export_profiles(stage_result, path: str):
    df = pd.DataFrame([c.__dict__ for c in stage_result.cells])
    df.to_parquet(path) if path.endswith(".parquet") else df.to_csv(path,index=False)
