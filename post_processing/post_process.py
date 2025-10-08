# heat_transfer/post_processing/post_process.py
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import pydot
from io import BytesIO
from heat_transfer.runner.ChainStages import GasStreamProfile, WaterStreamProfile
from common.units import Q_, ureg

def create_dataframes(gas_profile: GasStreamProfile, water_profile: WaterStreamProfile) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Convert profiles to Pandas DataFrames for easier calcs/plotting (strip units to magnitudes)."""
    gas_data = {
        'z': [p.z.magnitude for p in gas_profile.points],
        'temperature': [p.temperature.magnitude for p in gas_profile.points],
        'pressure': [p.pressure.magnitude for p in gas_profile.points],
        'heat_transfer_rate': [p.heat_transfer_rate.magnitude if p.heat_transfer_rate else np.nan for p in gas_profile.points],
        'velocity': [p.velocity.magnitude if p.velocity else np.nan for p in gas_profile.points],
        'reynolds_number': [p.reynolds_number.magnitude if p.reynolds_number else np.nan for p in gas_profile.points],
        'density': [p.density.magnitude if p.density else np.nan for p in gas_profile.points],
    }
    water_data = {
        'z': [p.z.magnitude for p in water_profile.points],
        'temperature': [p.temperature.magnitude for p in water_profile.points],
        'enthalpy': [p.enthalpy.magnitude if p.enthalpy else np.nan for p in water_profile.points],
        'heat_transfer_rate': [p.heat_transfer_rate.magnitude if p.heat_transfer_rate else np.nan for p in water_profile.points],
        'quality': [p.quality.magnitude if p.quality else np.nan for p in water_profile.points],
        'saturation_temperature': [p.saturation_temperature.magnitude if p.saturation_temperature else np.nan for p in water_profile.points],
    }
    gas_df = pd.DataFrame(gas_data)
    water_df = pd.DataFrame(water_data)
    return gas_df, water_df

def compute_derived_quantities(gas_df: pd.DataFrame, water_df: pd.DataFrame) -> dict:
    """Perform additional calculations on results (add your custom vars here)."""
    derived = {}
    
    # Total heat transfer (integrate Q_dot over z)
    mask = ~gas_df['heat_transfer_rate'].isna()
    total_heat = np.trapz(gas_df['heat_transfer_rate'][mask], gas_df['z'][mask])  # W
    derived['total_heat_transfer'] = Q_(total_heat, 'watt')
    
    # Average gas velocity (ignore NaNs)
    avg_velocity = gas_df['velocity'].mean(skipna=True)
    derived['avg_gas_velocity'] = Q_(avg_velocity, 'm/s')
    
    # Phase change point (first z where quality > 0)
    boiling_start = water_df[water_df['quality'] > 0]['z'].min()
    derived['boiling_start_z'] = Q_(boiling_start, 'm') if not np.isnan(boiling_start) else None
    
    # Efficiency example (e.g., heat absorbed / potential; customize)
    potential_heat = (gas_df['temperature'].max() - gas_df['temperature'].min()) * 1e3  # Placeholder
    derived['efficiency'] = total_heat / potential_heat * 100  # %
    
    # Gradients (e.g., dT/dz for gas)
    gas_df['dT_dz'] = np.gradient(gas_df['temperature'], gas_df['z'])
    derived['max_gas_dT_dz'] = Q_(gas_df['dT_dz'].max(), 'K/m')
    
    return derived

def plot_profiles(gas_df: pd.DataFrame, water_df: pd.DataFrame, save_path: str = 'plots/'):
    """Generate graphs with multiple data series (subplots, overlays)."""
    plt.figure(figsize=(12, 8))
    
    # Subplot 1: Temperature profiles (overlay gas and water)
    fig, ax1 = plt.subplots()
    ax1.plot(gas_df['z'], gas_df['temperature'], label='Gas T', color='red')
    ax1.set_xlabel('z (m)')
    ax1.set_ylabel('Temperature (K)')
    ax2 = ax1.twinx()
    ax2.plot(water_df['z'], water_df['temperature'], label='Water T', color='blue', linestyle='--')
    ax2.set_ylabel('Temperature (K)')
    fig.legend()
    plt.title('Temperature Profiles')
    plt.savefig(f'{save_path}temperature_profiles.png')
    
    # Subplot 2: Heat transfer rate
    plt.figure()
    plt.plot(gas_df['z'], gas_df['heat_transfer_rate'], label='Q_dot')
    plt.xlabel('z (m)')
    plt.ylabel('Heat Transfer Rate (W/m)')
    plt.title('Heat Transfer Rate')
    plt.legend()
    plt.savefig(f'{save_path}heat_transfer.png')
    
    # More plots: Reynolds number distribution (histogram)
    plt.figure()
    sns.histplot(gas_df['reynolds_number'].dropna(), kde=True)
    plt.title('Reynolds Number Distribution')
    plt.savefig(f'{save_path}re_distribution.png')
    
    # Quality vs. z
    plt.figure()
    plt.plot(water_df['z'], water_df['quality'], label='Quality')
    plt.xlabel('z (m)')
    plt.ylabel('Vapor Quality')
    plt.title('Vapor Quality Profile')
    plt.legend()
    plt.savefig(f'{save_path}quality_profile.png')
    
    plt.close('all')  # Clean up

def generate_flow_diagram(gas_profile: GasStreamProfile, water_profile: WaterStreamProfile, save_path: str = 'plots/flow_diagram.png'):
    """Simple flow diagram (stages with key values; use for distributions too)."""
    graph = pydot.Dot(graph_type='digraph')
    
    # Nodes for stages (simplified; add from config if needed)
    stages = ['Inlet', 'Pass1', 'Superheater', 'Reversal1', 'Pass2', 'Reversal2', 'Pass3', 'Economizer', 'Outlet']
    prev_node = None
    for i, stage in enumerate(stages):
        label = f"{stage}\nAvg T_gas: {np.mean([p.temperature.magnitude for p in gas_profile.points[i*10:(i+1)*10]]):.1f} K"  # Example slice
        node = pydot.Node(stage, label=label)
        graph.add_node(node)
        if prev_node:
            graph.add_edge(pydot.Edge(prev_node, node))
        prev_node = node
    
    # Save as PNG
    graph.write_png(save_path)

def run_post_processing(gas_profile: GasStreamProfile, water_profile: WaterStreamProfile):
    """Main entry: Convert to DFs, compute derived, plot, generate diagrams."""
    gas_df, water_df = create_dataframes(gas_profile, water_profile)
    derived = compute_derived_quantities(gas_df, water_df)
    print("Derived Quantities:", derived)
    
    plot_profiles(gas_df, water_df)
    generate_flow_diagram(gas_profile, water_profile)