"""
Calculate volumes for positive and negative air spring chambers.

Reads configuration from inputs.yaml and calculates the annular volume
for each chamber based on shaft and spring diameters.
"""

import yaml
import math
import argparse

# Try to import matplotlib for graph generation
try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Constants
PSI_TO_PA = 6894.76  # Conversion factor: 1 psi = 6894.76 Pa
PSI_MM_TO_N_MM = 0.00689476  # Conversion factor: 1 psi × mm = 0.00689476 N/mm
MM2_TO_M2 = 1e-6  # Conversion factor: 1 mm² = 1e-6 m²
MM3_TO_LITERS = 1_000_000  # Conversion factor: 1 L = 1,000,000 mm³
MM3_TO_M3 = 1_000_000_000  # Conversion factor: 1 m³ = 1,000,000,000 mm³
SEPARATOR = "=" * 50  # Separator line for output formatting
IRT_LENGTH = 77  # IRT chamber length in mm
SHAFT_OUTER_DIAMETER = 10  # Shaft outer diameter in mm
SPRING_INNER_DIAMETER = 33  # Spring inner diameter in mm
POSITIVE_CHAMBER_LENGTH_OFFSET = 16  # Positive chamber length = travel + 16 mm
NEGATIVE_CHAMBER_LENGTH_OFFSET = 120  # Negative chamber length = (travel - 120) * 2 mm
NEGATIVE_CHAMBER_LENGTH_MULTIPLIER = 2  # Multiplier for negative chamber length calculation


def calculate_annular_volume(outer_diameter: float, inner_diameter: float, length: float) -> float:
    """
    Calculate the volume of an annular (ring-shaped) chamber.
    
    Volume = π × (R_outer² - R_inner²) × length
    
    Parameters:
    -----------
    outer_diameter : float
        Outer diameter of the annulus (mm)
    inner_diameter : float
        Inner diameter of the annulus (mm)
    length : float
        Length of the chamber (mm)
    
    Returns:
    --------
    float
        Volume in mm³
    """
    outer_radius = outer_diameter / 2.0
    inner_radius = inner_diameter / 2.0
    
    volume = math.pi * (outer_radius**2 - inner_radius**2) * length
    
    return volume


def load_inputs(file_path: str = "inputs.yaml") -> dict:
    """
    Load configuration from YAML file.
    
    Parameters:
    -----------
    file_path : str
        Path to the YAML input file
    
    Returns:
    --------
    dict
        Configuration dictionary
    
    Raises:
    -------
    FileNotFoundError
        If the input file doesn't exist
    yaml.YAMLError
        If the YAML file is invalid
    """
    try:
        with open(file_path, 'r') as f:
            inputs = yaml.safe_load(f)
        if inputs is None:
            raise ValueError(f"Input file '{file_path}' is empty or invalid")
        return inputs
    except FileNotFoundError:
        raise FileNotFoundError(f"Input file '{file_path}' not found")
    except yaml.YAMLError as e:
        raise yaml.YAMLError(f"Error parsing YAML file '{file_path}': {e}")


def calculate_annular_area(outer_diameter: float, inner_diameter: float) -> float:
    """
    Calculate the cross-sectional area of an annulus.
    
    Area = π × (R_outer² - R_inner²)
    
    Parameters:
    -----------
    outer_diameter : float
        Outer diameter (mm)
    inner_diameter : float
        Inner diameter (mm)
    
    Returns:
    --------
    float
        Area in mm²
    """
    outer_radius = outer_diameter / 2.0
    inner_radius = inner_diameter / 2.0
    
    area = math.pi * (outer_radius**2 - inner_radius**2)
    
    return area


def calculate_positive_chamber_pressure(
    initial_pressure: float,
    initial_volume: float,
    initial_length: float,
    travel: float,
    annular_area: float,
    irt_volume: float = 0,
    irt_pressure: float = 0
) -> dict:
    """
    Calculate positive chamber pressure at a given travel point using Boyle's Law.
    
    Boyle's Law: P₁V₁ = P₂V₂
    As the spring compresses, the positive chamber length decreases, reducing volume.
    If positive pressure >= IRT pressure, IRT volume is added to positive chamber.
    
    Parameters:
    -----------
    initial_pressure : float
        Initial pressure (psi)
    initial_volume : float
        Initial volume at travel = 0 (mm³)
    initial_length : float
        Initial length of positive chamber (mm)
    travel : float
        Travel distance from initial position (mm, 0 to max travel)
    annular_area : float
        Cross-sectional annular area (mm²)
    irt_volume : float, optional
        IRT chamber volume to add when pressure threshold is reached (mm³)
    irt_pressure : float, optional
        IRT pressure threshold (psi)
    
    Returns:
    --------
    dict
        Dictionary containing pressure, volume, length, and whether IRT is combined
    """
    # Calculate current length (compression reduces length)
    current_length = initial_length - travel
    
    # Ensure length doesn't go negative
    if current_length < 0:
        current_length = 0
    
    # Calculate current volume (positive chamber only)
    pos_volume_only = annular_area * current_length
    
    # Calculate pressure with positive chamber only first
    if pos_volume_only > 0:
        pos_pressure_only = (initial_pressure * initial_volume) / pos_volume_only
    else:
        pos_pressure_only = float('inf')
    
    # Check if positive pressure >= IRT pressure, then combine volumes
    irt_combined = False
    if irt_volume > 0 and irt_pressure > 0 and pos_pressure_only >= irt_pressure:
        # IRT volume is added to positive chamber
        # Find the travel point where pressure equals IRT pressure for combination
        # P_irt = (P₀ × V₀) / V_combine
        # V_combine = (P₀ × V₀) / P_irt
        combine_volume = (initial_pressure * initial_volume) / irt_pressure
        combine_travel = initial_length - (combine_volume / annular_area)
        
        # If we're at or past the combination point
        if travel >= combine_travel:
            irt_combined = True
            current_volume = pos_volume_only + irt_volume
            
            # After combination, use Boyle's Law with combined system
            # Initial state at combination: P = IRT pressure, V = pos_volume_at_combine + irt_volume
            combine_pos_volume = combine_volume
            combined_initial_volume = combine_pos_volume + irt_volume
            combined_initial_pressure = irt_pressure
            
            # Current state: volume is current_volume, calculate pressure
            if current_volume > 0:
                current_pressure = (combined_initial_pressure * combined_initial_volume) / current_volume
            else:
                current_pressure = float('inf')
        else:
            # Not yet at combination point
            current_volume = pos_volume_only
            current_pressure = pos_pressure_only
    else:
        # IRT not combined, use positive chamber only
        current_volume = pos_volume_only
        current_pressure = pos_pressure_only
    
    return {
        'travel': travel,
        'pressure_psi': current_pressure,
        'volume_mm3': current_volume,
        'volume_liters': current_volume / MM3_TO_LITERS,
        'length': current_length,
        'irt_combined': irt_combined
    }


def calculate_negative_chamber_pressure(
    initial_pressure: float,
    initial_volume: float,
    initial_length: float,
    travel: float,
    annular_area: float
) -> dict:
    """
    Calculate negative chamber pressure at a given travel point using Boyle's Law.
    
    As the spring compresses, the negative chamber typically expands (length increases),
    increasing volume and decreasing pressure.
    
    Parameters:
    -----------
    initial_pressure : float
        Initial pressure (psi)
    initial_volume : float
        Initial volume at travel = 0 (mm³)
    initial_length : float
        Initial length of negative chamber (mm)
    travel : float
        Travel distance from initial position (mm, 0 to max travel)
    annular_area : float
        Cross-sectional annular area (mm²)
    
    Returns:
    --------
    dict
        Dictionary containing pressure, volume, and length at travel point
    """
    # Calculate current length (negative chamber expands as positive compresses)
    current_length = initial_length + travel
    
    # Calculate current volume
    current_volume = annular_area * current_length
    
    # Use Boyle's Law: P₁V₁ = P₂V₂
    # P₂ = P₁V₁ / V₂
    if current_volume > 0:
        current_pressure = (initial_pressure * initial_volume) / current_volume
    else:
        current_pressure = 0
    
    return {
        'travel': travel,
        'pressure_psi': current_pressure,
        'volume_mm3': current_volume,
        'volume_liters': current_volume / MM3_TO_LITERS,
        'length': current_length
    }


def calculate_irt_chamber_pressure(
    initial_pressure: float,
    initial_volume: float,
    initial_length: float,
    travel: float,
    annular_area: float
) -> dict:
    """
    Calculate IRT chamber pressure at a given travel point using Boyle's Law.
    
    The IRT chamber typically maintains a constant volume (length doesn't change with travel),
    so pressure remains constant.
    
    Parameters:
    -----------
    initial_pressure : float
        Initial pressure (psi)
    initial_volume : float
        Initial volume at travel = 0 (mm³)
    initial_length : float
        Initial length of IRT chamber (mm)
    travel : float
        Travel distance from initial position (mm, 0 to max travel)
    annular_area : float
        Cross-sectional annular area (mm²)
    
    Returns:
    --------
    dict
        Dictionary containing pressure, volume, and length at travel point
    """
    # IRT chamber maintains constant volume (length doesn't change)
    current_length = initial_length
    current_volume = initial_volume
    
    # Pressure remains constant since volume doesn't change
    current_pressure = initial_pressure
    
    return {
        'travel': travel,
        'pressure_psi': current_pressure,
        'volume_mm3': current_volume,
        'volume_liters': current_volume / MM3_TO_LITERS,
        'length': current_length
    }


def validate_travel(travel: float, max_travel: float) -> float:
    """
    Validate and clamp travel value to valid range.
    
    Parameters:
    -----------
    travel : float
        Travel distance in mm
    max_travel : float
        Maximum travel distance in mm
    
    Returns:
    --------
    float
        Validated travel value
    """
    if travel < 0:
        return 0.0
    elif travel > max_travel:
        return max_travel
    return travel


def calculate_spring_rate(
    travel: float,
    results: dict = None,
    inputs: dict = None
) -> dict:
    """
    Calculate spring rate at a given travel point considering both positive and negative chambers.
    
    Spring rate: k = dF/dx
    Total force: F = F_positive - F_negative = (P_pos - P_neg) × A
    Spring rate: k = A × (dP_pos/dx - dP_neg/dx)
    
    Using Boyle's Law derivatives:
    dP/dx = -P₀V₀/(V²) × (dV/dx)
    
    Parameters:
    -----------
    travel : float
        Travel distance in mm (0 to max travel)
    inputs : dict, optional
        Input dictionary. If None, loads from inputs.yaml
    
    Returns:
    --------
    dict
        Dictionary with spring rate in N/mm, forces in N, and pressures in psi
    """
    if results is None:
        if inputs is None:
            inputs = load_inputs()
        results = calculate_chamber_volumes(inputs)
    
    # Validate travel range
    max_travel = results['travel']
    travel = validate_travel(travel, max_travel)
    
    # Get chamber data
    annular_area = results['annular_area']  # mm²
    
    # Positive chamber
    pos_initial_pressure = results['positive_chamber']['starting_pressure']  # psi
    pos_initial_volume = results['positive_chamber']['initial_volume_mm3']  # mm³
    pos_initial_length = results['positive_chamber']['length']  # mm
    
    # Negative chamber
    neg_initial_pressure = results['negative_chamber']['starting_pressure']  # psi
    neg_initial_volume = results['negative_chamber']['volume_mm3']  # mm³
    neg_initial_length = results['negative_chamber']['length']  # mm
    
    # IRT chamber (if present)
    has_irt = 'irt_chamber' in results
    if has_irt:
        irt_initial_pressure = results['irt_chamber']['starting_pressure']  # psi
        irt_initial_volume = results['irt_chamber']['volume_mm3']  # mm³
        irt_initial_length = results['irt_chamber']['length']  # mm
    
    # Calculate current state
    # Pass IRT info to positive chamber calculation for conditional combining
    irt_vol = irt_initial_volume if has_irt else 0
    irt_pres = irt_initial_pressure if has_irt else 0
    
    pos_data = calculate_positive_chamber_pressure(
        pos_initial_pressure, pos_initial_volume, pos_initial_length,
        travel, annular_area, irt_volume=irt_vol, irt_pressure=irt_pres
    )
    neg_data = calculate_negative_chamber_pressure(
        neg_initial_pressure, neg_initial_volume, neg_initial_length,
        travel, annular_area
    )
    
    # IRT chamber (only if not combined with positive)
    irt_combined = pos_data.get('irt_combined', False)
    if has_irt and not irt_combined:
        irt_data = calculate_irt_chamber_pressure(
            irt_initial_pressure, irt_initial_volume, irt_initial_length,
            travel, annular_area
        )
    elif has_irt and irt_combined:
        # IRT is combined, so it has no separate pressure/volume
        irt_data = {
            'pressure_psi': pos_data['pressure_psi'],  # Same as combined pressure
            'volume_mm3': 0,  # Volume is now in positive chamber
            'volume_liters': 0,
            'length': irt_initial_length
        }
    
    # Current volumes
    pos_volume = pos_data['volume_mm3']
    neg_volume = neg_data['volume_mm3']
    
    # Volume change rates
    # Positive chamber: compresses, dV/dx = -A (negative)
    # If IRT is combined, volume change is still -A (IRT volume is constant)
    dV_pos_dx = -annular_area  # mm³/mm
    
    # Negative chamber: expands, dV/dx = +A (positive)
    dV_neg_dx = annular_area  # mm³/mm
    
    # Pressure change rates (Boyle's Law derivative) in psi/mm
    # dP/dx = -P₀V₀/(V²) × (dV/dx)
    # If IRT is combined, use combined initial volume and current volume
    if irt_combined:
        # Combined initial: pos_initial_volume + irt_initial_volume
        # Combined current: pos_volume
        combined_initial_volume = pos_initial_volume + irt_initial_volume
        # Initial pressure is the pressure when they combine (IRT pressure)
        combined_initial_pressure = irt_initial_pressure
        if pos_volume > 0:
            dP_pos_dx = -(combined_initial_pressure * combined_initial_volume) / (pos_volume ** 2) * dV_pos_dx
        else:
            dP_pos_dx = float('inf')
    else:
        # Not combined, use positive chamber only
        if pos_volume > 0:
            dP_pos_dx = -(pos_initial_pressure * pos_initial_volume) / (pos_volume ** 2) * dV_pos_dx
        else:
            dP_pos_dx = float('inf')
    
    if neg_volume > 0:
        dP_neg_dx = -(neg_initial_pressure * neg_initial_volume) / (neg_volume ** 2) * dV_neg_dx
    else:
        dP_neg_dx = 0
    
    # Spring rate: k = A × (dP_pos/dx - dP_neg/dx)
    # Calculate directly in N/mm
    # dP/dx is in psi/mm, A is in mm²
    # Result: (mm²) × (psi/mm) = psi × mm
    
    # To convert psi × mm to N/mm:
    # 1 psi = 6894.76 Pa = 6894.76 N/m²
    # 1 mm = 0.001 m
    # 1 psi × mm = 6894.76 N/m² × 0.001 m = 6.89476 N/m
    # 1 N/m = 0.001 N/mm
    # So: 1 psi × mm = 6.89476 N/m = 0.00689476 N/mm
    # Conversion factor: 0.00689476 N/mm per (psi × mm)
    
    # Spring rate in N/mm
    spring_rate_n_per_mm = annular_area * (dP_pos_dx - dP_neg_dx) * PSI_MM_TO_N_MM
    
    # Also calculate forces in N (keeping pressures in psi)
    # Convert area from mm² to m²
    area_m2 = annular_area * MM2_TO_M2
    # F = P × A, where P is in Pa and A is in m², result is in N
    force_positive_n = pos_data['pressure_psi'] * PSI_TO_PA * area_m2
    force_negative_n = neg_data['pressure_psi'] * PSI_TO_PA * area_m2
    force_irt_n = 0
    if has_irt:
        force_irt_n = irt_data['pressure_psi'] * PSI_TO_PA * area_m2
    # Net force: positive - negative (IRT typically doesn't contribute to net force if constant volume)
    force_net_n = force_positive_n - force_negative_n
    
    # Convert spring rate to N/m for reference (1 N/mm = 1000 N/m)
    spring_rate_n_per_m = spring_rate_n_per_mm * 1000
    
    result = {
        'travel': travel,
        'spring_rate_n_per_mm': spring_rate_n_per_mm,
        'spring_rate_n_per_m': spring_rate_n_per_m,
        'force_positive_n': force_positive_n,
        'force_negative_n': force_negative_n,
        'force_net_n': force_net_n,
        'pressure_positive_psi': pos_data['pressure_psi'],
        'pressure_negative_psi': neg_data['pressure_psi'],
        'volume_positive_mm3': pos_volume,
        'volume_negative_mm3': neg_volume,
        'dP_pos_dx_psi_per_mm': dP_pos_dx,
        'dP_neg_dx_psi_per_mm': dP_neg_dx
    }
    
    if has_irt:
        result['force_irt_n'] = force_irt_n
        result['pressure_irt_psi'] = irt_data['pressure_psi']
        result['volume_irt_mm3'] = irt_data['volume_mm3']
    
    return result


def calculate_chamber_volumes(inputs: dict) -> dict:
    """
    Calculate volumes for positive and negative chambers.
    
    Parameters:
    -----------
    inputs : dict
        Configuration dictionary from inputs.yaml
    
    Returns:
    --------
    dict
        Dictionary containing calculated volumes and input parameters
    """
    shaft_outer_diameter = SHAFT_OUTER_DIAMETER  # mm
    spring_inner_diameter = SPRING_INNER_DIAMETER  # mm
    starting_pressure = inputs['startingPressure']  # psi
    travel = inputs.get('travel', 0)  # mm
    positive_length = travel + POSITIVE_CHAMBER_LENGTH_OFFSET  # mm
    negative_length = (travel - NEGATIVE_CHAMBER_LENGTH_OFFSET) * NEGATIVE_CHAMBER_LENGTH_MULTIPLIER  # mm
    irt_pressure = inputs.get('irtPressure', starting_pressure)  # psi
    irt_length = IRT_LENGTH  # mm
    
    # Calculate annular cross-sectional area
    annular_area = calculate_annular_area(spring_inner_diameter, shaft_outer_diameter)
    
    # Calculate volumes
    positive_volume = calculate_annular_volume(
        spring_inner_diameter,
        shaft_outer_diameter,
        positive_length
    )
    
    negative_volume = calculate_annular_volume(
        spring_inner_diameter,
        shaft_outer_diameter,
        negative_length
    )
    
    # Calculate IRT volume (always calculated since IRT_LENGTH is constant and > 0)
    irt_volume = calculate_annular_volume(
        spring_inner_diameter,
        shaft_outer_diameter,
        irt_length
    )
    
    result = {
        'shaft_outer_diameter': shaft_outer_diameter,
        'spring_inner_diameter': spring_inner_diameter,
        'starting_pressure': starting_pressure,
        'travel': travel,
        'annular_area': annular_area,
        'positive_chamber': {
            'length': positive_length,
            'starting_pressure': starting_pressure,
            'initial_volume_mm3': positive_volume,
            'volume_mm3': positive_volume,
            'volume_liters': positive_volume / MM3_TO_LITERS,
            'volume_m3': positive_volume / MM3_TO_M3
        },
        'negative_chamber': {
            'length': negative_length,
            'starting_pressure': starting_pressure,
            'volume_mm3': negative_volume,
            'volume_liters': negative_volume / MM3_TO_LITERS,
            'volume_m3': negative_volume / MM3_TO_M3
        }
    }
    
    # Add IRT chamber (always present since IRT_LENGTH is constant and > 0)
    result['irt_chamber'] = {
        'length': irt_length,
        'starting_pressure': irt_pressure,
        'volume_mm3': irt_volume,
        'volume_liters': irt_volume / MM3_TO_LITERS,
        'volume_m3': irt_volume / MM3_TO_M3
    }
    
    return result


def create_force_graph(results: dict, inputs: dict, output_file: str = 'force_vs_travel.png'):
    """
    Create a graph showing force vs travel.
    
    Parameters:
    -----------
    results : dict
        Dictionary containing calculated chamber volumes and parameters
    inputs : dict
        Input dictionary from inputs.yaml
    output_file : str
        Output filename for the graph image (default: 'force_vs_travel.png')
    
    Raises:
    -------
    ImportError
        If matplotlib is not installed
    """
    if not HAS_MATPLOTLIB:
        raise ImportError(
            "matplotlib is required for graph generation. "
            "Install it with: pip install matplotlib"
        )
    max_travel = results['travel']
    # Generate travel points at 1mm increments for smoother curve
    travel_points = [x for x in range(0, int(max_travel) + 1, 1)]
    
    # Calculate force at each travel point
    forces = []
    for travel_point in travel_points:
        spring_data = calculate_spring_rate(travel_point, results, inputs)
        forces.append(spring_data['force_net_n'])
    
    # Get input values for display
    max_travel_input = inputs.get('travel', max_travel)
    starting_pressure = inputs.get('startingPressure', 0)
    irt_pressure = inputs.get('irtPressure', 0)
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(travel_points, forces, 'b-', linewidth=2)
    plt.xlabel('Travel (mm)', fontsize=12)
    plt.ylabel('Force (N)', fontsize=12)
    plt.title('Force vs Travel', fontsize=14, fontweight='bold')
    plt.grid(True, alpha=0.3)
    
    # Add input values as text box
    input_text = f'Travel: {max_travel_input} mm\n'
    input_text += f'Main Chamber Pressure: {starting_pressure} psi\n'
    input_text += f'IRT Pressure: {irt_pressure} psi'
    
    # Position text box in upper left corner
    plt.text(0.02, 0.98, input_text,
             transform=plt.gca().transAxes,
             fontsize=10,
             verticalalignment='top',
             horizontalalignment='left',
             bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8))
    
    plt.tight_layout()
    
    # Save the figure
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Graph saved to: {output_file}")


def main():
    """Main function to load inputs and calculate volumes."""
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description='Calculate air spring chamber volumes and spring rates')
    parser.add_argument('-v', '--verbose', action='store_true', 
                       help='Show detailed output including configuration, volumes, and pressures')
    parser.add_argument('-g', '--graph', nargs='?', const='force_vs_travel.png',
                       help='Generate a graph image of force vs travel (optionally specify output filename)')
    args = parser.parse_args()
    
    # Load inputs
    inputs = load_inputs()
    
    # Calculate volumes
    results = calculate_chamber_volumes(inputs)
    
    # Generate travel points at 10mm increments from 0 to max_travel
    max_travel = results['travel']
    travel_points = [x for x in range(0, int(max_travel) + 1, 10)]
    # Ensure max_travel is included even if not on a 10mm increment
    if travel_points[-1] != max_travel:
        travel_points.append(max_travel)
    
    if args.verbose:
        # Verbose output: show all details
        print("Air Spring Chamber Volume Calculator")
        print(SEPARATOR)
        print(f"\nConfiguration:")
        print(f"  Shaft outer diameter: {results['shaft_outer_diameter']} mm")
        print(f"  Spring inner diameter: {results['spring_inner_diameter']} mm")
        print(f"  Starting pressure: {results['starting_pressure']} psi")
        print(f"  Max travel: {results['travel']} mm")
        
        print(f"\nPositive Chamber (Initial):")
        print(f"  Length: {results['positive_chamber']['length']} mm")
        print(f"  Starting pressure: {results['positive_chamber']['starting_pressure']} psi")
        print(f"  Volume: {results['positive_chamber']['volume_mm3']:.2f} mm³")
        print(f"  Volume: {results['positive_chamber']['volume_liters']:.6f} L")
        print(f"  Volume: {results['positive_chamber']['volume_m3']:.9f} m³")
        
        print(f"\nNegative Chamber:")
        print(f"  Length: {results['negative_chamber']['length']} mm")
        print(f"  Starting pressure: {results['negative_chamber']['starting_pressure']} psi")
        print(f"  Volume: {results['negative_chamber']['volume_mm3']:.2f} mm³")
        print(f"  Volume: {results['negative_chamber']['volume_liters']:.6f} L")
        print(f"  Volume: {results['negative_chamber']['volume_m3']:.9f} m³")
        
        if 'irt_chamber' in results:
            print(f"\nIRT Chamber:")
            print(f"  Length: {results['irt_chamber']['length']} mm")
            print(f"  Starting pressure: {results['irt_chamber']['starting_pressure']} psi")
            print(f"  Volume: {results['irt_chamber']['volume_mm3']:.2f} mm³")
            print(f"  Volume: {results['irt_chamber']['volume_liters']:.6f} L")
            print(f"  Volume: {results['irt_chamber']['volume_m3']:.9f} m³")
        
        print(f"\nTotal Volume:")
        total_volume_mm3 = results['positive_chamber']['volume_mm3'] + results['negative_chamber']['volume_mm3']
        total_volume_liters = results['positive_chamber']['volume_liters'] + results['negative_chamber']['volume_liters']
        total_volume_m3 = results['positive_chamber']['volume_m3'] + results['negative_chamber']['volume_m3']
        if 'irt_chamber' in results:
            total_volume_mm3 += results['irt_chamber']['volume_mm3']
            total_volume_liters += results['irt_chamber']['volume_liters']
            total_volume_m3 += results['irt_chamber']['volume_m3']
        print(f"  Total: {total_volume_mm3:.2f} mm³")
        print(f"  Total: {total_volume_liters:.6f} L")
        print(f"  Total: {total_volume_m3:.9f} m³")
        
        # Calculate and display pressure at different travel points
        print(f"\n{SEPARATOR}")
        print("Chamber Pressures vs Travel (Boyle's Law)")
        print(SEPARATOR)
        
        has_irt = 'irt_chamber' in results
        irt_vol = results['irt_chamber']['volume_mm3'] if has_irt else 0
        irt_pres = results['irt_chamber']['starting_pressure'] if has_irt else 0
        
        for travel_point in travel_points:
            pos_data = calculate_positive_chamber_pressure(
                initial_pressure=results['positive_chamber']['starting_pressure'],
                initial_volume=results['positive_chamber']['initial_volume_mm3'],
                initial_length=results['positive_chamber']['length'],
                travel=travel_point,
                annular_area=results['annular_area'],
                irt_volume=irt_vol,
                irt_pressure=irt_pres
            )
            
            neg_data = calculate_negative_chamber_pressure(
                initial_pressure=results['negative_chamber']['starting_pressure'],
                initial_volume=results['negative_chamber']['volume_mm3'],
                initial_length=results['negative_chamber']['length'],
                travel=travel_point,
                annular_area=results['annular_area']
            )
            
            irt_combined = pos_data.get('irt_combined', False)
            if has_irt and not irt_combined:
                irt_data = calculate_irt_chamber_pressure(
                    initial_pressure=results['irt_chamber']['starting_pressure'],
                    initial_volume=results['irt_chamber']['volume_mm3'],
                    initial_length=results['irt_chamber']['length'],
                    travel=travel_point,
                    annular_area=results['annular_area']
                )
            elif has_irt and irt_combined:
                irt_data = {
                    'pressure_psi': pos_data['pressure_psi'],
                    'volume_mm3': 0,
                    'volume_liters': 0,
                    'length': results['irt_chamber']['length']
                }
            
            percent_travel = (travel_point / max_travel) * 100
            print(f"\nTravel: {travel_point:.1f} mm ({percent_travel:.1f}% through travel)")
            print(f"  Positive Chamber:")
            if irt_combined:
                print(f"    Pressure: {pos_data['pressure_psi']:.2f} psi (IRT combined)")
                print(f"    Volume: {pos_data['volume_mm3']:.2f} mm³ ({pos_data['volume_liters']:.6f} L) [includes IRT]")
            else:
                print(f"    Pressure: {pos_data['pressure_psi']:.2f} psi")
                print(f"    Volume: {pos_data['volume_mm3']:.2f} mm³ ({pos_data['volume_liters']:.6f} L)")
            print(f"    Length: {pos_data['length']:.2f} mm")
            print(f"  Negative Chamber:")
            print(f"    Pressure: {neg_data['pressure_psi']:.2f} psi")
            print(f"    Volume: {neg_data['volume_mm3']:.2f} mm³ ({neg_data['volume_liters']:.6f} L)")
            print(f"    Length: {neg_data['length']:.2f} mm")
            if has_irt:
                if irt_combined:
                    print(f"  IRT Chamber: Combined with Positive Chamber")
                else:
                    print(f"  IRT Chamber:")
                    print(f"    Pressure: {irt_data['pressure_psi']:.2f} psi")
                    print(f"    Volume: {irt_data['volume_mm3']:.2f} mm³ ({irt_data['volume_liters']:.6f} L)")
                    print(f"    Length: {irt_data['length']:.2f} mm")
        
        # Calculate and display spring rate at different travel points
        print(f"\n{SEPARATOR}")
        print("Spring Rate vs Travel (Both Chambers)")
        print(SEPARATOR)
        
        for travel_point in travel_points:
            spring_data = calculate_spring_rate(travel_point, results, inputs)
            percent_travel = (travel_point / max_travel) * 100
            
            print(f"\nTravel: {travel_point:.1f} mm ({percent_travel:.1f}% through travel)")
            print(f"  Spring Rate: {spring_data['spring_rate_n_per_mm']:.2f} N/mm ({spring_data['spring_rate_n_per_m']:.2f} N/m)")
            print(f"  Net Force: {spring_data['force_net_n']:.2f} N")
            print(f"  Positive Force: {spring_data['force_positive_n']:.2f} N (from {spring_data['pressure_positive_psi']:.2f} psi)")
            print(f"  Negative Force: {spring_data['force_negative_n']:.2f} N (from {spring_data['pressure_negative_psi']:.2f} psi)")
            if 'force_irt_n' in spring_data:
                # Get pos_data for this travel point to check IRT combined status
                pos_data = calculate_positive_chamber_pressure(
                    initial_pressure=results['positive_chamber']['starting_pressure'],
                    initial_volume=results['positive_chamber']['initial_volume_mm3'],
                    initial_length=results['positive_chamber']['length'],
                    travel=travel_point,
                    annular_area=results['annular_area'],
                    irt_volume=irt_vol,
                    irt_pressure=irt_pres
                )
                irt_combined_status = " (combined)" if pos_data.get('irt_combined', False) else ""
                print(f"  IRT Force: {spring_data['force_irt_n']:.2f} N (from {spring_data['pressure_irt_psi']:.2f} psi){irt_combined_status}")
            print(f"  Pressure Difference: {spring_data['pressure_positive_psi'] - spring_data['pressure_negative_psi']:.2f} psi")
        
        # Calculate travel point where positive chamber pressure equals irtPressure
        if 'irtPressure' in inputs:
            irt_pressure = inputs['irtPressure']
            print(f"\n{SEPARATOR}")
            print(f"Travel Point at IRT Pressure ({irt_pressure} psi)")
            print(SEPARATOR)
            
            irt_data = find_travel_at_pressure(irt_pressure, results, inputs)
            percent_travel_irt = (irt_data['travel'] / max_travel) * 100
            
            print(f"\nTravel: {irt_data['travel']:.2f} mm ({percent_travel_irt:.1f}% through travel)")
    else:
        # Simplified output: only spring rate
        print("Spring Rate vs Travel")
        print(SEPARATOR)
        
        for travel_point in travel_points:
            spring_data = calculate_spring_rate(travel_point, results, inputs)
            percent_travel = (travel_point / max_travel) * 100
            print(f"{travel_point:.1f} mm ({percent_travel:.1f}%): {spring_data['spring_rate_n_per_mm']:.2f} N/mm")
    
    # Generate graph if requested
    if args.graph:
        create_force_graph(results, inputs, args.graph)
    
    return results


def find_travel_at_pressure(target_pressure: float, results: dict = None, inputs: dict = None) -> dict:
    """
    Find the travel point where positive chamber pressure equals target pressure.
    
    Uses inverse of Boyle's Law to solve for travel:
    P₂ = P₁V₁ / V₂
    Where V₂ = A × (L₀ - travel)
    Solving for travel: travel = L₀ - (P₁V₁) / (P₂A)
    
    Parameters:
    -----------
    target_pressure : float
        Target pressure in psi
    results : dict, optional
        Pre-calculated chamber volumes dictionary
    inputs : dict, optional
        Input dictionary. If None, loads from inputs.yaml
    
    Returns:
    --------
    dict
        Dictionary with travel point, pressure, volume, and length
    """
    if results is None:
        if inputs is None:
            inputs = load_inputs()
        results = calculate_chamber_volumes(inputs)
    
    # Get initial values
    initial_pressure = results['positive_chamber']['starting_pressure']  # psi
    initial_volume = results['positive_chamber']['initial_volume_mm3']  # mm³
    initial_length = results['positive_chamber']['length']  # mm
    annular_area = results['annular_area']  # mm²
    max_travel = results['travel']  # mm
    
    # Solve for travel using inverse of Boyle's Law
    # P₂ = P₁V₁ / (A × (L₀ - travel))
    # travel = L₀ - (P₁V₁) / (P₂A)
    travel = initial_length - (initial_pressure * initial_volume) / (target_pressure * annular_area)
    
    # Validate travel range
    travel = validate_travel(travel, max_travel)
    
    # Calculate actual pressure at this travel point to verify
    pressure_data = calculate_positive_chamber_pressure(
        initial_pressure=initial_pressure,
        initial_volume=initial_volume,
        initial_length=initial_length,
        travel=travel,
        annular_area=annular_area
    )
    
    return {
        'travel': travel,
        'target_pressure_psi': target_pressure,
        'actual_pressure_psi': pressure_data['pressure_psi'],
        'volume_mm3': pressure_data['volume_mm3'],
        'volume_liters': pressure_data['volume_liters'],
        'length': pressure_data['length']
    }


def get_pressure_at_travel(travel: float, results: dict = None, inputs: dict = None) -> dict:
    """
    Calculate positive chamber pressure at a specific travel point.
    
    Parameters:
    -----------
    travel : float
        Travel distance in mm (0 to max travel)
    results : dict, optional
        Pre-calculated chamber volumes dictionary
    inputs : dict, optional
        Input dictionary. If None, loads from inputs.yaml
    
    Returns:
    --------
    dict
        Dictionary with pressure, volume, and length at travel point
    """
    if results is None:
        if inputs is None:
            inputs = load_inputs()
        results = calculate_chamber_volumes(inputs)
    
    # Validate travel range
    max_travel = results['travel']
    travel = validate_travel(travel, max_travel)
    
    return calculate_positive_chamber_pressure(
        initial_pressure=results['positive_chamber']['starting_pressure'],
        initial_volume=results['positive_chamber']['initial_volume_mm3'],
        initial_length=results['positive_chamber']['length'],
        travel=travel,
        annular_area=results['annular_area']
    )


if __name__ == "__main__":
    results = main()

