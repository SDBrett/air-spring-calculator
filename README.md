# Air Spring Calculator

This is based on the Manitou Mezzer Pro fork, can be used for other Manitou IRT forks but the force numbers will be incorrect.

A Python tool for calculating air spring chamber volumes, pressures, and spring rates using Boyle's Law for positive, negative, and IRT (Infinite Rate Tune) chambers.

## Overview

This project calculates volumes, pressures, and spring rates for air spring chambers using annular (ring-shaped) geometry. The calculations account for:
- Positive chamber compression (volume decreases with travel)
- Negative chamber expansion (volume increases with travel)
- IRT chamber that combines with positive chamber when pressure threshold is reached
- Spring rate calculations considering both chambers

## Features

- Calculates volumes, pressures, and spring rates at different travel points
- Supports positive, negative, and IRT (Infinite Rate Tune) chambers
- Uses Boyle's Law (P₁V₁ = P₂V₂) for pressure calculations
- Displays results in multiple units (mm³, liters, m³, psi, N/mm, N)
- Simplified output (default) or verbose detailed output
- Graph generation for force vs travel visualization
- Reads configuration from YAML input file

## Requirements

- Python 3.6+
- PyYAML
- matplotlib (optional, for graph generation)

Install dependencies:
```bash
pip install -r requirements.txt
```

Or install individually:
```bash
pip install pyyaml matplotlib
```

## Input File Format

The `inputs.yaml` file contains the configuration parameters:

```yaml
travel: 160              # mm - Maximum travel distance
startingPressure: 56    # psi - Initial pressure for all chambers
irtPressure: 86         # psi - IRT chamber pressure threshold
```

### Parameters

- **travel**: Maximum travel distance in millimeters (mm)
- **startingPressure**: Initial pressure in psi for positive and negative chambers
- **irtPressure**: IRT chamber pressure threshold in psi (when reached, IRT volume combines with positive chamber)

### Constants (Defined in Code)

The following values are defined as constants in the code and can be modified there:
- **SHAFT_OUTER_DIAMETER**: 10 mm
- **SPRING_INNER_DIAMETER**: 33 mm
- **IRT_LENGTH**: 77 mm
- **POSITIVE_CHAMBER_LENGTH_OFFSET**: 16 mm (positive length = travel + 16)
- **NEGATIVE_CHAMBER_LENGTH_OFFSET**: 120 mm (negative length = (travel - 120) × 2)
- **NEGATIVE_CHAMBER_LENGTH_MULTIPLIER**: 2

## Chamber Length Calculations

- **Positive Chamber**: `length = travel + 16` mm
  - As travel increases, positive chamber compresses (length decreases)
- **Negative Chamber**: `length = (travel - 120) × 2` mm
  - As travel increases, negative chamber expands (length increases)
- **IRT Chamber**: Constant length of 77 mm
  - Combines with positive chamber when positive pressure ≥ IRT pressure

## Volume Calculation Formula

The volume of each annular chamber is calculated using:

```
Volume = π × (R_outer² - R_inner²) × length
```

Where:
- **R_outer** = SPRING_INNER_DIAMETER / 2
- **R_inner** = SHAFT_OUTER_DIAMETER / 2
- **length** = chamber length (varies with travel for positive/negative chambers)

## Pressure Calculation (Boyle's Law)

Boyle's Law states: **P₁V₁ = P₂V₂** (at constant temperature)

As the spring compresses:
- Positive chamber volume decreases → pressure increases
- Negative chamber volume increases → pressure decreases
- Net force = (P_positive - P_negative) × Area

## Spring Rate Calculation

Spring rate (k) is calculated as the derivative of force with respect to travel:

```
k = dF/dx = A × (dP_pos/dx - dP_neg/dx)
```

Where:
- **A** = annular cross-sectional area
- **dP/dx** = pressure change rate (from Boyle's Law derivative)

## Usage

### Basic Usage (Simplified Output)

```bash
python calculate_volumes.py
```

Output shows spring rate at 10mm travel increments:
```
Spring Rate vs Travel
==================================================
0.0 mm (0.0%): 5.45 N/mm
10.0 mm (6.2%): 4.88 N/mm
...
```

### Verbose Output

```bash
python calculate_volumes.py --verbose
```

Shows detailed information including:
- Configuration parameters
- Initial chamber volumes
- Pressures at each travel point
- Spring rates and forces
- IRT combination status

### Generate Graph

```bash
python calculate_volumes.py --graph
```

Creates `force_vs_travel.png` showing force vs travel curve.

```bash
python calculate_volumes.py --graph my_graph.png
```

Creates a custom filename.

### Combine Options

```bash
python calculate_volumes.py --verbose --graph output.png
```

## Command-Line Arguments

- `-v, --verbose`: Show detailed output including configuration, volumes, and pressures
- `-g, --graph [FILENAME]`: Generate a graph image of force vs travel (default: `force_vs_travel.png`)
- `-h, --help`: Show help message

## Example Output

### Simplified Output
```
Spring Rate vs Travel
==================================================
0.0 mm (0.0%): 5.45 N/mm
10.0 mm (6.2%): 4.88 N/mm
20.0 mm (12.5%): 4.57 N/mm
...
160.0 mm (100.0%): 13.89 N/mm
```

### Verbose Output (excerpt)
```
Air Spring Chamber Volume Calculator
==================================================

Configuration:
  Shaft outer diameter: 10 mm
  Spring inner diameter: 33 mm
  Starting pressure: 56 psi
  Max travel: 160 mm

Positive Chamber (Initial):
  Length: 176 mm
  Starting pressure: 56 psi
  Volume: 136709.55 mm³
  Volume: 0.136710 L
  Volume: 0.000136710 m³

Negative Chamber:
  Length: 80 mm
  Starting pressure: 56 psi
  Volume: 62140.70 mm³
  Volume: 0.062141 L
  Volume: 0.000062141 m³

IRT Chamber:
  Length: 77 mm
  Starting pressure: 86 psi
  Volume: 59810.43 mm³
  Volume: 0.059810 L
  Volume: 0.000059810 m³

Chamber Pressures vs Travel (Boyle's Law)
==================================================

Travel: 0.0 mm (0.0% through travel)
  Positive Chamber:
    Pressure: 56.00 psi
    Volume: 136709.55 mm³ (0.136710 L)
    Length: 176.00 mm
  Negative Chamber:
    Pressure: 56.00 psi
    Volume: 62140.70 mm³ (0.062141 L)
    Length: 80.00 mm
  IRT Chamber:
    Pressure: 86.00 psi
    Volume: 59810.43 mm³ (0.059810 L)
    Length: 77.00 mm
...
```

## Project Structure

```
air-spring-calc/
├── calculate_volumes.py  # Main calculation script
├── inputs.yaml            # Configuration file
├── requirements.txt       # Python dependencies
└── README.md             # This file
```

## Key Functions

### `calculate_annular_volume(outer_diameter, inner_diameter, length)`
Calculates the volume of an annular chamber.

### `calculate_annular_area(outer_diameter, inner_diameter)`
Calculates the cross-sectional area of an annulus.

### `calculate_positive_chamber_pressure(...)`
Calculates positive chamber pressure at a given travel point using Boyle's Law. Handles IRT combination when pressure threshold is reached.

### `calculate_negative_chamber_pressure(...)`
Calculates negative chamber pressure at a given travel point using Boyle's Law.

### `calculate_spring_rate(travel, results, inputs)`
Calculates spring rate considering both positive and negative chambers, including IRT effects.

### `calculate_chamber_volumes(inputs)`
Calculates initial volumes for all chambers based on configuration.

### `create_force_graph(results, inputs, output_file)`
Generates a graph showing force vs travel.

## IRT Chamber Behavior

The IRT (Infinite Rate Tune) chamber has special behavior:
- Maintains constant volume (77 mm length)
- Has its own pressure (86 psi by default)
- When positive chamber pressure reaches or exceeds IRT pressure, the IRT volume is added to the positive chamber
- After combination, the combined system follows Boyle's Law with the combined volume

## Notes

- All dimensions are in millimeters (mm)
- Pressures are in psi (pounds per square inch)
- Forces are in Newtons (N)
- Spring rates are in N/mm (Newtons per millimeter)
- The calculation assumes perfect annular geometry
- Travel range is from 0 to the configured maximum travel
- Calculations are performed at 10mm increments for display, 1mm increments for graphs

## Error Handling

The script includes error handling for:
- Missing or invalid YAML input files
- Missing required parameters
- Invalid travel values (automatically clamped to valid range)
