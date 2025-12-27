"""
API service for Air Spring Calculator.

Provides REST API endpoints for calculating spring rates, forces, and pressures
for air spring chambers using Boyle's Law.
"""

from flask import Flask, request, jsonify, send_from_directory
from flask_cors import CORS
import sys
import os

# Add the current directory to the path to import calculate_volumes
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from calculate_volumes import (
    calculate_spring_rate,
    calculate_chamber_volumes,
    load_inputs,
    calculate_positive_chamber_pressure,
    calculate_negative_chamber_pressure,
    calculate_irt_chamber_pressure
)

app = Flask(__name__, static_folder='.')
CORS(app)  # Enable CORS for web applications


@app.route('/', methods=['GET'])
def index():
    """Serve the main webpage."""
    return send_from_directory('.', 'index.html')


@app.route('/health', methods=['GET'])
def health():
    """Health check endpoint."""
    return jsonify({
        'status': 'healthy',
        'service': 'Air Spring Calculator API'
    })


@app.route('/api/calculate', methods=['POST'])
def calculate():
    """
    Calculate spring rate and forces at a specific travel point.
    
    Request body (JSON):
    {
        "travel": 160,           # mm
        "startingPressure": 56,   # psi
        "irtPressure": 86,        # psi
        "travelPoint": 80         # mm (optional, defaults to 0)
    }
    
    Returns:
    {
        "travel": 80,
        "springRateNPerMm": 4.83,
        "springRateNPerM": 4830.0,
        "forceNetN": 1234.56,
        "forceMainN": 2000.0,
        "forceNegativeN": 765.44,
        "pressureMainPsi": 65.5,
        "pressureNegativePsi": 45.2,
        ...
    }
    """
    try:
        data = request.get_json()
        
        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400
        
        # Validate required fields
        required_fields = ['travel', 'startingPressure', 'irtPressure']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing required field: {field}'}), 400
        
        travel = float(data['travel'])
        starting_pressure = float(data['startingPressure'])
        irt_pressure = float(data['irtPressure'])
        travel_point = float(data.get('travelPoint', 0))
        
        # Validate values
        if travel <= 0 or starting_pressure < 0 or irt_pressure < 0:
            return jsonify({'error': 'Invalid parameter values'}), 400
        
        if travel_point < 0 or travel_point > travel:
            return jsonify({'error': 'travelPoint must be between 0 and travel'}), 400
        
        # Create inputs dictionary
        inputs = {
            'travel': travel,
            'startingPressure': starting_pressure,
            'irtPressure': irt_pressure
        }
        
        # Calculate chamber volumes
        results = calculate_chamber_volumes(inputs)
        
        # Calculate spring rate at the specified travel point
        spring_data = calculate_spring_rate(travel_point, results, inputs)
        
        # Prepare response
        response = {
            'travel': travel_point,
            'maxTravel': travel,
            'percentThroughTravel': (travel_point / travel) * 100 if travel > 0 else 0,
            'springRateNPerMm': round(spring_data['spring_rate_n_per_mm'], 2),
            'springRateNPerM': round(spring_data['spring_rate_n_per_m'], 2),
            'forceNetN': round(spring_data['force_net_n'], 2),
            'forceMainN': round(spring_data['force_positive_n'], 2),
            'forceNegativeN': round(spring_data['force_negative_n'], 2),
            'pressureMainPsi': round(spring_data['pressure_positive_psi'], 2),
            'pressureNegativePsi': round(spring_data['pressure_negative_psi'], 2),
            'volumeMainMm3': round(spring_data['volume_positive_mm3'], 2),
            'volumeNegativeMm3': round(spring_data['volume_negative_mm3'], 2)
        }
        
        # Add IRT data if present
        if 'force_irt_n' in spring_data:
            response['forceIrtN'] = round(spring_data['force_irt_n'], 2)
            response['pressureIrtPsi'] = round(spring_data['pressure_irt_psi'], 2)
            response['volumeIrtMm3'] = round(spring_data['volume_irt_mm3'], 2)
        
        return jsonify(response)
        
    except ValueError as e:
        return jsonify({'error': f'Invalid number format: {str(e)}'}), 400
    except Exception as e:
        return jsonify({'error': f'Calculation error: {str(e)}'}), 500


@app.route('/api/calculate-range', methods=['POST'])
def calculate_range():
    """
    Calculate spring rates and forces at multiple travel points.
    
    Request body (JSON):
    {
        "travel": 160,
        "startingPressure": 56,
        "irtPressure": 86,
        "increment": 10  # mm (optional, defaults to 10)
    }
    
    Returns:
    {
        "results": [
            {
                "travel": 0,
                "springRateNPerMm": 5.45,
                "forceNetN": 1000.0,
                ...
            },
            ...
        ]
    }
    """
    try:
        data = request.get_json()
        
        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400
        
        # Validate required fields
        required_fields = ['travel', 'startingPressure', 'irtPressure']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing required field: {field}'}), 400
        
        travel = float(data['travel'])
        starting_pressure = float(data['startingPressure'])
        irt_pressure = float(data['irtPressure'])
        increment = float(data.get('increment', 10))
        
        # Validate values
        if travel <= 0 or starting_pressure < 0 or irt_pressure < 0 or increment <= 0:
            return jsonify({'error': 'Invalid parameter values'}), 400
        
        # Create inputs dictionary
        inputs = {
            'travel': travel,
            'startingPressure': starting_pressure,
            'irtPressure': irt_pressure
        }
        
        # Calculate chamber volumes
        results = calculate_chamber_volumes(inputs)
        
        # Generate travel points
        travel_points = []
        t = 0
        while t <= travel:
            travel_points.append(t)
            t += increment
        if travel_points[-1] != travel:
            travel_points.append(travel)
        
        # Calculate results for each travel point
        results_list = []
        for travel_point in travel_points:
            spring_data = calculate_spring_rate(travel_point, results, inputs)
            
            result_item = {
                'travel': round(travel_point, 1),
                'percentThroughTravel': round((travel_point / travel) * 100 if travel > 0 else 0, 1),
                'springRateNPerMm': round(spring_data['spring_rate_n_per_mm'], 2),
                'springRateNPerM': round(spring_data['spring_rate_n_per_m'], 2),
                'forceNetN': round(spring_data['force_net_n'], 2),
                'forceMainN': round(spring_data['force_positive_n'], 2),
                'forceNegativeN': round(spring_data['force_negative_n'], 2),
                'pressureMainPsi': round(spring_data['pressure_positive_psi'], 2),
                'pressureNegativePsi': round(spring_data['pressure_negative_psi'], 2)
            }
            
            if 'force_irt_n' in spring_data:
                result_item['forceIrtN'] = round(spring_data['force_irt_n'], 2)
                result_item['pressureIrtPsi'] = round(spring_data['pressure_irt_psi'], 2)
            
            results_list.append(result_item)
        
        return jsonify({
            'maxTravel': travel,
            'increment': increment,
            'results': results_list
        })
        
    except ValueError as e:
        return jsonify({'error': f'Invalid number format: {str(e)}'}), 400
    except Exception as e:
        return jsonify({'error': f'Calculation error: {str(e)}'}), 500


@app.route('/api/chamber-volumes', methods=['POST'])
def get_chamber_volumes():
    """
    Get initial chamber volumes and configuration.
    
    Request body (JSON):
    {
        "travel": 160,
        "startingPressure": 56,
        "irtPressure": 86
    }
    
    Returns chamber volumes and configuration.
    """
    try:
        data = request.get_json()
        
        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400
        
        required_fields = ['travel', 'startingPressure', 'irtPressure']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing required field: {field}'}), 400
        
        inputs = {
            'travel': float(data['travel']),
            'startingPressure': float(data['startingPressure']),
            'irtPressure': float(data['irtPressure'])
        }
        
        results = calculate_chamber_volumes(inputs)
        
        response = {
            'travel': results['travel'],
            'startingPressure': results['starting_pressure'],
            'annularAreaMm2': round(results['annular_area'], 2),
            'mainChamber': {
                'length': results['positive_chamber']['length'],
                'volumeMm3': round(results['positive_chamber']['volume_mm3'], 2),
                'volumeLiters': round(results['positive_chamber']['volume_liters'], 6),
                'volumeM3': round(results['positive_chamber']['volume_m3'], 9)
            },
            'negativeChamber': {
                'length': results['negative_chamber']['length'],
                'volumeMm3': round(results['negative_chamber']['volume_mm3'], 2),
                'volumeLiters': round(results['negative_chamber']['volume_liters'], 6),
                'volumeM3': round(results['negative_chamber']['volume_m3'], 9)
            }
        }
        
        if 'irt_chamber' in results:
            response['irtChamber'] = {
                'length': results['irt_chamber']['length'],
                'startingPressure': results['irt_chamber']['starting_pressure'],
                'volumeMm3': round(results['irt_chamber']['volume_mm3'], 2),
                'volumeLiters': round(results['irt_chamber']['volume_liters'], 6),
                'volumeM3': round(results['irt_chamber']['volume_m3'], 9)
            }
        
        return jsonify(response)
        
    except ValueError as e:
        return jsonify({'error': f'Invalid number format: {str(e)}'}), 400
    except Exception as e:
        return jsonify({'error': f'Calculation error: {str(e)}'}), 500


@app.route('/api/pressures', methods=['POST'])
def get_pressures():
    """
    Get chamber pressures at a specific travel point.
    
    Request body (JSON):
    {
        "travel": 160,
        "startingPressure": 56,
        "irtPressure": 86,
        "travelPoint": 80
    }
    
    Returns detailed pressure and volume information.
    """
    try:
        data = request.get_json()
        
        if not data:
            return jsonify({'error': 'No JSON data provided'}), 400
        
        required_fields = ['travel', 'startingPressure', 'irtPressure', 'travelPoint']
        for field in required_fields:
            if field not in data:
                return jsonify({'error': f'Missing required field: {field}'}), 400
        
        inputs = {
            'travel': float(data['travel']),
            'startingPressure': float(data['startingPressure']),
            'irtPressure': float(data['irtPressure'])
        }
        
        travel_point = float(data['travelPoint'])
        results = calculate_chamber_volumes(inputs)
        
        # Get IRT info
        has_irt = 'irt_chamber' in results
        irt_vol = results['irt_chamber']['volume_mm3'] if has_irt else 0
        irt_pres = results['irt_chamber']['starting_pressure'] if has_irt else 0
        
        # Calculate pressures
        main_data = calculate_positive_chamber_pressure(
            results['positive_chamber']['starting_pressure'],
            results['positive_chamber']['initial_volume_mm3'],
            results['positive_chamber']['length'],
            travel_point,
            results['annular_area'],
            irt_volume=irt_vol,
            irt_pressure=irt_pres
        )
        
        neg_data = calculate_negative_chamber_pressure(
            results['negative_chamber']['starting_pressure'],
            results['negative_chamber']['volume_mm3'],
            results['negative_chamber']['length'],
            travel_point,
            results['annular_area']
        )
        
        response = {
            'travel': travel_point,
            'mainChamber': {
                'pressurePsi': round(main_data['pressure_psi'], 2),
                'volumeMm3': round(main_data['volume_mm3'], 2),
                'volumeLiters': round(main_data['volume_liters'], 6),
                'length': round(main_data['length'], 2),
                'irtCombined': main_data.get('irt_combined', False)
            },
            'negativeChamber': {
                'pressurePsi': round(neg_data['pressure_psi'], 2),
                'volumeMm3': round(neg_data['volume_mm3'], 2),
                'volumeLiters': round(neg_data['volume_liters'], 6),
                'length': round(neg_data['length'], 2)
            }
        }
        
        if has_irt:
            irt_combined = main_data.get('irt_combined', False)
            if not irt_combined:
                irt_data = calculate_irt_chamber_pressure(
                    results['irt_chamber']['starting_pressure'],
                    results['irt_chamber']['volume_mm3'],
                    results['irt_chamber']['length'],
                    travel_point,
                    results['annular_area']
                )
                response['irtChamber'] = {
                    'pressurePsi': round(irt_data['pressure_psi'], 2),
                    'volumeMm3': round(irt_data['volume_mm3'], 2),
                    'volumeLiters': round(irt_data['volume_liters'], 6),
                    'length': round(irt_data['length'], 2)
                }
            else:
                response['irtChamber'] = {
                    'pressurePsi': round(main_data['pressure_psi'], 2),
                    'volumeMm3': 0,
                    'volumeLiters': 0,
                    'length': results['irt_chamber']['length'],
                    'combined': True
                }
        
        return jsonify(response)
        
    except ValueError as e:
        return jsonify({'error': f'Invalid number format: {str(e)}'}), 400
    except Exception as e:
        return jsonify({'error': f'Calculation error: {str(e)}'}), 500


@app.errorhandler(404)
def not_found(error):
    return jsonify({'error': 'Endpoint not found'}), 404


@app.errorhandler(500)
def internal_error(error):
    return jsonify({'error': 'Internal server error'}), 500


if __name__ == '__main__':
    # Run the API service
    # Default: http://localhost:5000
    import os
    debug_mode = os.getenv('FLASK_ENV') != 'production'
    app.run(debug=debug_mode, host='0.0.0.0', port=5000)

