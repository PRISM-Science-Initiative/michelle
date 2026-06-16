import json
from z3 import *

# ==========================================
# 1. DATA LOADING & PREPARATION
# ==========================================

def load_parts_library(filepath="felix_database.json"):
    """Loads the JSON database and creates a fast-lookup dictionary."""
    with open(filepath, 'r') as file:
        data = json.load(file)
    
    parts_lookup = {part["name"]: part for part in data["parts"]}
    return data, parts_lookup


# ==========================================
# 2. CORE BIOLOGICAL & PHYSICAL RULES
# ==========================================

def apply_golden_gate_rules(solver, part_overhangs, construct_names=None):
    """
    Applies Golden Gate physical constraints with Unsat Core tracking labels.
    """
    junctions = []
    
    # 1. Linear Assembly Rule: Right overhang equals next Left overhang
    for i in range(len(part_overhangs) - 1):
        junction = part_overhangs[i][1]
        next_left = part_overhangs[i+1][0]
        
        # Create a dynamic label based on the parts being glued together
        if construct_names:
            label = f"Assembly_Mismatch:_Cannot_ligate_{construct_names[i]}_to_{construct_names[i+1]}"
        else:
            label = f"Linearity_Error_at_junction_{i}"
            
        solver.assert_and_track(junction == next_left, label)
        junctions.append(junction)

    # Add the flanking backbone overhangs
    junctions.append(part_overhangs[0][0])
    junctions.append(part_overhangs[-1][1])

    # 3. Length Rule: Overhangs must be exactly 4 base pairs
    for idx, j in enumerate(junctions):
        label = f"Invalid_Overhang_Length_at_part_index_{idx}"
        solver.assert_and_track(Length(j) == 4, label)

    # 4. Orthogonality Rule: All junctions must be strictly unique
    if len(junctions) > 1:
        solver.assert_and_track(Distinct(*junctions), "Orthogonality_Violation:_Duplicate_overhangs_found")

def calculate_circuit_parameters(construct, parts_lookup):
    """
    Calculates max_expr and basal leakiness using the Contextual Penalties model.
    """
    p_data, rbs_data, cds_data, t_data = None, None, None, None
    
    for part_name in construct:
        part = parts_lookup.get(part_name)
        if not part: continue
        
        if part["type"] == "promoter": p_data = part
        elif part["type"] == "rbs": rbs_data = part
        elif part["type"] == "cds": cds_data = part
        elif part["type"] == "terminator": t_data = part
        
    # Safety Check: If parts are missing, expression effectively fails to calculate
    if not all([p_data, rbs_data, cds_data, t_data]):
        return 0.0, 0.0, 0.0
        
    # Multiplicative Base
    base_expr = p_data["strength"] * rbs_data["strength"]
    
    # Contextual Penalties (Burden & Read-through)
    cds_burden_multiplier = 1000.0 / cds_data["length_bp"] if cds_data["length_bp"] > 1000 else 1.0
    terminator_multiplier = t_data["strength"] if t_data["strength"] < 0.9 else 1.0
        
    # Final Calculation
    predicted_expression = base_expr * cds_burden_multiplier * terminator_multiplier
    max_expr = predicted_expression * 100.0
    basal_leak = ((1.0 - p_data["strength"]) * 0.05) * 100.0
    
    return max_expr, basal_leak, predicted_expression


# ==========================================
# 3. PIPELINE FUNCTIONS
# ==========================================

def verify_experiment(experiment, parts_lookup):
    """Verifies if an experiment physically assembles and matches predicted math."""
    construct = experiment["construct"]
    print(f"--- Verifying {experiment.get('id', 'Unknown')} ---")
    print(f"Construct: {construct}")
    
    s = Solver()
    s.set(unsat_core=True)  # Turn on the tracking engine!
    
    # Extract hardcoded string values for the Golden Gate helper
    try:
        part_overhangs = [
            (StringVal(parts_lookup[p]["overhang_left"]), StringVal(parts_lookup[p]["overhang_right"])) 
            for p in construct
        ]
        # Pass the construct names so the labels are readable
        apply_golden_gate_rules(s, part_overhangs, construct)
    except KeyError as e:
        print(f"Result: INVALID \nConstruct contains undefined parts: {e}\n" + "-" * 50)
        return

    # Evaluate Biology
    max_expr, basal_leak, predicted_exp = calculate_circuit_parameters(construct, parts_lookup)
    dataset_expression = experiment.get("expression", 0.0)
    
    Z3_expr = Real('calculated_expression')
    s.add(Z3_expr == predicted_exp)
    
    # Track the expression bounds so we know if the math is what caused the failure
    s.assert_and_track(Z3_expr >= dataset_expression - 0.10, f"Expression_Too_Low_(Expected_~{dataset_expression})")
    s.assert_and_track(Z3_expr <= dataset_expression + 0.10, f"Expression_Too_High_(Expected_~{dataset_expression})")
    
    if s.check() == sat:
        print(f"Result: VALID")
        print(f"Physical assembly: SUCCESS")
        print(f"Expression match: YES (Predicted: {predicted_exp:.2f}, Actual: {dataset_expression})")
    else:
        print("Result: INVALID")
        core = s.unsat_core()
        print("--- Conflict Analysis ---")
        for conflict in core:
            # Clean up the Z3 label formatting (replaces underscores with spaces)
            reason = str(conflict).replace('_', ' ')
            print(f"  • {reason}")
    print("-" * 50)

def automate_circuit_design(target_expression, library, tolerance=0.10):
    """Uses Z3 to search the library for a physically valid circuit hitting a target expression."""
    print(f"--- Automated Design: Target Expression = {target_expression} ---")
    s = Solver()

    # Define Z3 identity variables
    p_name, rbs_name, cds_name, t_name = Strings('Promoter RBS CDS Terminator')
    p_str, rbs_str, t_str = Reals('p_strength rbs_strength t_strength')
    cds_len = Real('cds_length')

    # Define Z3 physical overhang variables
    p_L, p_R = Strings('p_L p_R')
    rbs_L, rbs_R = Strings('rbs_L rbs_R')
    cds_L, cds_R = Strings('cds_L cds_R')
    t_L, t_R = Strings('t_L t_R')

    # Apply centralized physical constraints
    apply_golden_gate_rules(s, [(p_L, p_R), (rbs_L, rbs_R), (cds_L, cds_R), (t_L, t_R)])

    # Map the JSON library to Z3 Constraints
    def encode_parts(name_var, str_var, left_var, right_var, part_list, length_var=None):
        part_constraints = []
        for part in part_list:
            conditions = [
                name_var == StringVal(part['name']),
                left_var == StringVal(part['overhang_left']),
                right_var == StringVal(part['overhang_right'])
            ]
            # CDS uses length for penalties; others use strength
            if length_var is not None:
                conditions.append(length_var == RealVal(part['length_bp']))
            if str_var is not None:
                conditions.append(str_var == RealVal(part['strength']))
                
            part_constraints.append(And(conditions))
        s.add(Or(part_constraints))

    encode_parts(p_name, p_str, p_L, p_R, library['promoter'])
    encode_parts(rbs_name, rbs_str, rbs_L, rbs_R, library['rbs'])
    encode_parts(cds_name, None, cds_L, cds_R, library['cds'], length_var=cds_len)
    encode_parts(t_name, t_str, t_L, t_R, library['terminator'])

    # Add Biological Math Constraints
    base_expr = p_str * rbs_str
    cds_penalty = If(cds_len > 1000.0, 1000.0 / cds_len, RealVal(1.0))
    t_penalty = If(t_str < 0.9, t_str, RealVal(1.0))
    
    final_expr = base_expr * cds_penalty * t_penalty

    # Set the Goal State
    s.add(final_expr >= RealVal(target_expression - tolerance))
    s.add(final_expr <= RealVal(target_expression + tolerance))

    # Solve and Extract
    if s.check() == sat:
        model = s.model()
        print("Result: DESIGN FOUND!")
        
        # Clean the strings before putting them in the f-string
        promoter_str = model[p_name].as_string().strip('"')
        rbs_str = model[rbs_name].as_string().strip('"')
        cds_str = model[cds_name].as_string().strip('"')
        t_str = model[t_name].as_string().strip('"')
        
        print(f"  • Promoter:   {promoter_str}")
        print(f"  • RBS:        {rbs_str}")
        print(f"  • CDS:        {cds_str}")
        print(f"  • Terminator: {t_str}")
        
        calculated_val = model.evaluate(final_expr)
        # Handle decimal extraction cleanly
        if str(calculated_val).isdigit():
            dec_val = float(calculated_val.as_long())
        else:
            dec_val = float(calculated_val.numerator_as_long()) / float(calculated_val.denominator_as_long())
            
        print(f"  • Predicted Expression: {dec_val:.3f}")
    else:
        print("Result: NO VALID DESIGN FOUND for this target with current parts.")
    print("-" * 50)


# ==========================================
# 4. EXECUTION DEMO
# ==========================================

if __name__ == "__main__":
    db_data, parts_dict = load_parts_library()
    
    # Pre-process the library for automation
    library = {
        "promoter": [p for p in db_data["parts"] if p["type"] == "promoter"],
        "rbs": [p for p in db_data["parts"] if p["type"] == "rbs"],
        "cds": [p for p in db_data["parts"] if p["type"] == "cds"],
        "terminator": [p for p in db_data["parts"] if p["type"] == "terminator"]
    }

    print("\n[ VERIFICATION ENGINE: DATASET EXPERIMENTS ]")
    print("="*50)
    
    actual_vals = []
    predicted_vals = []
    experiment_labels = []

    # Loop through ALL experiments in the dataset
    for exp in db_data["simulated_experiments"]:
        # Run the Z3 Verification
        verify_experiment(exp, parts_dict)
        
        # Calculate the raw numbers again just to save them for our graph
        _, _, predicted_exp = calculate_circuit_parameters(exp["construct"], parts_dict)
        actual_exp = exp.get("expression", 0.0)
        
        actual_vals.append(actual_exp)
        predicted_vals.append(predicted_exp)
        experiment_labels.append(exp["id"])

    print("\n[ VERIFICATION ENGINE: STRUCTURAL FAILURES ]")
    print("="*50)
    
    # 1. Missing Part (No RBS)
    # Fails because P_strong (TACT) cannot ligate directly to GFP (AATG)
    verify_experiment({
        "id": "exp_structural_missing_rbs",
        "construct": ["P_strong", "GFP", "T_strong"], 
        "expression": 0.82
    }, parts_dict)

    # 2. Out of Order Assembly
    # Fails because RBS (TACT->AATG) cannot precede Promoter (GGAG->TACT)
    verify_experiment({
        "id": "exp_structural_out_of_order",
        "construct": ["RBS_strong", "P_strong", "GFP", "T_strong"], 
        "expression": 0.82
    }, parts_dict)

    # 3. Position Mismatch (Mixing Position 1 and Position 2 parts)
    # Fails because P_strong leaves a 'TACT' overhang, but RBS2_strong expects a 'TACG' overhang!
    verify_experiment({
        "id": "exp_structural_position_mismatch",
        "construct": ["P_strong", "RBS2_strong", "GFP", "T_strong"], 
        "expression": 0.82
    }, parts_dict)


    print("\n[ AUTOMATED DESIGN ENGINE ]")
    print("="*50)
    
    # Target 1: Low Expression
    automate_circuit_design(0.20, library, tolerance=0.05)
    
    # Target 2: Medium Expression
    automate_circuit_design(0.50, library, tolerance=0.05)
    
    # Target 3: High Expression
    automate_circuit_design(0.85, library, tolerance=0.05)
    
    # Target 4: Impossible Expression (Proves the solver knows the system limits)
    automate_circuit_design(1.50, library, tolerance=0.05)


import matplotlib.pyplot as plt
import numpy as np

def plot_model_accuracy(actual, predicted, labels):
    """Generates a scatter plot comparing dataset expression vs. model predictions."""
    plt.figure(figsize=(12, 8))
    
    # Create the scatter plot
    plt.scatter(actual, predicted, color='dodgerblue', edgecolor='black', s=100, alpha=0.8, zorder=3)
    
    # Draw the "Perfect Accuracy" y = x line
    min_val = 0.0
    max_val = 1.0
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label="Perfect Match (y=x)", zorder=2)
    
    # Calculate R-squared to mathematically prove accuracy
    correlation_matrix = np.corrcoef(actual, predicted)
    correlation_xy = correlation_matrix[0,1]
    r_squared = correlation_xy**2
    
    # Formatting the Graph
    plt.title("Model Accuracy: Predicted vs. Actual Gene Expression", fontsize=20, pad=15)
    plt.xlabel("Actual Expression (Dataset Ground Truth)", fontsize=20)
    plt.ylabel("Predicted Expression (Z3 Contextual Model)", fontsize=20)
    plt.xlim(0, 1.0)
    plt.ylim(0, 1.0)
    plt.grid(True, linestyle=':', alpha=0.7, zorder=1)
    
    # Add an R^2 text box
    plt.text(0.05, 0.90, f"$R^2 = {r_squared:.3f}$", fontsize=15, 
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    
    plt.legend(loc="lower right")
    
    # Save the figure so you can put it on your poster
    plt.savefig("model_accuracy_plot.png", dpi=300, bbox_inches='tight')
    print("\nGraph saved as 'model_accuracy_plot.png'")
    plt.show()

# --- Call this at the very end of your script! ---
plot_model_accuracy(actual_vals, predicted_vals, experiment_labels)