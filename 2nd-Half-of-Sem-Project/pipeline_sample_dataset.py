import csv
from z3 import *
import matplotlib.pyplot as plt
import numpy as np

# mock up running on a sample dataset. Split up further into rules and z3_pipeline

# ==========================================
# 1. DATA LOADING & PREPARATION (CSV ENGINE)
# ==========================================

def load_parts_library(filepath="felix_sim_testdata.csv"):
    """
    Loads the CSV database, extracts unique genetic parts from sequences, 
    and prepares the simulated experiment dataset.
    """
    print(f"Loading data from CSV: {filepath}")
    
    parts_list = []
    simulated_experiments = []
    
    unique_promoters = set()
    unique_utr5s = set()
    unique_cdss = set()
    unique_terms = set()
    
    with open(filepath, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file)
        for i, row in enumerate(reader):
            p_seq = row['promoter'].strip()
            r_seq = row['utr5'].strip()  # 5' UTR acts as the RBS
            c_seq = row['cds'].strip()
            t_seq = row['term'].strip()
            
            # 1. Process Unique Promoters
            if p_seq not in unique_promoters:
                unique_promoters.add(p_seq)
                parts_list.append({
                    "name": p_seq,
                    "type": "promoter",
                    "strength": round(0.3 + (len(unique_promoters) % 7) * 0.1, 2),  # Create a gradient of strengths
                    "overhang_left": "GGAG",
                    "overhang_right": "TACT",
                    "length_bp": int(row['promoter_len'])
                })
                
            # 2. Process Unique UTR5s (RBS)
            if r_seq not in unique_utr5s:
                unique_utr5s.add(r_seq)
                parts_list.append({
                    "name": r_seq,
                    "type": "rbs",
                    "strength": round(0.4 + (len(unique_utr5s) % 6) * 0.1, 2),  # Create a gradient of strengths
                    "overhang_left": "TACT",
                    "overhang_right": "AATG",
                    "length_bp": int(row['utr5_len'])
                })
                
            # 3. Process Unique CDSs
            if c_seq not in unique_cdss:
                unique_cdss.add(c_seq)
                parts_list.append({
                    "name": c_seq,
                    "type": "cds",
                    "strength": 1.0,  # strengths are all 1
                    "overhang_left": "AATG",
                    "overhang_right": "GCTT",
                    "length_bp": int(row['cds_len'])
                })
                
            # 4. Process Unique Terminators
            if t_seq not in unique_terms:
                unique_terms.add(t_seq)
                parts_list.append({
                    "name": t_seq,
                    "type": "terminator",
                    "strength": round(0.1 + (len(unique_terms) % 5) * 0.1, 2),  # Create a gradient of strengths
                    "overhang_left": "GCTT",
                    "overhang_right": "CGGA",
                    "length_bp": int(row['term_len'])
                })
            
            # 5. Extract the Experiment
            simulated_experiments.append({
                "id": row['construct_id'].strip(),
                "construct": [p_seq, r_seq, c_seq, t_seq],
                "expression": float(row['expr_mean'])
            })
            
    parts_lookup = {part["name"]: part for part in parts_list}
    print(f"Successfully cataloged {len(parts_list)} unique parts.")
    print(f"Loaded {len(simulated_experiments)} experimental constructs.\n")
    return {"parts": parts_list, "simulated_experiments": simulated_experiments}, parts_lookup


# ==========================================
# 2. CORE BIOLOGICAL & PHYSICAL RULES
# ==========================================

def apply_golden_gate_rules(solver, part_overhangs, construct_names=None):
    """Applies Golden Gate physical constraints with Unsat Core tracking labels."""
    junctions = []
    
    # Linear Assembly Rule: Right overhang equals next Left overhang
    for i in range(len(part_overhangs) - 1):
        junction = part_overhangs[i][1]
        next_left = part_overhangs[i+1][0]
        
        if construct_names:
            label = f"Assembly_Mismatch:_Cannot_ligate_Part_{i}_to_Part_{i+1}"
        else:
            label = f"Linearity_Error_at_junction_{i}"
            
        solver.assert_and_track(junction == next_left, label)
        junctions.append(junction)

    # Add the flanking backbone overhangs
    junctions.append(part_overhangs[0][0])
    junctions.append(part_overhangs[-1][1])

    # Length Rule: Overhangs must be exactly 4 base pairs
    for idx, j in enumerate(junctions):
        label = f"Invalid_Overhang_Length_at_part_index_{idx}"
        solver.assert_and_track(Length(j) == 4, label)

    # Orthogonality Rule: All junctions must be strictly unique
    if len(junctions) > 1:
        solver.assert_and_track(Distinct(*junctions), "Orthogonality_Violation:_Duplicate_overhangs_found")

def calculate_circuit_parameters(construct, parts_lookup):
    """Calculates max_expr and basal leakiness using the Contextual Penalties model."""
    p_data, rbs_data, cds_data, t_data = None, None, None, None
    
    for part_name in construct:
        part = parts_lookup.get(part_name)
        if not part: continue
        
        if part["type"] == "promoter": p_data = part
        elif part["type"] == "rbs": rbs_data = part
        elif part["type"] == "cds": cds_data = part
        elif part["type"] == "terminator": t_data = part
        
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
    
    s = Solver()
    s.set(unsat_core=True)  # Turn on tracking!
    
    try:
        part_overhangs = [
            (StringVal(parts_lookup[p]["overhang_left"]), StringVal(parts_lookup[p]["overhang_right"])) 
            for p in construct
        ]
        apply_golden_gate_rules(s, part_overhangs, construct)
    except KeyError as e:
        print(f"Result: INVALID \nConstruct contains undefined sequence parts!\n" + "-" * 50)
        return

    max_expr, basal_leak, predicted_exp = calculate_circuit_parameters(construct, parts_lookup)
    dataset_expression = experiment.get("expression", 0.0)
    
    Z3_expr = Real('calculated_expression')
    s.add(Z3_expr == predicted_exp)
    
    s.assert_and_track(Z3_expr >= dataset_expression - 0.25, f"Expression_Too_Low_(Expected_~{dataset_expression:.3f})")
    s.assert_and_track(Z3_expr <= dataset_expression + 0.25, f"Expression_Too_High_(Expected_~{dataset_expression:.3f})")
    
    if s.check() == sat:
        print(f"Result: VALID")
        print(f"Physical assembly: SUCCESS")
        print(f"Expression match: YES (Predicted: {predicted_exp:.3f}, Actual: {dataset_expression:.3f})")
    else:
        print("Result: INVALID")
        core = s.unsat_core()
        print("--- Conflict Analysis ---")
        for conflict in core:
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

    apply_golden_gate_rules(s, [(p_L, p_R), (rbs_L, rbs_R), (cds_L, cds_R), (t_L, t_R)])

    def encode_parts(name_var, str_var, left_var, right_var, part_list, length_var=None):
        part_constraints = []
        for part in part_list:
            conditions = [
                name_var == StringVal(part['name']),
                left_var == StringVal(part['overhang_left']),
                right_var == StringVal(part['overhang_right'])
            ]
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

    base_expr = p_str * rbs_str
    cds_penalty = If(cds_len > 1000.0, 1000.0 / cds_len, RealVal(1.0))
    t_penalty = If(t_str < 0.9, t_str, RealVal(1.0))
    final_expr = base_expr * cds_penalty * t_penalty

    s.add(final_expr >= RealVal(target_expression - tolerance))
    s.add(final_expr <= RealVal(target_expression + tolerance))

    if s.check() == sat:
        model = s.model()
        print("Result: DESIGN FOUND!")
        
        p_seq = model[p_name].as_string().strip('"')
        r_seq = model[rbs_name].as_string().strip('"')
        c_seq = model[cds_name].as_string().strip('"')
        t_seq = model[t_name].as_string().strip('"')
        
        # Displaying truncated segments for clean terminal output
        print(f"  • Promoter:   {p_seq[:15]}... ({len(p_seq)} bp)")
        print(f"  • RBS/UTR5:   {r_seq[:15]}... ({len(r_seq)} bp)")
        print(f"  • CDS:        {c_seq[:15]}... ({len(c_seq)} bp)")
        print(f"  • Terminator: {t_seq[:15]}... ({len(t_seq)} bp)")
        
        calculated_val = model.evaluate(final_expr)
        if str(calculated_val).isdigit():
            dec_val = float(calculated_val.as_long())
        else:
            dec_val = float(calculated_val.numerator_as_long()) / float(calculated_val.denominator_as_long())
            
        print(f"  • Predicted Expression Value: {dec_val:.3f}")
    else:
        print("Result: NO VALID DESIGN FOUND for this target with current parts library rules.")
    print("-" * 50)


def plot_model_accuracy(actual, predicted, labels):
    """Generates a scatter plot comparing dataset expression vs. model predictions."""
    plt.figure(figsize=(10, 6))
    plt.scatter(actual, predicted, color='dodgerblue', edgecolor='black', s=40, alpha=0.6, zorder=3)
    
    min_val, max_val = min(actual + predicted), max(actual + predicted)
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', linewidth=2, label="Perfect Match (y=x)", zorder=2)
    
    correlation_matrix = np.corrcoef(actual, predicted)
    r_squared = correlation_matrix[0,1]**2 if not np.isnan(correlation_matrix[0,1]) else 0.0
    
    plt.title("Model Accuracy Check: CSV Mean Expression vs Z3 Model Predictions", fontsize=14, pad=15)
    plt.xlabel("Actual Expression (CSV 'expr_mean')", fontsize=12)
    plt.ylabel("Predicted Expression (Z3 Multiplicative Model)", fontsize=12)
    plt.grid(True, linestyle=':', alpha=0.7, zorder=1)
    
    plt.text(0.05, 0.85, f"$R^2 = {r_squared:.3f}$", fontsize=12, transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.5'))
    plt.legend(loc="lower right")
    
    plt.savefig("model_accuracy_plot.png", dpi=300, bbox_inches='tight')
    print("\nGraph saved successfully as 'model_accuracy_plot.png'")
    plt.show()


# ==========================================
# 4. EXECUTION DEMO
# ==========================================

if __name__ == "__main__":
    # Make sure this matches your local path to the file!
    db_data, parts_dict = load_parts_library("felix_sim_testdata.csv")
    
    library = {
        "promoter": [p for p in db_data["parts"] if p["type"] == "promoter"],
        "rbs": [p for p in db_data["parts"] if p["type"] == "rbs"],
        "cds": [p for p in db_data["parts"] if p["type"] == "cds"],
        "terminator": [p for p in db_data["parts"] if p["type"] == "terminator"]
    }

    print("[ VERIFICATION ENGINE: DATASET EXPERIMENTS (Sampled Subset) ]")
    print("="*60)
    
    actual_vals = []
    predicted_vals = []
    experiment_labels = []

    # Process first 30 experimental entries for validation profiling
    for exp in db_data["simulated_experiments"][:30]:
        verify_experiment(exp, parts_dict)
        
    # Process full list for graph coordinates
    for exp in db_data["simulated_experiments"]:
        _, _, predicted_exp = calculate_circuit_parameters(exp["construct"], parts_dict)
        actual_vals.append(exp["expression"])
        predicted_vals.append(predicted_exp)
        experiment_labels.append(exp["id"])

    print("\n[ VERIFICATION ENGINE: STRUCTURAL FAILURES DEMO ]")
    print("="*60)
    
    # Extract dynamic samples from your loaded elements
    sample_p = library["promoter"][0]["name"]
    sample_r = library["rbs"][0]["name"]
    sample_c = library["cds"][0]["name"]
    sample_t = library["terminator"][0]["name"]

    # 1. Missing Part (No RBS)
    verify_experiment({
        "id": "exp_structural_missing_rbs",
        "construct": [sample_p, sample_c, sample_t], 
        "expression": 0.5
    }, parts_dict)

    # 2. Out of Order Assembly
    verify_experiment({
        "id": "exp_structural_out_of_order",
        "construct": [sample_r, sample_p, sample_c, sample_t], 
        "expression": 0.5
    }, parts_dict)

    # 3. Position Mismatch Simulation
    mismatch_rbs = "MOCK_SEQUENCE_BREAKING_OVERHANG_RULES_TACT_EXPECTED"
    parts_dict[mismatch_rbs] = {
        "name": mismatch_rbs, "type": "rbs", "strength": 0.5,
        "overhang_left": "AAAA",  # Explicitly breaks ligation junction with promoter (TACT)
        "overhang_right": "AATG", "length_bp": 30
    }
    verify_experiment({
        "id": "exp_structural_position_mismatch",
        "construct": [sample_p, mismatch_rbs, sample_c, sample_t], 
        "expression": 0.5
    }, parts_dict)


    print("\n[ AUTOMATED DESIGN ENGINE ]")
    print("="*60)
    
    # Request combinations that hit different mathematical thresholds from the loaded library
    automate_circuit_design(0.15, library, tolerance=0.02)
    automate_circuit_design(0.35, library, tolerance=0.02)
    automate_circuit_design(0.55, library, tolerance=0.02)
    automate_circuit_design(1.95, library, tolerance=0.02)  # Intentional impossible calculation edge-case

    # Generate accuracy distribution plot
    plot_model_accuracy(actual_vals, predicted_vals, experiment_labels)