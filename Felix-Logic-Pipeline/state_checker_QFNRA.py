from z3 import *

# This updated code models continuous protein concentrations (Reals) 
# and steady-state invariants using QF_NRA arithmetic.

# here is the circuit that is being modeled
    #    IPTG (Inducer)                         aTc (Inducer)
    #         |                                      |
    #         v                                      v
    #   (deactivates)                          (deactivates)
    #         |                                      |
    #   +-----v-----+                          +-----v-----+
    #   |   LacI    | --------represses------| |   TetR    |
    #   | Repressor |                          | Repressor |
    #   +-----------+ |--represses--------     +-----------+

theoretical_circuits = [
    {
        "circuit_id": "Toggle_Switch_Stable",
        "description": "A well-designed genetic toggle switch. Tight promoters with low leakiness.",
        "parts": ["LacI", "TetR"],
        "interactions": [
            # LacI represses TetR tightly (basal expression is low)
            {"source": "LacI", "type": "represses", "target": "TetR", "max_expr": 100.0, "kd": 10.0, "basal": 1.0},
            
            # TetR represses LacI tightly
            {"source": "TetR", "type": "represses", "target": "LacI", "max_expr": 100.0, "kd": 10.0, "basal": 1.0}
        ],
        # Our goal: We want the switch to be capable of a state where LacI is dominant.
        "goal_state": {"LacI_min": 70.0, "TetR_max": 15.0} 
    },
    {
        "circuit_id": "Toggle_Switch_Leaky_Failure",
        "description": "A flawed toggle switch. The promoter for TetR is too leaky (high basal).",
        "parts": ["LacI", "TetR"],
        "interactions": [
            # promoter for TetR is leaky; even when fully repressed, it produces 30 units.
            {"source": "LacI", "type": "represses", "target": "TetR", "max_expr": 100.0, "kd": 10.0, "basal": 30.0},
            
            # TetR represses LacI normally
            {"source": "TetR", "type": "represses", "target": "LacI", "max_expr": 100.0, "kd": 10.0, "basal": 1.0}
        ],
        # ask the solver: is it physically possible to achieve the same dominant state?
        "goal_state": {"LacI_min": 70.0, "TetR_max": 15.0}
    },
    {
        "circuit_id": "Inducible_Toggle_Switch",
        "description": "A toggle switch on a medium-copy plasmid. LacI and TetR form dimers (n=2).",
        "parts": ["LacI", "TetR"],
        "plasmid_copy_number": 15, # Medium-copy plasmid (e.g., p15A origin)
        "interactions": [
            {
                "source": "LacI", "type": "represses", "target": "TetR", 
                "max_expr": 5.0, "kd": 10.0, "basal": 0.1, 
                "hill_coeff": 2, 
                "inducer": "IPTG", "k_inducer": 5.0 # IPTG binds to LacI, preventing it from repressing TetR
            },
            {
                "source": "TetR", "type": "represses", "target": "LacI", 
                "max_expr": 5.0, "kd": 10.0, "basal": 0.1, 
                "hill_coeff": 2, 
                "inducer": "aTc", "k_inducer": 5.0 # aTc binds to TetR, preventing it from repressing LacI
            }
        ],
        # Goal: We flush the cell with IPTG to force TetR to turn ON.
        "goal_state": {"IPTG": 50.0, "aTc": 0.0, "TetR_min": 50.0, "LacI_max": 10.0} 
    }
]

def run_continuous_verification(circuit):
    print(f"--- Verifying: {circuit['circuit_id']} ---")
    print(f"Description: {circuit['description']}")

    s = SolverFor('QF_NRA')
    copy_num = circuit.get('plasmid_copy_number', 1.0)
    
    # create Z3 Real vars for all parts to represent concentration levels
    vars = {name: Real(name) for name in circuit['parts']}
    
    # biological invariant: concentrations cannot be negative
    for name, var in vars.items():
        s.assert_and_track(var >= 0.0, f"BIOLOGY_INVARIANT: {name}_is_positive")
    
    # translate interactions --> constraints
    for i, inter in enumerate(circuit['interactions']):
        src = vars[inter['source']]
        tgt = vars[inter['target']]
        n = inter.get('hill_coeff', 1)

        max_expr = inter.get('max_expr', 50.0) * copy_num
        basal = inter.get('basal', 0.0) * copy_num
        kd = inter.get('kd', 5.0)
        
        inducer_name = inter.get('inducer')
        if inducer_name:
            # add inducer to the solver if it doesn't exist yet
            if inducer_name not in vars:
                vars[inducer_name] = Real(inducer_name)
                s.assert_and_track(vars[inducer_name] >= 0.0, f"BIOLOGY_INVARIANT: {inducer_name}_is_positive")
            
            inducer_var = vars[inducer_name]
            k_ind = inter.get('k_inducer', 1.0)
            
            # inducer binds the repressor. The "effective" repressor concentration drops.
            # effective_Src = Total_Src * (K_ind / (K_ind + Inducer_Conc))
            effective_src = src * (k_ind / (k_ind + inducer_var))
        else:
            effective_src = src

        label = f"BIOLOGY_RULE: {inter['source']} {inter['type']} {inter['target']}"
        
        if inter['type'] == "represses":
            # (Tgt - Basal) * (Kd^n + Effective_Src^n) == Max * Kd^n
            s.assert_and_track((tgt - basal) * ((kd**n) + (effective_src**n)) == (max_expr * (kd**n)), label)
            
        elif inter['type'] == "activates":
            # (Tgt - Basal) * (Kd^n + Effective_Src^n)  == (Max * Effective_Src^n)
            s.assert_and_track((tgt - basal) * ((kd**n) + (effective_src**n)) == (max_expr * effective_src**n), label)

    # apply goal state as concentration thresholds rather than strict equality
    if "goal_state" in circuit:
        for key, val in circuit['goal_state'].items():
                if "_" in key and ("min" in key or "max" in key):
                    part = key.split('_')[0]
                    if "min" in key:
                        s.assert_and_track(vars[part] >= val, f"GOAL: {part} >= {val}")
                    elif "max" in key:
                        s.assert_and_track(vars[part] <= val, f"GOAL: {part} <= {val}")
                else:
                    # extract concentrations (used for flushing inducers)
                    s.assert_and_track(vars[key] == val, f"ENVIRONMENT: {key} == {val}")

    result = s.check()
    
    if result == sat:
        model = s.model()
        print("Result: Consistent Design. State is reachable.")
        print("Steady-State Concentrations:")
        for var_name, var_obj in vars.items():
            val = model.eval(var_obj)
            if is_algebraic_value(val):
                val = val.approx(10) 
            elif is_rational_value(val):
                val = float(val.numerator_as_long()) / float(val.denominator_as_long())
            else:
                val = str(val)
            print(f"  • {var_name}: {float(val):.2f} units")
            
    elif result == unknown:
        print("Result: UNKNOWN.")
        print("Reason: The combination of tracking labels and non-linear math (exponents/division) is too complex for the general solver.")
        print("Tip: Switch to s = SolverFor('QF_NRA') and remove tracking to evaluate this circuit.")
        
    else:
        print("Result: INCONSISTENT. The circuit cannot sustain this state.")
        core = s.unsat_core()
        print("\n--- Conflict Analysis ---")
        for c in core:
            print(f"    • {c}")
        
    print("-" * 50 + "\n")

for circuit in theoretical_circuits:
    run_continuous_verification(circuit)