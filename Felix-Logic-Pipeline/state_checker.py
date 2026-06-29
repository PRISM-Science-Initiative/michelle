from z3 import *

# This code provides logical circuits/interactions within a cell. It also contains a 
# "test state" and determines if the "test state" is allowed/reachable within the
# constraints of the circuit/biological interaction rules.

# if it is not reachable, it should print what the conflict is. Otherwise, the state
# is allowed/reachable and is printed out as is.
 
theoretical_circuits = [
    {
        "circuit_id": "C01_Inverter",
        "description": "Simple NOT gate. Input A represses Output B.",
        "parts": ["Protein_A", "Promoter_B", "Output_B"],
        "interactions": [
            {"source": "Protein_A", "type": "represses", "target": "Promoter_B"},
            {"source": "Promoter_B", "type": "produces", "target": "Output_B"}
        ],
        "test_state": {"Protein_A": True, "Output_B": False} # should be good!
    },
    {
        "circuit_id": "C02_FeedForward",
        "description": "Coherent Feed-Forward Loop. A activates B and C; B activates C.",
        "parts": ["A", "B", "C"],
        "interactions": [
            {"source": "A", "type": "activates", "target": "B"},
            {"source": "B", "type": "activates", "target": "C"},
            {"source": "A", "type": "activates", "target": "C"}
        ],
        "test_state": {"A": True, "C": True} # should also be good
    },
    {
        "circuit_id": "C03_Paradox_Loop",
        "description": "A flawed circuit for Figure 3C. A represses B, B represses A, both required.",
        "parts": ["A", "B"],
        "interactions": [
            {"source": "A", "type": "represses", "target": "B"},
            {"source": "B", "type": "represses", "target": "A"}
        ],
        "test_state": {"A": True, "B": True} # unsat!
    },
    {
        "circuit_id": "C04_Negative_Feedback",
        "description": "Protein A starts the process, but the end product C shuts it down.",
        "parts": ["Protein_A", "Gene_B", "Protein_C"],
        "interactions": [
            {"source": "Protein_A", "type": "activates", "target": "Gene_B"},
            {"source": "Gene_B", "type": "produces", "target": "Protein_C"},
            {"source": "Protein_C", "type": "represses", "target": "Protein_A"}
        ],
        "test_state": {"Protein_A": True} # unsat!
    },
    {
        "circuit_id": "C05_OR_Gate",
        "description": "Either A OR B can turn on the Output.",
        "parts": ["A", "B", "Output"],
        "interactions": [
            {"source": ["A", "B"], "type": "or_activates", "target": "Output"}
        ],
        "test_state": {"A": False, "B": False, "Output": True} # unsat!
    },
    {
        "circuit_id": "C06_AND_Gate",
        "description": "Both A AND B are required to turn on the Output.",
        "parts": ["A", "B", "Output"],
        "interactions": [
            {"source": ["A", "B"], "type": "and_activates", "target": "Output"}
        ],
        "test_state": {"A": True, "B": False, "Output": True} # unsat!
    },
    {
        "circuit_id": "C07_Complex_Controller",
        "description": "Feed-forward loop with a master inhibitor. High-level decision logic.",
        "parts": ["Master_Switch", "Protein_A", "Protein_B", "Inhibitor", "Final_Output"],
        "interactions": [
            # the master switch turns on the first worker
            {"source": "Master_Switch", "type": "activates", "target": "Protein_A"},
            
            # protein A turns on protein B
            {"source": "Protein_A", "type": "activates", "target": "Protein_B"},
            
            # output needs both A and B (coherent feed-forward)
            {"source": ["Protein_A", "Protein_B"], "type": "and_activates", "target": "Final_Output"},
            
            # inhibitor shuts down the whole chain
            {"source": "Inhibitor", "type": "represses", "target": "Protein_A"}
        ],
        "test_state": {"Master_Switch": True, "Inhibitor": True, "Final_Output": True}
    }
]

def run_verification(circuit):
    print(f"--- Verifying: {circuit['circuit_id']} ---")
    print(f"Description: {circuit['description']}")
    s = Solver()
    
    # create Z3 Variables for all parts
    # use a dictionary to store Z3 objects: {'A': Bool('A'), ...}
    vars = {name: Bool(name) for name in circuit['parts']}
    
    # translate interactions
    for i, inter in enumerate(circuit['interactions']):
        src_names = inter['source'] if isinstance(inter['source'], list) else [inter['source']]
        src_vars = [vars[name] for name in src_names]
        tgt_var = vars[inter['target']]

        sources_str = ", ".join(src_names)
        label = f"BIOLOGY_RULE: {sources_str} {inter['type']} {inter['target']}"
        
        if inter['type'] == "represses":
            s.assert_and_track(Implies(Or(*src_vars), Not(tgt_var)), label)
            # if any source is true, target is false
        
        elif inter['type'] == "activates":
            s.assert_and_track(Implies(Or(*src_vars), tgt_var), label)
            # if any source is true, target is true
            
        elif inter['type'] == "produces":
            # promoter -> protein
            s.assert_and_track(src_vars[0] == tgt_var, label)
        
        elif inter['type'] == "and_activates":
            # target is true IF AND ONLY IF all sources are true
            s.assert_and_track(tgt_var == And(*src_vars), label)

        elif inter['type'] == "or_activates":
            # target is true IF AND ONLY IF at least one source is true
            s.assert_and_track(tgt_var == Or(*src_vars), label)

    # apply the specific scenario/data to test
    if "test_state" in circuit:
        for part, state in circuit['test_state'].items():
            # track test states too so we can see which one fails
            s.assert_and_track(vars[part] == state, f"Constraint_{part}_{state}")

    # check consistency of the given state in our cell
    if s.check() == sat:
        print("Result: Consistent Design.")
        print("Stable Model:", s.model())
    else:
        print("Result: INCONSISTENT (Paradox Found).")
        core = s.unsat_core()
        print("\n--- Conflict Analysis ---")
        print("The solver found a contradiction between these specific elements:")
        
        rules = [str(c) for c in core if "BIOLOGY_RULE" in str(c)]
        constraints = [str(c) for c in core if "Constraint" in str(c)]
        
        print("  From the Circuit Design:")
        for r in rules: print(f"    • {r}")
            
        print("  From the Requested State:")
        for c in constraints: print(f"    • {c}")
        
        print("\nConclusion: The physical rules of the circuit do not allow this combination of states.")
    print("\n")

# run the simulation on the dataset
for circuit in theoretical_circuits:
    run_verification(circuit)