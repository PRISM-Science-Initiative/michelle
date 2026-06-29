from z3 import *

# This code provides logical circuits/interactions within a cell. It then takes a 
# "goal state" and determines if the "goal state" is reachable within the
# constraints of the circuit/biological interaction rules.

# if it is not reachable, it ____ oops idk yet still in progress ____ . Otherwise, the state
# that allows us to reach our goals is printed out
 
theoretical_circuits = [
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
        "test_state": {"Master_Switch": True, "Inhibitor": True, "Final_Output": True}, # unsat!
        "goal_state": {"Final_Output": True}
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
            # changed implies to == so proteins don't spontaneously turn on
            s.assert_and_track((Or(*src_vars) == tgt_var), label)
            
        elif inter['type'] == "produces":
            # promoter -> protein
            s.assert_and_track(src_vars[0] == tgt_var, label)
        
        elif inter['type'] == "and_activates":
            # target is true IF AND ONLY IF all sources are true
            s.assert_and_track(tgt_var == And(*src_vars), label)

        elif inter['type'] == "or_activates":
            # target is true IF AND ONLY IF at least one source is true
            s.assert_and_track(tgt_var == Or(*src_vars), label)

    # set the GOAL (what we want the cell to do)
    if "goal_state" in circuit:
        for part, state in circuit['goal_state'].items():
            s.add(vars[part] == state)
            print(f"Targeting Goal: {part} must be {state}")

    # solve for the required inputs
    if s.check() == sat:
        model = s.model()
        print("Result: Solution Found! To reach the goal, the cell needs:")
        
        # We only care about the "Inputs" (Master Switch and Inhibitor)
        inputs = ["Master_Switch", "Inhibitor"] # CHANGE THESE WHEN YOU ADD MORE EXAMPLE CIRCUITS!!!!!!
        for var in inputs:
            print(f"  • {var}: {model[vars[var]]}")
    else:
        print("Result: Impossible. No biological state satisfies this goal.")
    print("-" * 40 + "\n")
    print("\n")

# run the simulation on the dataset
for circuit in theoretical_circuits:
    run_verification(circuit)