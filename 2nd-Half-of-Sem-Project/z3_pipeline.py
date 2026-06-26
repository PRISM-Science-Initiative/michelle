import pandas as pd
import yaml
from z3 import *

def dynamic_pipeline(csv_path, yaml_path):
    # 1. Load the Configuration File
    print(f"Loading configuration framework from {yaml_path}...")
    with open(yaml_path, 'r') as file:
        config = yaml.safe_load(file)
    
    # 2. Dynamically Generate Z3 Variables
    print("Initializing dynamic Z3 symbolic variables...")
    vars_dict = {}
    for f in config['features']:
        name = f['name']
        t = f['type'].lower()
        if t == 'int':       vars_dict[name] = Int(name)
        elif t == 'real':     vars_dict[name] = Real(name)
        elif t == 'bool':     vars_dict[name] = Bool(name)
        elif t == 'string':   vars_dict[name] = String(name)
        else:
            raise ValueError(f"Unsupported data type '{t}' defined for feature '{name}'")

    # 3. Load Dataset and Dynamically Encode Rows
    print(f"Loading dataset from {csv_path}...")
    df = pd.read_csv(csv_path)
    
    print(f"Encoding {len(df)} rows into Z3 space using dynamic feature extraction...")
    row_constraints = []

    for _, row in df.iterrows():
        feature_assertions = []
        
        eval_env = {'row': row, 'pd': pd, 'len': len, 'any': any, 'str': str}
        
        for f in config['features']:
            name = f['name']
            t = f['type'].lower()
            
            # Dynamically compute the feature value
            computed_value = eval(f['source'], eval_env)
            
            # Construct equivalent Z3 constant literal
            if t == 'int':       z3_val = IntVal(int(computed_value))
            elif t == 'real':     z3_val = RealVal(str(computed_value))
            elif t == 'bool':     z3_val = BoolVal(bool(computed_value))
            elif t == 'string':   z3_val = StringVal(str(computed_value))
            
            feature_assertions.append(vars_dict[name] == z3_val)
            
        row_constraints.append(And(feature_assertions))
        
    dataset_constraint = Or(row_constraints)

    # 4. Parse Dynamic Rules from YAML
    print("Compiling SMT logical rules...")
    parsed_rules = []
    z3_env = {**globals(), **vars_dict} 
    
    for r in config.get('rules', []):
        name = r.get('name', 'Unnamed Rule')
        rule_type = r.get('type', 'simple')
        
        try:
            if rule_type == 'simple':
                z3_expr = eval(r['expr'], z3_env)
            elif rule_type == 'conditional':
                if_expr = eval(r['if'], z3_env)
                then_expr = eval(r['then'], z3_env)
                z3_expr = Implies(if_expr, then_expr)
            else:
                print(f"⚠️ Skipping unknown rule type '{rule_type}' for: {name}")
                continue
            parsed_rules.append({'name': name, 'expr': z3_expr})
        except Exception as e:
            print(f"  Error compiling logical rule '{name}': {e}")

    # 5. Execute Solvers
    solver = Solver()
    for rule in parsed_rules:
        solver.push()
        solver.add(dataset_constraint)
        solver.add(Not(rule['expr'])) # Look for a dataset violation counterexample
        
        result = solver.check()
        print(f"\n Verification: {rule['name']}")
        if result == unsat:
            print("     CONSISTENT: Dataset complies perfectly with this rule.")
        elif result == sat:
            print("     INCONSISTENT: Counterexample found violating this rule!")
            m = solver.model()
            for var in m.decls():
                print(f"      ↳ {var.name()}: {m[var]}")
        else:
            print("     UNKNOWN: Solver could not resolve consistency.")
        solver.pop()

if __name__ == "__main__":
    csv_file = "felix_sim_testdata.csv"
    yaml_file = "rules.yaml"
    dynamic_pipeline(csv_file, yaml_file)