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
    
    # MODIFIED: Track row metadata (IDs) paired with their specific Z3 assertions
    row_records = []

    for _, row in df.iterrows():
        feature_assertions = []
        eval_env = {'row': row, 'pd': pd, 'len': len, 'any': any, 'str': str}
        
        for f in config['features']:
            name = f['name']
            t = f['type'].lower()
            
            computed_value = eval(f['source'], eval_env)
            
            if t == 'int':       z3_val = IntVal(int(computed_value))
            elif t == 'real':     z3_val = RealVal(str(computed_value))
            elif t == 'bool':     z3_val = BoolVal(bool(computed_value))
            elif t == 'string':   z3_val = StringVal(str(computed_value))
            
            feature_assertions.append(vars_dict[name] == z3_val)
            
        # Store metadata along with the concrete row constraints
        row_records.append({
            'gene_id': row.get('gene_id', 'Unknown'),
            'gene_name': row.get('gene_name', 'Unknown'),
            'constraints': And(feature_assertions)
        })
        
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

    # 5. Execute Solvers (MODIFIED to find ALL violations per rule)
    print("\nEvaluating rules row-by-row to collect all violations...")
    
    for rule in parsed_rules:
        violating_sequences = []
        
        # Test every sequence constraint individually against the negated rule
        for record in row_records:
            solver = Solver()
            solver.add(record['constraints'])
            solver.add(Not(rule['expr'])) # Look for a violation counterexample
            
            if solver.check() == sat:
                # If SAT, this specific sequence broke the rule
                violating_sequences.append(f"{record['gene_id']} ({record['gene_name']})")
                
        print(f"\n Verification: {rule['name']}")
        if len(violating_sequences) == 0:
            print("     CONSISTENT: Dataset complies perfectly with this rule.")
        else:
            print(f"     INCONSISTENT: Found {len(violating_sequences)} violating sequences!")
            print(f"     ↳ Violators: {', '.join(violating_sequences)}")

    # ADDED TASK: Specific analysis of non-standard start codons using Pandas
    print("\n" + "="*50)
    print("ADDITIONAL ANALYSIS: Non-Standard Start Codon Breakdown")
    print("="*50)
    
    if 'cds_seq' in df.columns:
        # Extract the first 3 nucleotides of each CDS sequence
        df['extracted_start_codon'] = df['cds_seq'].astype(str).str.upper().str[:3]
        
        # Define standard bacterial start codons based on your rules configuration
        standard_start_codons = ['ATG', 'GTG', 'TTG']
        
        # Filter for rows that don't start with standard codons
        non_standard_df = df[~df['extracted_start_codon'].isin(standard_start_codons)]
        
        print(f"Total sequences with a non-standard start codon: {len(non_standard_df)}")
        
        if len(non_standard_df) > 0:
            print("\nBreakdown of what those non-standard codons are:")
            codon_counts = non_standard_df['extracted_start_codon'].value_counts()
            for codon, count in codon_counts.items():
                display_name = "MISSING/EMPTY" if codon.strip() == "" or codon == "NAN" else codon
                print(f"  - Codon '{display_name}': {count} genes")
    else:
        print("Error: 'cds_seq' column not found in dataset to analyze start codons.")

if __name__ == "__main__":
    csv_file = "ecoli_full_tu_expression.csv"
    yaml_file = "rules copy.yaml"
    dynamic_pipeline(csv_file, yaml_file)