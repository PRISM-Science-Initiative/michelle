import re
import os
from felix_schema import BioPart, Construct, Role

# role mapping logic
# expanded mapping dictionary for fuzzy matching labels

# THIS PART NEEDS A LOT OF WORK... mostly because I don't know biology so this was
# mostly gemini's description of what goes where

# wondering if we can use vectorization/ML/
# other smarter logical mappings to make this parser choose the correct role
role_keywords = {
    # LOGIC: Essential for DIO/Double-Floxed constructs
    Role.recombinase_site: [
        "lox", "frt", "attp", "attb", "dre", "roxa", "recomb", "protein bind", "recognition"
    ],
    
    # STRUCTURE: Defines viral packaging boundaries
    Role.structural: [
        "itr", "ltr", "packaging", "sin-ltr", "inverted terminal", "rre", "wpre"
    ],
    
    # ORIGIN: Replication start points
    Role.origin: [
        "ori", "origin", "pbr322", "puc", "f1", "sv40ori", "col1", "replication"
    ],
    
    # PROMOTERS: Pol II and Pol III "On Switches"
    Role.promoter: [
        "promoter", "enhancer", "cmv", "u6", "cag", "ef1a", "sv40", "tre", "h1", "syn"
    ],
    
    # TERMINATORS: Transcription stop signals
    Role.terminator: [
        "polya", "terminator", "stop", "hgh", "bgh", "t7 terminator", "transcription stop"
    ],
    
    # RBS: Translation initiation
    Role.rbs: [
        "rbs", "kozak", "shine-dalgarno", "ires", "ribosome"
    ],
    
    # VERIFICATION: Non-functional lab tools
    Role.verification: [
        "primer", "bind", "sequencing", "seq", "pcr", "probe", "universal"
    ],
    
    # MAINTENANCE: Lab survival genes
    Role.selection_marker: [
        "ampr", "kanr", "neor", "puro", "bla", "resistance", "hygro", "zeo"
    ]
}

def get_smart_role(f_type, label):
    label_clean = label.lower()
    
    # RULE 1: primers are always Verification
    if f_type == 'primer_bind':
        return Role.verification

    # RULE 2: standard genbank types (highest accuracy)
    type_map = {
        'CDS': Role.cds,
        'rep_origin': Role.origin,
        'repeat_region': Role.structural,
        'terminator': Role.terminator,
        'polyA_signal': Role.terminator
    }
    if f_type in type_map:
        return type_map[f_type]

    # RULE 3: Fuzzy keywords for ambiguous features (promoter, enhancer, lox sites)
    for role, keywords in role_keywords.items():
        if any(key in label_clean for key in keywords):
            return role
            
    return Role.selection_marker
# Actual parser

def parse_genbank(file_path):
    with open(file_path, 'r') as f:
        content = f.read()

    # extract name and length from LOCUS for Lean4 (no hyphens)
    locus_match = re.search(r'LOCUS\s+(\S+)', content)
    raw_name = locus_match.group(1) if locus_match else os.path.basename(file_path)
    # make it appropriate for Lean syntax
    name = raw_name.replace('-', '_').strip('_')

    # extract full DNA sequence from ORIGIN
    origin_start = content.find("ORIGIN")
    full_seq = re.sub(r'[\d\s/]+', '', content[origin_start:].replace("ORIGIN", "").replace("//", "")).upper()

    # initialize resulting construct
    construct = Construct({'name': name, 'intent': 'Extracted from GenBank'})

    # track coordinates to prevent duplicates: (start, end, orientation)
    seen_features = set()

    # process features (actually makes the parts)
    feature_section = re.search(r'FEATURES\s+Location/Qualifiers(.*?)\nORIGIN', content, re.DOTALL)  

    if feature_section:
        raw_features = feature_section.group(1)
        pattern = re.compile(r'^ {5}(\w+) +(?:(complement)\()?([\d\.\,]+|join\([\d\.\,\w \(\)]+\))\)?', re.MULTILINE)
        feature_matches = list(pattern.finditer(raw_features))

        for i, match in enumerate(feature_matches):
            f_type, is_comp, loc = match.groups()
            # skip the 'source' feature (the whole plasmid)
            if f_type == 'source':
                continue
            try:
                coords = re.findall(r'\d+', loc)
                start, end = int(coords[0]), int(coords[-1])
            except: 
                continue

            # de-duplication
            orientation = "reverse" if is_comp else "forward"
            feature_key = (start, end, orientation)
            if feature_key in seen_features: continue
            seen_features.add(feature_key)

            # extract label
            end_q = feature_matches[i+1].start() if i+1 < len(feature_matches) else len(raw_features)
            block = raw_features[match.end():end_q]
            label_match = re.search(r'/(?:label|note)="([^"]+)"', block)
            final_label = label_match.group(1).replace("\n", " ").strip() if label_match else f_type

            # overhang extraction 
            #  pull the 4bp flanking the part in the actual plasmid sequence
            left_idx = start - 1
            # handle circular wrap-around for the left overhang
            l_oh = full_seq[left_idx-4 : left_idx] if left_idx >= 4 else full_seq[-4:]
            # handle circular wrap-around for the right overhang
            r_oh = full_seq[end : end+4] if end+4 <= len(full_seq) else full_seq[:4]

            role = get_smart_role(f_type, final_label)

            p_seq = full_seq[start-1:end]
            part = BioPart(
                part_id=f"ID_{final_label[:8].replace(' ', '_')}",
                name=final_label,
                roles=[role.name], # passes the string name of the Enum
                left_oh=l_oh, 
                right_oh=r_oh,
                sequence=p_seq,
                metadata=f"Type: {f_type}, Loc: {loc}"
            )
            construct.add_part(part, orientation=orientation)

    return construct