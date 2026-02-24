### This file parses GenBank files using RegEx and then turns it into ###
###    our own internal representation of each plasmid part.

### To use this for your own "BioPart" representation, go to the last section
###    of code and modify where it is commented to modify

import re
import os
from felix_schema import BioPart, Construct, Role

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
    if origin_start == -1: return None
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

            # extract block of metadata for this feature
            # if you want to change what metadata you extract, here is where you would do it
            next_start = feature_matches[i+1].start() if i+1 < len(feature_matches) else len(raw_features)
            block = raw_features[match.end():next_start]
            label_match = re.search(r'/(?:label|note|gene)="([^"]+)"', block)
            final_label = label_match.group(1).replace("\n", " ").strip() if label_match else f_type
            p_seq = full_seq[start-1:end]
            # overhang extraction 
            #  pull the 4bp flanking the part in the actual plasmid sequence
            left_idx = start - 1
            # handle circular wrap-around for the left overhang
            l_oh = full_seq[left_idx-4 : left_idx] if left_idx >= 4 else full_seq[-4:]
            # handle circular wrap-around for the right overhang
            r_oh = full_seq[end : end+4] if end+4 <= len(full_seq) else full_seq[:4]


            ### THIS IS WHERE YOU WOULD MODIFY TO MAKE THIS RELEVANT TO YOUR OWN BIOPART ###
            ### please reach out to me if you want to work on this :)))
            part = BioPart(
                part_id=f"ID_{final_label[:8].replace(' ', '_')}",
                name=final_label,
                roles=["unknown"],
                left_oh=l_oh, 
                right_oh=r_oh,
                sequence=p_seq,
                metadata=f"Type: {f_type}, Loc: {loc}",
                confidence_score=0.0 
            )
            construct.add_part(part, orientation=orientation)

    return construct