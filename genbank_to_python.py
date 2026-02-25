### This file parses GenBank files using RegEx and then turns it into ###
###    our own internal representation of each plasmid part.

### To use this for your own "BioPart" representation, go to the last section
###    of code and modify where it is commented to modify

import re
import os
from felix_schema import BioPart, Construct, Role

role_keywords = {
    # LOGIC: Essential for DIO/Double-Floxed constructs
    Role.recombinase_site: [
        "lox", "frt", "attp", "attb", "attl", "attr", "dre", "roxa", "recombination"
    ],
    
    # STRUCTURE: Defines viral packaging boundaries and post-transcriptional elements
    Role.structural: [
        "itr", "ltr", "packaging", "sin-ltr", "inverted terminal", "rre", "wpre", "scaffold"
    ],
    
    # ORIGIN: Replication start points (both standard and transfer)
    Role.origin: [
        "ori", "origin", "pbr322", "puc", "f1", "sv40ori", "col1", "repa", "p15a", "cole1"
    ],
    Role.origin_transfer: [
        "orit", "mob", "conjugation", "transfer"
    ],
    
    # PROMOTERS: Pol II and Pol III "On Switches"
    Role.promoter: [
        "promoter", "prm", "ptrc", "ptac", "plac", "pbad", "cmv", "u6", "cag", 
        "ef1a", "sv40", "tre", "h1", "syn", "t7", "sp6"
    ],
    Role.operator: [
        "operator", "laco", "teto", "repressor binding", "uas"
    ],
    Role.enhancer: [
        "enhancer", "enh"
    ],
    
    # TERMINATORS: Transcription stop signals
    Role.terminator: [
        "terminator", "term", "stop", "rrnb", "t0", "t1", "t7te", "transcription stop"
    ],
    Role.polyA_signal: [
        "polya", "bgh", "hgh", "sv40 poly(a)", "poly(a)"
    ],
    
    # RBS: Translation initiation
    Role.rbs: [
        "rbs", "kozak", "shine-dalgarno", "ires", "ribosome binding", "tir"
    ],
    
    # VERIFICATION & TAGS: Lab tools and visual markers
    Role.verification: [
        "primer", "bind", "sequencing", "seq", "pcr", "probe", "universal", "m13"
    ],
    Role.protein_tag: [
        "tag", "gfp", "rfp", "yfp", "cfp", "mcherry", "his", "flag", "myc", "ha", 
        "v5", "gst", "mbp", "luciferase", "nlux"
    ],
    
    # SPECIALIZED RNA & GENOMIC ELEMENTS
    Role.non_coding_rna: [
        "sgrna", "crrna", "tracr", "mirna", "shrna", "trna", "rrna", "guide rna"
    ],
    Role.insulator: [
        "insulator", "chs4", "gypsy", "boundary element"
    ],
    Role.intron: [
        "intron", "ivs", "splice"
    ],
    Role.scar: [
        "scar", "assembl", "golden gate", "rfc", "biobrick"
    ],
    
    # MAINTENANCE: Lab survival genes
    Role.selection_marker: [
        "ampr", "kanr", "neor", "puro", "bla", "resistance", "hygro", "zeo", 
        "spec", "strep", "cat", "chloramphenicol", "tetr", "ampicillin", "kanamycin"
    ]
}

def get_smart_role(f_type, label):
    f_type_lower = f_type.lower()
    label_clean = label.lower()
    
    # RULE 1: Strict GenBank types that don't need keyword scanning
    if f_type_lower == 'primer_bind':
        return Role.verification
    if f_type_lower == 'terminator':
        return Role.terminator

    # RULE 2: Fuzzy keywords for ambiguous features 
    # (We do this BEFORE Rule 3 so a CDS labeled "AmpR" gets caught as a selection_marker, not just a generic CDS)
    for role, keywords in role_keywords.items():
        if any(key in label_clean for key in keywords):
            return role

    # RULE 3: Broad GenBank Type fallbacks (if keywords didn't catch it)
    type_map = {
        'cds': Role.cds,
        'rep_origin': Role.origin,
        'repeat_region': Role.repeat_region, 
        'polya_signal': Role.polyA_signal
    }
    if f_type_lower in type_map:
        return type_map[f_type_lower]

    # RULE 4: Safe Catch-all 
    # Changed from Role.selection_marker to prevent fatal assembly errors!
    return Role.unknown

def parse_genbank(file_path=None, mock_content=None):
    content = mock_content # using mock_content for testing
    if file_path and not mock_content:
        with open(file_path, 'r') as f:
            content = f.read()

    # extract name and length from LOCUS for Lean4 (no hyphens)
    locus_match = re.search(r'LOCUS\s+(\S+)', content)
    raw_name = locus_match.group(1) if locus_match else os.path.basename(file_path)
    # make it appropriate for Lean syntax
    name = re.sub(r'\W+', '_', raw_name).strip('_')

    # extract full DNA sequence from ORIGIN
    origin_start = content.find("ORIGIN")
    if origin_start == -1: return None
    # cleanup
    seq_block = content[origin_start:].replace("ORIGIN", "").replace("//", "")
    full_seq = re.sub(r'[^acgtACGT]', '', seq_block).upper()
    # initialize resulting construct
    construct = Construct({'name': name, 'intent': 'GenBank imported construct'})
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

            # structured metadata string for semantic / LLM analysis
            metadata = f"""
                feature_type: {f_type}
                location: {loc}
                orientation: {orientation}
                label: {final_label}
                """
            
            # allows us to see boundary
            part_id = f"{name}_{start}_{end}"


            ### THIS IS WHERE YOU WOULD MODIFY TO MAKE THIS RELEVANT TO YOUR OWN BIOPART ###
            ### please reach out to me if you want to work on this :)))

            role = get_smart_role(f_type, final_label)
            part = BioPart(
                part_id=part_id,
                name=final_label,
                roles=[role.name],
                left_oh=l_oh, 
                right_oh=r_oh,
                sequence=p_seq,
                metadata=metadata,
                confidence_score=0.0
            )
            construct.add_part(part, orientation=orientation)

    return construct