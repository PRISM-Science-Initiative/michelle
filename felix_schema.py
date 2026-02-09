# this class represents the "lego brick" functional chunks of dna
class BioPart:
    # complement table for base pairs
    _COMPLEMENT_TABLE = str.maketrans("ATCGatcg", "TAGCtagc")

    def __init__(self, part_id, name, roles, left_oh, right_oh, sequence, metadata):
        self.part_id = part_id  # e.g., "U_WPRE"
        self.name = name        # e.g., "WPRE"
        self.sequence = sequence
        self.roles = roles      # e.g., ["promoter", "rna_stabilizer"...]
        self.left_oh = left_oh # overhang
        self.right_oh = right_oh 
        self.metadata = metadata 

    @staticmethod 
    def rc(seq):
    # static utility string -> string reverse complement
        return seq.translate(BioPart._COMPLEMENT_TABLE)[::-1].upper()

    def reverse_complement(self):
    # returns reverse complement of this part
        return self.rc(self.sequence)

    @property
    def primary_type(self):
        return self.roles[0] if self.roles else "unknown"

    def __repr__(self):
        return f"<{self.name} ({self.primary_type}) | Roles: {self.roles[1:]}>"

# this class represents the entire construct made out of chunks
class Construct:
    def __init__(self, registry_info, is_circular=True):
        self.name = registry_info.get('name', 'Unnamed_Construct')
        self.addgene_id = registry_info.get('id')
        self.intent = registry_info.get('intent', 'general_assembly')
        self.is_circular = is_circular
        # parts_layout stores tuples of (BioPart, orientation)
        self.parts_layout = []

    def add_part(self, part: BioPart, orientation: str = "forward", index=None):
    # allows us to insert parts into our construct
    # input: part to insert, orientation (forward or backwards), and index to insert at
        if orientation not in ["forward", "reverse"]:
            raise ValueError("Orientation must be 'forward' or 'reverse'")
        if index is None:
            self.parts_layout.append((part,orientation))
        else:
            self.parts_layout.insert(index, (part, orientation))

    def _get_active_overhangs(self, part: BioPart, orientation):
    # calculates effective overhangs based on orientation
        if orientation == "forward":
            return part.left_oh, part.right_oh
        else:
            # reverse orientation flips part: left is reverse comp of right, 
            # right is reverse comp of left
            return part.rc(part.right_oh), part.rc(part.left_oh)

    def find_incompatibility_gaps(self):
    # returns a list of indices where part overhangs don't match
        gaps = []
        num_parts = len(self.parts_layout)
        if num_parts < 2:
            return gaps
        # check internal overlaps
        for i in range(num_parts - 1):
            j = (i+1) % num_parts
            if not self.is_circular and i == num_parts - 1:
                break
            part_a, orient_a = self.parts_layout[i]
            part_b, orient_b = self.parts_layout[i+1]
            
            _, effective_right_a = self._get_active_overhangs(part_a, orient_a)
            effective_left_b, _ = self.get_active_overhangs(part_b, orient_b)

            if effective_right_a != effective_left_b:
                gaps.append({
                    "junction": (i, j),
                    "expected": effective_right_a,
                    "actual": effective_left_b,
                    "part_a_name": part_a.name,
                    "part_b_name": part_b.name
                })

        return gaps

    def suggest_patches(self):
    # this is thinking ahead to the logic behind fixing gaps
        gaps = self.find_incompatibility_gaps()
        if not gaps:
            print(f"Construct '{self.name}' is physically compatible (MoClo valid).")
            return
        for gap in gaps:
            i, j = gap["junction"]
            print(f"Gap found between {i} and {j} ({gap['part_a_name']} to {gap['part_b_name']})")
            print(f"   Current junction: {gap['expected']} | {gap['actual']}")
            print(f"   Required: A part with Left [{gap['expected']}] and Right [{gap['actual']}]")

    def get_full_sequence(self):
    # assemble sequences, using biopart reverse complement logic
        full_seq = []
        for part, orientation in self.parts_layout:
            if orientation == "reverse":
                full_seq.append(part.reverse_complement())
            else:
                full_seq.append(part.sequence)
        return "".join(full_seq)
    
    def to_lean_definition(self):
        lean_parts = []
        for part, orientation in self.parts_layout:
            roles_list = ", ".join([f"Role.{r}" for r in part.roles])
            # if we end up changing lean implementation, we have to go back and
            # change what's here
            part_def = (
                f"{{ part := {{ "
                f"sequence := \"{part.sequence}\", "
                f"roles := [{roles_list}], "
                f"leftOverhang := \"{part.left_oh}\", "
                f"rightOverhang := \"{part.right_oh}\", "
                f"metadata := \"{part.metadata}\" }}, "
                f"orientation := Orientation.{orientation} }}"
            )
            lean_parts.append(part_def)
        parts_list_str = ",\n    ".join(lean_parts)
        return (f"def {self.name.replace(' ', '_')} : Plasmid := {{\n"
            f"  elements := [\n    {parts_list_str}\n  ],\n"
            f"  is_valid := by simp; exact Nat.zero_lt_succ _ \n"
            f"}}")
    
    def __repr__(self):
        layout_str = " -> ".join([f"{p.name}({o[0].upper()})" for p, o in self.parts_layout])
        return f"Construct: {self.name} | Layout: [{layout_str}]"