class BioPart:
    def __init__(self, part_id, name, roles, metadata=None):
        self.part_id = part_id  # e.g., "U_WPRE"
        self.name = name        # e.g., "WPRE"
        self.roles = roles      # e.g., ["utility", "rna_stabilizer"...]
        self.metadata = metadata or {}

    @property # shorthand
    def primary_type(self):
        return self.roles[0] if self.roles else "unknown"

    def __repr__(self):
        return f"<{self.name} ({self.primary_type}) | Roles: {self.roles[1:]}>"