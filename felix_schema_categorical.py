class BioPart:
    def __init__(self, part_id, name, category, roles=None, metadata=None):
        self.part_id = part_id # i.e. P_CAG
        self.name = name # i.e. CAG Promoter
        self.category = category # i.e. promoter
        self.roles = roles if roles else []
        self.metadata = metadata if metadata else {}
    
    def has_role(self, role_name):
        return role_name in self.roles
        
    def __repr__(self):
        return f"<{self.part_id}:{self.name} | Roles: {self.roles}>"
    
class Construct:
    def __init__(self, registry_info):
        self.name = registry_info['name']
        self.addgene_id = registry_info['id']
        self.intent = registry_info['intent']
        self.parts_layout = [] # change this later

    