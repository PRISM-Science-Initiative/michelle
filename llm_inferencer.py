from openai import OpenAI
from felix_schema import Role

# this is chopped as hell rn please don't use this, it's just framework for what
#    I will do later

class LLMInferencer:
    def __init__(self, api_key):
        self.client = OpenAI(api_key=api_key)
    
    # takes the metadata and returns a (list of roles, confidence number)
    def get_role_list(self, metadata_text):
        prompt = f"""
        Act as a computational biologist. Analyze this DNA part metadata and
        assign one or more roles from this list: {[r.name for r in Role]}.
        Metadata: {metadata_text}
        Return ONLY a comma-separated list of roles.
        """
        try:
            response = self.client.chat.completions.create(
                model="gpt-4o",
                messages=[{"role": "user", "content": prompt}]
            )
            roles = [r.strip() for r in response.choices[0].message.content.split(",")]
            return roles, 0.85
            # MODIFY THIS TO CHANGE HOW CONFIDENT WE ARE IN THE LLM'S REASONING
        except:
            return ["selection_marker"], 0.1