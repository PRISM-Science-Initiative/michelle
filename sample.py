# Create your construct
p = Construct({'name': 'MyExperiment'})

# Add a CDS first
p.add_part(gfp_part, orientation="forward")

# Now 'retroactively' add a promoter at the beginning (index 0) 
# because a constraint requires it
p.add_part(promoter_part, orientation="forward", index=0)

# Check the layout
print(p) # Output: Construct: MyExperiment | Layout: [pCMV(F) -> GFP(F)]