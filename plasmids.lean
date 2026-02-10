section
open List

inductive Role where
  | promoter
  | rbs
  | cds
  | terminator
  | recombinase_site
  | structural
  | origin
  | selection_marker
  | repeat_region
  | verification
  | non_coding_rna
  | insulator
  | intron
  -- add more roles here when we come across them
  deriving BEq, Repr, Inhabited

structure Part where
  sequence : String
  roles : List Role
  leftOverhang : String
  rightOverhang : String
  metadata : String
  deriving Repr, Inhabited

inductive Orientation where
  | forward -- 5' to 3' on reference
  | reverse -- 3' to 5' on complement
  deriving BEq, Repr, Inhabited

structure DirectedPart where
  part : Part
  orientation : Orientation
  deriving Repr, Inhabited

structure Plasmid where
  elements : List DirectedPart
  -- constraint: must have at least one part to be circular
  is_valid : elements.length > 0
  -- add the other constraints here later

def complement (c : Char) : Char :=
  match c with
  | 'A' => 'T'
  | 'T' => 'A'
  | 'C' => 'G'
  | 'G' => 'C'
  | 'a' => 't'
  | 't' => 'a'
  | 'c' => 'g'
  | 'g' => 'c'
  | _   => c

def reverseComplement (s : String) : String :=
  let charList := s.toList.reverse
  let complementedList := charList.map complement
  String.ofList complementedList

def getEffectiveRightOH (p : DirectedPart) : String :=
  match p.orientation with
  | Orientation.forward => p.part.rightOverhang
  | Orientation.reverse => reverseComplement p.part.leftOverhang

def getEffectiveLeftOH (p : DirectedPart) : String :=
  match p.orientation with
  | Orientation.forward => p.part.leftOverhang
  | Orientation.reverse => reverseComplement p.part.rightOverhang

def getNextPart (p : Plasmid) (idx : Nat) : DirectedPart :=
  -- just represent circularity with mods
  p.elements[idx % p.elements.length]!

def overhangsMatch (p1 p2 : DirectedPart) : Prop :=
  getEffectiveRightOH p1 = getEffectiveLeftOH p2

def isMoCloValid (p : Plasmid) : Prop :=
  ∀ (i : Nat), i < p.elements.length →
    let current := p.elements[i]!
    let next := getNextPart p (i+1)
    overhangsMatch current next

def circularDist (p : Plasmid) (i j : Nat) : Nat :=
  let n := p.elements.length
  (j+n-i) % n

def isUpstream (p : Plasmid) (i j : Nat) : Prop :=
  -- dist from i to j is smaller  than full revolution
  let d := circularDist p i j
  d > 0 ∧ d < p.elements.length
-- this is just an example of a constraint we could encode

def hasValidTranscriptionUnit (p : Plasmid) : Prop :=
  ∃ (i j k : Nat),
    i < p.elements.length ∧ j < p.elements.length ∧ k < p.elements.length ∧
    let promoter := p.elements[i]!
    let rbs := p.elements[j]!
    let cds := p.elements[k]!
    promoter.part.roles.contains Role.promoter ∧
    rbs.part.roles.contains Role.rbs ∧
    cds.part.roles.contains Role.cds ∧
    -- orientation and ordering constraints
    promoter.orientation = rbs.orientation ∧
    rbs.orientation = cds.orientation ∧
    isUpstream p i j ∧ isUpstream p j k
end
