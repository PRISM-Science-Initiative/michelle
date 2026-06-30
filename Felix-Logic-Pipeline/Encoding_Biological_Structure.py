# this script contains classes and functions for encoding biological structure information
# as described in Dr. Naik's summer plan document.

class CircularSequence:
    # this class represents a circular sequence, such as a plasmid, using the duplicated string technique

    def __init__(self, sequence_string):
        self.seq = sequence_string.upper()
        self.N = len(self.seq)
        self.duplicated_sequence_string = self.seq + self.seq

    def display_info(self):
        print(f"Sequence: {self.seq}")
        print(f"Length: {self.N}")
        print(f"Duplicated Sequence: {self.duplicated_sequence_string}")
    

    def get_slice(self, start_index, length):
        # returns a slice of the circular sequence from start_index of length 'length'
        if length > self.N:
            raise ValueError("Requested length exceeds sequence length.")
        end_index = start_index + length
        return self.duplicated_sequence_string[start_index:end_index]
    
    def find_feature(self, feature_sequence):
        # returns the start index of the feature_sequence in the circular sequence, or -1 if not found
        feature_sequence = feature_sequence.upper()
        index = self.duplicated_sequence_string.find(feature_sequence)
        if index == -1 or index >= self.N:
            return -1
        return index
        # feature to think about adding: flag if it wraps around the origin (i.e., if the feature starts near the end of the sequence and continues at the beginning)
    
    def ordinary_overlap(self, linear_sequence):
        # returns true if the linear_sequence is found in the circular sequence without wrapping around the origin, false otherwise
        linear_sequence = linear_sequence.upper()
        index = self.seq.find(linear_sequence)
        if index == -1:
            return False
        return True
    
    def circular_overlap(self, linear_sequence):
        # returns true if the linear_sequence spans accross the junction where the sequence loops back on itself, false otherwise
        linear_sequence = linear_sequence.upper()
        index = self.duplicated_sequence_string.find(linear_sequence)
        if index == -1 or index >= self.N:
            return False
        else: 
            end_index = index + len(linear_sequence)
            if end_index > self.N:
                return True
            else:
                return False
    
    def has_internal_cut_site(self, enzyme_recognition_sequences):
        # returns true if the enzyme recognition sequence is found in the circular sequence, false otherwise
        for enzyme_seq in enzyme_recognition_sequences:
            if self.find_feature(enzyme_seq) == -1:
                return False
        return True
        

def test_circular_sequence():
    # test the CircularSequence class with a sample sequence
    seq = "ATGCGTACGTA"
    circular_seq = CircularSequence(seq)
    circular_seq.display_info()
    
    # test get_slice
    print("Slice from index 8 of length 5:", circular_seq.get_slice(8, 5))  # should wrap around
    
    # test find_feature
    print("Index of feature 'CGT':", circular_seq.find_feature("CGT"))  # should return index 3
    
    # test ordinary_overlap
    print("Ordinary overlap with 'GTA':", circular_seq.ordinary_overlap("GTA"))  # should return True
    
    # test circular_overlap
    print("Circular overlap with 'ACGTA':", circular_seq.circular_overlap("ACGTA"))  # should return True
    
    # test has_internal_cut_site
    enzymes = ["CGT", "TAC"]
    print("Has internal cut site for enzymes:", circular_seq.has_internal_cut_site(enzymes))  # should return True

test_circular_sequence()