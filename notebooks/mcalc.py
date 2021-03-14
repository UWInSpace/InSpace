def get_prot_mass(AASeq):
    
    '''Tabulate approximate total mass of input amino acid string based on the sum of
    the average masses for every amino acid in the string.'''
    
    # Create dict with amino acid sybols as keys and associated average masses as values
    # Average masses taken from http://proteomicsresource.washington.edu/protocols06/masses.php
    molecular_mass_dict = {
    'A':71.0779, 
    'R':156.18568, 
    'G':57.05132, 
    'S':87.0773, 
    'P':97.11518, 
    'V':99.13106, 
    'T':101.10388, 
    'C':103.1429, 
    'L':113.15764, 
    'I':113.15764, 
    'N':114.10264, 
    'D':115.0874, 
    'Q':128.12922, 
    'K':128.17228, 
    'E':129.11398, 
    'M':131.19606, 
    'H':137.13928, 
    'F':147.17386, 
    'U':150.3079, 
    'Y':163.17326, 
    'W':186.2099, 
    'O':237.29816
    }
    
    # Calculate approx mass for amino acid sequence
    mass_list = []
    
    # Capitalize all letters in amino acid string prior to taking mass
    AASeq = AASeq.upper()
    
    # Raise error if amino acid sequence contains characters not in the dict
    acceptable_amino_acids_list = list(molecular_mass_dict.keys())

    for char in AASeq:
        if char not in acceptable_amino_acids_list:
            if char == ' ':
                continue
            else:
                raise ValueError('AASeq contains unacceptable characters')
    
    for amino_acid, mw in molecular_mass_dict.items():
        aa_weight = AASeq.count(amino_acid)*mw
        mass_list.append(aa_weight)
    total_mass = sum(mass_list)
    
    return total_mass