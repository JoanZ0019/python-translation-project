#!/usr/bin/env
import sys
def translate_sequence(rna_sequence, genetic_code):
    protein = ""
    if len(rna_sequence) >= 3:
        for i in range(0, len(rna_sequence), 3):
            codon = rna_sequence[i:i+3]
            codon = codon.upper()
            if len(codon) ==3:
                if codon == 'UAA' or codon == 'UGA' or codon == 'UAG':
                    break
                else:
                    protein += genetic_code[codon]
            elif len(codon) < 3:
                break
        
    return protein

    """Translates a sequence of RNA into a sequence of amino acids.

    Translates `rna_sequence` into string of amino acids, according to the
    `genetic_code` given as a dict. Translation begins at the first position of
    the `rna_sequence` and continues until the first stop codon is encountered
    or the end of `rna_sequence` is reached.

    If `rna_sequence` is less than 3 bases long, or starts with a stop codon,
    an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the translated amino acids.
    """

def get_all_translations(rna_sequence, genetic_code):
    """Get a list of all amino acid sequences encoded by an RNA sequence.

    All three reading frames of `rna_sequence` are scanned from 'left' to
    'right', and the generation of a sequence of amino acids is started
    whenever the start codon 'AUG' is found. The `rna_sequence` is assumed to
    be in the correct orientation (i.e., no reverse and/or complement of the
    sequence is explored).

    The function returns a list of all possible amino acid sequences that
    are encoded by `rna_sequence`.

    If no amino acids can be translated from `rna_sequence`, an empty list is
    returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    list
        A list of strings; each string is an sequence of amino acids encoded by
        `rna_sequence`.
    """
    rna_sequence = rna_sequence.upper()
    if len(rna_sequence) < 3:
        return []
    amino_acid_seq_list = []
    
    for i in range(len(rna_sequence)-2):
        codon = rna_sequence[i: i + 3]
        if codon == "AUG":
            aa_seq = translate_sequence(
                    rna_sequence = rna_sequence[i:],
                    genetic_code = genetic_code)
            amino_acid_seq_list.append(aa_seq)
    return amino_acid_seq_list

    pass

def get_reverse(sequence):
    """Reverse orientation of `sequence`.

    Returns a string with `sequence` in the reverse order.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_reverse('AUGC')
    'CGUA'
    """
    if len(sequence)>=1:
        seq = sequence.upper()
        reversed_seq = seq[::-1]
        return reversed_seq
    else:
        return ''
    
    pass

def get_complement(sequence):
    """Get the complement of a `sequence` of nucleotides.

    Returns a string with the complementary sequence of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> get_complement('AUGC')
    'UACG'
    """
    if len(sequence)>=1:
        seq = list(sequence.upper())
        complement = {'A':'U', 'C':'G', 'U':'A', 'G':'C'}
        seq=[complement[base] for base in seq]
        return ''.join(seq)
    else:
        return ''

    pass

def reverse_and_complement(sequence):
    """Get the reversed and complemented form of a `sequence` of nucleotides.

    Returns a string that is the reversed and complemented sequence
    of `sequence`.

    If `sequence` is empty, an empty string is returned.

    Examples
    --------
    >>> reverse_and_complement('AUGC')
    'GCAU'
    """
    if len(sequence)>=1:
        reversed_seq = get_reverse(sequence)
        reversed_complement_seq = get_complement(reversed_seq)
        return reversed_complement_seq
    else:
        return''

def get_longest_peptide(rna_sequence, genetic_code):
    """Get the longest peptide encoded by an RNA sequence.

    Explore six reading frames of `rna_sequence` (the three reading frames of
    `rna_sequence`, and the three reading frames of the reverse and complement
    of `rna_sequence`) and return (as a string) the longest sequence of amino
    acids that it encodes, according to the `genetic_code`.

    If no amino acids can be translated from `rna_sequence` nor its reverse and
    complement, an empty string is returned.

    Parameters
    ----------
    rna_sequence : str
        A string representing an RNA sequence (upper or lower-case).

    genetic_code : dict
        A dictionary mapping all 64 codons (strings of three RNA bases) to
        amino acids (string of single-letter amino acid abbreviation). Stop
        codons should be represented with asterisks ('*').

    Returns
    -------
    str
        A string of the longest sequence of amino acids encoded by
        `rna_sequence`.
    """
    rna_sequence = rna_sequence.upper()
    start_pos = 0
    longest = ""
    amino_acids = []

    def translated(start_pos, rna_sequence, genetic_code):
        amino_acids = ""
        for i in range(start_pos, len(rna_sequence), 3):
            codon = rna_sequence[i:i + 3]
            if codon in ["UAG", "UAA", "UGA"] or len(codon) != 3:
                break
            else:
                amino_acids += genetic_code[codon]
        return amino_acids
    
    def valid_seqs(start_pos, rna_sequence, genetic_code, amino_acids):
        while start_pos < len(rna_sequence):
            start_codon = rna_sequence[start_pos:start_pos + 3]
            if start_codon == "AUG":
                amino_acid = translated(start_pos, rna_sequence, genetic_code)
                amino_acids.append(amino_acid)
            start_pos += 1
        return amino_acids

    rc_sequence = reverse_and_complement(rna_sequence)
    amino_acids = valid_seqs(start_pos, rna_sequence, genetic_code, amino_acids)
    amino_acids = valid_seqs(start_pos, rc_sequence, genetic_code, amino_acids)
    max_len = 0
    for seq in amino_acids:
        if len(seq) > max_len:
            max_len = len(seq)
            longest = seq
    return longest



if __name__ == '__main__':
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    rna_seq = ("AUG"
            "UAC"
            "UGG"
            "CAC"
            "GCU"
            "ACU"
            "GCU"
            "CCA"
            "UAU"
            "ACU"
            "CAC"
            "CAG"
            "AAU"
            "AUC"
            "AGU"
            "ACA"
            "GCG")
    longest_peptide = get_longest_peptide(rna_sequence = rna_seq,
            genetic_code = genetic_code)
    assert isinstance(longest_peptide, str), "Oops: the longest peptide is {0}, not a string".format(longest_peptide)
    message = "The longest peptide encoded by\n\t'{0}'\nis\n\t'{1}'\n".format(
            rna_seq,
            longest_peptide)
    sys.stdout.write(message)
    if longest_peptide == "MYWHATAPYTHQNISTA":
        sys.stdout.write("Indeed.\n")
