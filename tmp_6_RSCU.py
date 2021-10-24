from Bio import SeqIO
from itertools import chain
import Bio.Data.CodonTable as ct
from collections import Counter

# modified based on https://github.com/gavieira/python_bioinfo/blob/1e1cf7786c75f3ce87ac9af94f176cd5e6b254dd/codon_usage/CAI.py


def _synonymous_codons(genetic_code_dict):

    # invert the genetic code dictionary to map each amino acid to its codons
    codons_for_amino_acid = {}
    for codon, amino_acid in genetic_code_dict.items():
        codons_for_amino_acid[amino_acid] = codons_for_amino_acid.get(amino_acid, [])
        codons_for_amino_acid[amino_acid].append(codon)

    # create dictionary of synonymous codons
    # Example: {'CTT': ['CTT', 'CTG', 'CTA', 'CTC', 'TTA', 'TTG'], 'ATG': ['ATG']...}
    return {
        codon: codons_for_amino_acid[genetic_code_dict[codon]]
        for codon in genetic_code_dict.keys()
    }


_synonymous_codons = {k: _synonymous_codons(v.forward_table) for k, v in ct.unambiguous_dna_by_id.items()}
_non_synonymous_codons = {k: {codon for codon in v.keys() if len(v[codon]) == 1} for k, v in _synonymous_codons.items()}


def RSCU(sequences, genetic_code=11):
    r"""Calculates the relative synonymous codon usage (RSCU) for a set of sequences.

    RSCU is 'the observed frequency of [a] codon divided by the frequency
    expected under the assumption of equal usage of the synonymous codons for an
    amino acid' (page 1283).

    In math terms, it is

    .. math::

        \frac{X_{ij}}{\frac{1}{n_i}\sum_{j=1}^{n_i}x_{ij}}

    "where :math:`X` is the number of occurrences of the :math:`j` th codon for
    the :math:`i` th amino acid, and :math:`n` is the number (from one to six)
    of alternative codons for the :math:`i` th amino acid" (page 1283).

    Args:
        sequences (list): The reference set of sequences.
        genetic_code (int, optional): The translation table to use. Defaults to 11, the standard genetic code.

    Returns:
        dict: The relative synonymous codon usage.

    Raises:
        ValueError: When an invalid sequence is provided or a list is not provided.
    """

    if not isinstance(sequences, (list, tuple)):
        raise ValueError(
            "Be sure to pass a list of sequences, not a single sequence. "
            "To find the RSCU of a single sequence, pass it as a one element list."
        )

    # ensure all input sequences are divisible by three
    for sequence in sequences:
        if len(sequence) % 3 != 0:
            raise ValueError("Input sequence not divisible by three")
        if not sequence:
            raise ValueError("Input sequence cannot be empty")

    # count the number of each codon in the sequences
    sequences = (
        (sequence[i : i + 3].upper() for i in range(0, len(sequence), 3))
        for sequence in sequences
    )
    codons = chain.from_iterable(
        sequences
    )  # flat list of all codons (to be used for counting)
    counts = Counter(codons)

    # "if a certain codon is never used in the reference set... assign [its
    # count] a value of 0.5" (page 1285)
    for codon in ct.unambiguous_dna_by_id[genetic_code].forward_table:
        if counts[codon] == 0:
            counts[codon] = 0.5

    # determine the synonymous codons for the genetic code
    synonymous_codons = _synonymous_codons[genetic_code]

    # hold the result as it is being calulated
    result = {}

    # calculate RSCU values
    for codon in ct.unambiguous_dna_by_id[genetic_code].forward_table:
        result[codon] = counts[codon] / (
            (len(synonymous_codons[codon]) ** -1)
            * (sum((counts[_codon] for _codon in synonymous_codons[codon])))
        )

    return result


gene_1_nc_seq = 'ATGGCTCGCGCGACCCAGGGTTCGAACAATGTCTGGGTGATCTGGCTAAGCCTGCTGGTCGGGCTGGTGCTGGCCGTAGCGCCGATGCCTACCTTCACCGAGATCGGCCGCCCGCTGTGGCTGGCGATGCTGCTGACCTACTGGGTGTTGCTGCTGCCCGAGCGCGTCGGGATGGTCACCGCCTGGCTGCTCGGCCTGGCCCAGGATGTGCTCTACGGCAACCTGCTCGGGCAGAACGCGCTGATCCTCGGCCTGATCACCTTCCTCGTGCTGTCCCTGCACCAGCGCCTGCGCATGTTCCCGGCCTGGCAGCAGTGCCTGGTGCTGGTGGTGGTCTACGGCCTGGCGCAGCTGCTCCAGCTCTGGCTCAACGCCCTGACCGGCAATCGCCCGCCGACCCTGACCTTCCTGCTTCCGGCCCTGGTCAGCGCGCTGCTCTGGCCCTGGGTGTACGCCCTGCTGCAGTTCGTCCGTCTGCGTCTCAACGTGAGATAG'
gene_2_nc_seq = 'ATGGCTCGCGCGACCCAGGGTTCGAACAATGTCTGGGTGATCTGGCTGAGCCTGCTGGTCGGGCTGGTGCTGGCCGTAGCGCCGATGCCTACCTTCACCGAGATCGGCCGCCCGCTGTGGCTGGCGATGCTGCTGACCTACTGGGTGTTGCTGCTGCCCGAGCGCGTCGGGATGGTCACCGCCTGGCTGCTCGGCCTGGCCCAGGATGTGCTCTACGGCAACCTGCTCGGGCAGAACGCGCTGATCCTCGGCCTGATCACCTTCCTCGTGCTGTCCCTGCACCAGCGCCTGCGCATGTTCCCGGCCTGGCAGCAGTGCCTGGTGCTGGTGGTGGTCTACGGCCTGGCGCAGCTGCTCCAGCTCTGGCTCAACGCCCTGACCGGCAATCGCCCGCCGACCCTGACCTTCCTGCTTCCGGCCCTGGTCAGCGCGCTGCTCTGGCCCTGGGTGTACGCCCTGCTGCAGTTCGTCCGTCTGCGTCTCAACGTGAGATAG'
gene_1_RSCU_dict = RSCU([gene_1_nc_seq])
gene_2_RSCU_dict = RSCU([gene_2_nc_seq])
print(gene_1_RSCU_dict['CTA'])
print(gene_2_RSCU_dict['CTG'])
print(gene_1_nc_seq[47])

