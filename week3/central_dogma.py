from constants import rna_genetic_code

RNA_GENETIC_CODE = rna_genetic_code()

def complement(pattern: str) -> str:
    """Return the complement of a DNA string."""
    COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join([COMPLEMENT_MAP[s] for s in pattern])

def reverse_complement(pattern: str) -> str:
    """Return the reverse complement of a DNA string."""
    COMPLEMENT_MAP = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    seq = ''
    for s in pattern:
        seq = COMPLEMENT_MAP[s] + seq
    return seq

def transcribe(pattern: str) -> str:
    """Return an RNA string from a given DNA string."""
    return pattern.replace('T', 'U')

def translate(pattern: str, genetic_code: dict[str, str]) -> str:
    """Return a peptide sequence from a given RNA string."""
    if 'T' in pattern:
        pattern = pattern.replace('T', 'U')
    peptide = ''
    for i in range(0, len(pattern), 3):
        codon = pattern[i:i+3]
        if genetic_code[codon] == "*":
            return peptide
        peptide += genetic_code[codon]
    return peptide

def generate_reading_frames(pattern: str) -> list[str]:
    frames = []
    for i in range(3):
        frames.append(pattern[i:])
    return [f[:len(f)-(len(f)%3)] for f in frames]

def find_coding_region(pattern: str, peptide: str, genetic_code: dict[str, str] = RNA_GENETIC_CODE) -> list[str]:
    coding_regions = []
    # convert DNA to RNA
    l = len(peptide) * 3
    for i in range(len(pattern)-l+1):
        seq = pattern[i:i+l]
        rc = reverse_complement(seq)
        seq_pep = translate(seq, genetic_code)
        rc_pep = translate(rc, genetic_code)
        if seq_pep == peptide or rc_pep == peptide:
            coding_regions.append(seq)
    return coding_regions

def main() -> None:
    pattern = 'ATGGCCATGGCCCCCAGAACTGAGATCAATAGTACCCGTATTAACGGGTGA'
    peptide = 'MA'
    result = find_coding_region(pattern=pattern, peptide=peptide)
    print(result)

if __name__ == '__main__':
    main()
