from constants import integer_mass_table, reversed_integer_mass_table

def get_subpeptides(peptide: str) -> list[str]:
    l = len(peptide)
    subpeptides = [peptide]
    cyclic_peptide = peptide * 2
    for start in range(l):
        for length in range(1, l):
            subpeptides.append(cyclic_peptide[start:start+length])
    return sorted(subpeptides, key=lambda x: len(x))

def theoretical_spectrum(peptide: str) -> list[int]:
    tspectrum = [0]
    IMT = integer_mass_table()
    subpeptides = get_subpeptides(peptide)
    for peptide in subpeptides:
        mass =  0
        for aa in peptide:
            mass += IMT[aa]
        tspectrum.append(mass)
    return sorted(tspectrum)

def extend_peptide(peptide):
    IMT = integer_mass_table()
    return [peptide + aa for aa in IMT.keys()]

def compute_mass(peptide: str) -> int:
    IMT = integer_mass_table()
    mass = 0
    for aa in peptide:
        mass += IMT[aa]
    return mass

def is_consistent(subspectrum: list[int], spectrum: list[int]):
    for mass in subspectrum:
        if mass not in spectrum:
            return False
    return True


from itertools import accumulate
def cyclopeptide_sequencing(spectrum: list[int]) -> list[str]:
    mass_list = [
        57, 71, 87, 97, 99, 101, 103, 113, 114, 115,
        128, 129, 131, 137, 147, 156, 163, 186
    ]
    acceptable_aa = [m for m in mass_list if m in spectrum]
    candidates = [[]]
    final = []
    while candidates:
        _candidates = []
        for m in acceptable_aa:
            for c in candidates:
                curr_c = c + [m]
                if sum(curr_c) in spectrum:
                    if cyclic_spectrum(curr_c) == spectrum:
                        final.append(curr_c)
                    else:
                        _candidates.append(curr_c)
        candidates = _candidates
    return final

def cyclic_spectrum(peptide: list) -> list:
    prefix = [0] + list(accumulate(peptide))
    cyclic = [0]
    for (i, n) in enumerate(prefix[:-1]):
        for (j, m) in enumerate(prefix[i+1:], start=i+1):
            cyclic.append(m-n)
            if i > 0 and j < len(peptide):
                cyclic.append(prefix[-1] - (prefix[j] - prefix[i]))
    return sorted(cyclic)

spectrum = '0 113 128 186 241 299 314 427'.split(' ')
spectrum = list(map(int, spectrum))

print(cyclopeptide_sequencing(spectrum))

Alphabet = {57: 'G', 71: 'A', 87: 'S', 97: 'P',
            99: 'V', 101: 'T', 103: 'C', 113:'I/L',
            114: 'N', 115: 'D', 128: 'K/Q', 129: 'E',
            131: 'M', 137: 'H', 147: 'F', 156: 'R', 163: 'Y', 186: 'W'}

def CountingMass(Mass, masslist):
    if Mass == 0: return 1, masslist
    if Mass < 57: return 0, masslist
    if Mass in masslist: return masslist[Mass], masslist
    n = 0
    for i in Alphabet:
        k, masslist = CountingMass(Mass - i, masslist)
        n += k
    masslist[Mass] = n
    return n, masslist


