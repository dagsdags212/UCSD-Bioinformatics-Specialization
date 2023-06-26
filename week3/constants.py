from pathlib import Path

RNA_CODON_TABLE = Path('/home/dagsdags/home/courses/ucsd-bioinformatics-specialization/course-2/datasets/tables/RNA_codon_table_1.txt')
INTEGER_MASS_TABLE = Path('/home/dagsdags/home/courses/ucsd-bioinformatics-specialization/course-2/datasets/tables/integer_mass_table.txt')

def rna_genetic_code(path: Path = RNA_CODON_TABLE):
    RGC = {}
    with open(path, 'r') as fh:
        for line in fh:
            res = line.rstrip().split()
            if len(res) == 1:
                RGC[res[0]] = "*"
            else:
                RGC[res[0]] = res[1]
    return RGC

def integer_mass_table(path: Path = INTEGER_MASS_TABLE):
    IMT = {}
    with open(path, 'r') as fh:
        for line in fh:
            k, v = line.rstrip().split()
            IMT[k] = int(v)
    fh.close()
    return IMT

def reversed_integer_mass_table(path: Path = INTEGER_MASS_TABLE):
    IMT = {}
    with open(path, 'r') as fh:
        for line in fh:
            v, k = line.rstrip().split()
            IMT[int(k)] = v
    fh.close()
    return IMT

