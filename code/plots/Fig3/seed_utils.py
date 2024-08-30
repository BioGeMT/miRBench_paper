from Bio.Seq import Seq

def reverse_complement(seq):
    return str(Seq(seq).reverse_complement())

def find_seed_match(target, mirna):
    rc_target = reverse_complement(target)
    seed6 = mirna[1:7]
    seed7 = mirna[1:8]
    seed8 = 'A' + mirna[1:8]

    # Seed8mer
    if seed8 in rc_target:
        return 'Seed8mer'

    # Seed7mer
    if seed7 in rc_target or ('A' + seed6) in rc_target:
        return 'Seed7mer'

    # Seed6mer (including 6mer_m8/A1)
    if seed6 in rc_target or mirna[2:8] in rc_target or ('A' + mirna[1:6]) in rc_target:
        return 'Seed6mer'

    # SeedNonCanonical (6mer with bulge or mismatch for all 3 subcategories)
    seed_types = [seed6, mirna[2:8], 'A' + mirna[1:6]]  # 6mer, 6mer-m8, 6mer-A1
    for seed in seed_types:
        for pos in range(len(seed)):
            for nt in ['A', 'C', 'G', 'T']:
                # bulges
                if (seed[:pos] + nt + seed[pos:]) in rc_target:
                    return 'SeedNonCanonical'
                # mismatches
                if pos < len(seed) and (seed[:pos] + nt + seed[pos+1:]) in rc_target:
                    return 'SeedNonCanonical'

    return 'none'
