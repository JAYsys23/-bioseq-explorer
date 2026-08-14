"""Small, dependency-free utilities for educational DNA sequence analysis."""

from collections import Counter
import re

DNA_BASES = set("ACGTN")
COMPLEMENT = str.maketrans("ACGTN", "TGCAN")

# Standard genetic code (DNA codons). Stop codons are represented by '*'.
CODON_TABLE = {
    "TTT":"F", "TTC":"F", "TTA":"L", "TTG":"L", "CTT":"L", "CTC":"L", "CTA":"L", "CTG":"L",
    "ATT":"I", "ATC":"I", "ATA":"I", "ATG":"M", "GTT":"V", "GTC":"V", "GTA":"V", "GTG":"V",
    "TCT":"S", "TCC":"S", "TCA":"S", "TCG":"S", "CCT":"P", "CCC":"P", "CCA":"P", "CCG":"P",
    "ACT":"T", "ACC":"T", "ACA":"T", "ACG":"T", "GCT":"A", "GCC":"A", "GCA":"A", "GCG":"A",
    "TAT":"Y", "TAC":"Y", "TAA":"*", "TAG":"*", "CAT":"H", "CAC":"H", "CAA":"Q", "CAG":"Q",
    "AAT":"N", "AAC":"N", "AAA":"K", "AAG":"K", "GAT":"D", "GAC":"D", "GAA":"E", "GAG":"E",
    "TGT":"C", "TGC":"C", "TGA":"*", "TGG":"W", "CGT":"R", "CGC":"R", "CGA":"R", "CGG":"R",
    "AGT":"S", "AGC":"S", "AGA":"R", "AGG":"R", "GGT":"G", "GGC":"G", "GGA":"G", "GGG":"G",
}

# Recognition sites for a small set of widely used restriction enzymes.
RESTRICTION_ENZYMES = {
    "EcoRI": "GAATTC",
    "BamHI": "GGATCC",
    "HindIII": "AAGCTT",
    "HaeIII": "GGCC",
}


def clean_sequence(sequence: str) -> str:
    """Return an uppercase sequence with whitespace and FASTA headers removed."""
    body = (line for line in sequence.splitlines() if not line.lstrip().startswith(">"))
    return re.sub(r"\s+", "", "".join(body)).upper()


def invalid_bases(sequence: str) -> list[str]:
    return sorted(set(sequence) - DNA_BASES)


def gc_content(sequence: str) -> float:
    """GC percentage, excluding ambiguous N bases from the denominator."""
    concrete = sum(sequence.count(base) for base in "ACGT")
    return (100 * (sequence.count("G") + sequence.count("C")) / concrete) if concrete else 0.0


def melting_temperature(sequence: str) -> float:
    """Estimate DNA melting temperature using the Wallace rule for short oligos."""
    return 2 * (sequence.count("A") + sequence.count("T")) + 4 * (sequence.count("G") + sequence.count("C"))


def reverse_complement(sequence: str) -> str:
    return sequence.translate(COMPLEMENT)[::-1]


def transcribe(sequence: str) -> str:
    return sequence.replace("T", "U")


def translate(sequence: str, frame: int = 0) -> str:
    """Translate a DNA sequence from a zero-based reading frame; unknown codons are X."""
    codons = (sequence[index:index + 3] for index in range(frame, len(sequence) - 2, 3))
    return "".join(CODON_TABLE.get(codon, "X") for codon in codons)


def codon_counts(sequence: str, frame: int = 0) -> Counter:
    return Counter(sequence[index:index + 3] for index in range(frame, len(sequence) - 2, 3))


def find_orfs(sequence: str, minimum_aa: int = 10) -> list[dict]:
    """Find forward-strand ATG-to-stop ORFs in all three reading frames."""
    found = []
    for frame in range(3):
        for start in range(frame, len(sequence) - 2, 3):
            if sequence[start:start + 3] != "ATG":
                continue
            for end in range(start + 3, len(sequence) - 2, 3):
                if sequence[end:end + 3] in {"TAA", "TAG", "TGA"}:
                    protein = translate(sequence[start:end + 3])
                    if len(protein) - 1 >= minimum_aa:
                        found.append({
                            "Frame": frame + 1,
                            "Start (1-based)": start + 1,
                            "End (1-based)": end + 3,
                            "Length (aa)": len(protein) - 1,
                            "Protein": protein,
                        })
                    break
    return found


def restriction_sites(sequence: str) -> list[dict]:
    """Return one-based positions for selected restriction-enzyme recognition sites."""
    sites = []
    for enzyme, recognition_site in RESTRICTION_ENZYMES.items():
        start = 0
        while True:
            position = sequence.find(recognition_site, start)
            if position == -1:
                break
            sites.append({
                "Enzyme": enzyme,
                "Recognition site": recognition_site,
                "Position (1-based)": position + 1,
            })
            start = position + 1
    return sorted(sites, key=lambda item: item["Position (1-based)"])
