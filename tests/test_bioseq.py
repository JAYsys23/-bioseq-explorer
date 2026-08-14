from bioseq import clean_sequence, find_orfs, gc_content, melting_temperature, restriction_sites, reverse_complement, translate


def test_clean_sequence_removes_fasta_header_and_whitespace():
    assert clean_sequence(">sequence one\n atg c\n") == "ATGC"


def test_reverse_complement():
    assert reverse_complement("ATGCN") == "NGCAT"


def test_gc_content_excludes_n():
    assert gc_content("GCCATN") == 60.0


def test_translate_standard_codon_and_stop():
    assert translate("ATGTAA") == "M*"


def test_find_orf():
    results = find_orfs("ATGAAATAA", minimum_aa=2)
    assert results[0]["Protein"] == "MK*"


def test_melting_temperature_wallace_rule():
    assert melting_temperature("ATGC") == 12


def test_restriction_site_position_is_one_based():
    sites = restriction_sites("TTGAATTCAA")
    assert sites == [{"Enzyme": "EcoRI", "Recognition site": "GAATTC", "Position (1-based)": 3}]
