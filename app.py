import streamlit as st

from bioseq import (
    clean_sequence, codon_counts, find_orfs, gc_content, invalid_bases,
    melting_temperature, restriction_sites, reverse_complement, transcribe, translate,
)

st.set_page_config(page_title="BioSeq Explorer", page_icon="🧬", layout="wide")

SAMPLE = ">Example coding sequence\nATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

st.title("🧬 BioSeq Explorer")
st.caption("An educational DNA sequence analyzer built with Python and Streamlit.")

with st.sidebar:
    st.header("Analysis settings")
    frame = st.selectbox("Reading frame", [1, 2, 3], index=0)
    minimum_orf = st.slider("Minimum ORF length (amino acids)", 1, 100, 10)
    if st.button("Load example sequence"):
        st.session_state.sequence_input = SAMPLE

sequence_input = st.text_area(
    "Paste a DNA sequence or FASTA entry", value=st.session_state.get("sequence_input", SAMPLE), height=180,
    help="Allowed bases: A, C, G, T, and N. FASTA headers and whitespace are ignored.",
)
sequence = clean_sequence(sequence_input)
invalid = invalid_bases(sequence)

if not sequence:
    st.info("Paste a DNA sequence to begin.")
    st.stop()
if invalid:
    st.error(f"Invalid base(s) found: {', '.join(invalid)}. Please use only A, C, G, T, and N.")
    st.stop()

counts = {base: sequence.count(base) for base in "ACGTN"}
metrics = st.columns(4)
metrics[0].metric("Length", f"{len(sequence):,} bp")
metrics[1].metric("GC content", f"{gc_content(sequence):.2f}%")
metrics[2].metric("Ambiguous bases", counts["N"])
metrics[3].metric("Selected frame", frame)

overview, transforms, coding, enzymes, orfs = st.tabs(["Overview", "Transforms", "Translation & codons", "Restriction sites", "ORFs"])

with overview:
    st.subheader("Nucleotide composition")
    st.bar_chart({"Count": {base: counts[base] for base in "ACGTN"}})
    st.code(sequence, language=None)

with transforms:
    st.subheader("Reverse complement")
    st.code(reverse_complement(sequence), language=None)
    st.subheader("RNA transcript")
    st.code(transcribe(sequence), language=None)

with coding:
    st.subheader(f"Translation — frame {frame}")
    st.code(translate(sequence, frame - 1), language=None)
    st.subheader("Codon frequencies")
    codons = codon_counts(sequence, frame - 1)
    rows = [{"Codon": codon, "Count": count, "Amino acid": translate(codon)} for codon, count in sorted(codons.items())]
    st.dataframe(rows, use_container_width=True, hide_index=True)

with enzymes:
    st.subheader("Restriction enzyme recognition sites")
    st.caption("Scans the forward strand for EcoRI, BamHI, HindIII, and HaeIII recognition sites.")
    sites = restriction_sites(sequence)
    if sites:
        st.dataframe(sites, use_container_width=True, hide_index=True)
    else:
        st.info("No recognition sites for the selected enzymes were found.")

with orfs:
    st.subheader("Forward-strand open reading frames")
    st.caption("Searches each forward reading frame for ATG start codons followed by an in-frame stop codon.")
    results = find_orfs(sequence, minimum_orf)
    if results:
        st.dataframe(results, use_container_width=True, hide_index=True)
    else:
        st.info("No ORFs met the selected minimum length.")

report = f"""BioSeq Explorer report
======================
Sequence length: {len(sequence)} bp
GC content: {gc_content(sequence):.2f}%
Estimated melting temperature (Wallace rule): {melting_temperature(sequence):.1f} °C
Selected reading frame: {frame}
Ambiguous bases (N): {counts['N']}

DNA sequence
------------
{sequence}

Reverse complement
------------------
{reverse_complement(sequence)}

RNA transcript
--------------
{transcribe(sequence)}

Protein translation (frame {frame})
------------------------------
{translate(sequence, frame - 1)}
"""
st.download_button("Download analysis report (.txt)", report, file_name="bioseq_report.txt", mime="text/plain")

st.divider()
st.caption("For learning and exploratory use only; melting temperature is a simple Wallace-rule estimate and this is not a clinical annotation tool.")
