# BioSeq Explorer

A beginner-friendly web application for exploring DNA sequences. Paste a sequence, select a reading frame, and get useful biological summaries in seconds.

![Python](https://img.shields.io/badge/Python-3.10%2B-blue) ![Streamlit](https://img.shields.io/badge/UI-Streamlit-red)

## Features

- Validates DNA and reports ambiguous bases
- Calculates sequence length, GC content, and nucleotide composition
- Generates the reverse complement and RNA transcript
- Translates all three forward reading frames
- Lists codon frequencies for a selected reading frame
- Finds open reading frames (ORFs) beginning with `ATG`
- Estimates melting temperature with the Wallace rule
- Detects recognition sites for EcoRI, BamHI, HindIII, and HaeIII
- Downloads a plain-text analysis report

## Run locally

```bash
git clone https://github.com/JAYsys23/-bioseq-explorer.git
cd bioseq-explorer
python -m venv .venv
# Windows
.venv\\Scripts\\activate
# macOS/Linux
source .venv/bin/activate
pip install -r requirements.txt
streamlit run app.py
```

## Deploy free with Streamlit Community Cloud

1. Create a new GitHub repository named `bioseq-explorer` under `JAYsys23` and upload these files.
2. Sign in to [Streamlit Community Cloud](https://share.streamlit.io/) using GitHub.
3. Choose the repository, branch, and `app.py` as the entry point.
4. Click **Deploy** and add the resulting link to your GitHub profile README.

## Biological notes

This is an educational sequence-exploration tool, not a clinical or research-grade annotation pipeline. ORFs are reported on the forward strand only and use the standard genetic code.

## Suggested GitHub portfolio description

> Interactive DNA sequence analysis web app built by Jayachandhran with Python and Streamlit. Includes GC-content calculation, reverse complement, transcription, translation, codon analysis, and ORF discovery.

## License

MIT
