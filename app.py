import streamlit as st
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from fpdf import FPDF
import os
from utils import (
    read_fasta, gc_content, codon_frequency, motif_search, compare_sequences,
    find_ambiguous_bases, sliding_gc, gc_outliers, premature_stop_flags, simple_snp_diff
)

st.set_page_config(page_title="DNA Analysis & Comparison Tool", layout="wide")
st.markdown("<h1 style='text-align:center; color:#2E86C1;'>🧬 DNA Analysis & Comparison Tool</h1>", unsafe_allow_html=True)
st.markdown("---")

col_u1, col_u2 = st.columns(2)
with col_u1:
    file1 = st.file_uploader("📂 Upload DNA Sequence File 1 (FASTA)", type=["fasta", "fa"])
with col_u2:
    file2 = st.file_uploader("📂 Upload DNA Sequence File 2 (FASTA)", type=["fasta", "fa"])

motifs = ["ATG", "TATA", "GGC", "TTGACA", "Custom"]
motif_choice = st.selectbox("🔍 Select motif to search", motifs)
if motif_choice == "Custom":
    motif_input = st.text_input("✏️ Enter custom motif:", value="ATG")
else:
    motif_input = motif_choice

download_report = st.checkbox("📄 Generate Detailed PDF Report")

if file1 and file2:
    id1, seq1 = read_fasta(file1)
    id2, seq2 = read_fasta(file2)

    gc1, gc2 = gc_content(seq1), gc_content(seq2)
    codon1, codon2 = codon_frequency(seq1), codon_frequency(seq2)
    motif_pos1, motif_pos2 = motif_search(seq1, motif_input), motif_search(seq2, motif_input)
    similarity = compare_sequences(seq1[:2000], seq2[:2000])  # fast subset

    st.markdown("## 📊 Sequence Analysis")
    col1, col2 = st.columns(2)
    with col1:
        st.success(f"**{id1}**\n\nLength: {len(seq1)}\nGC%: {gc1:.2f}\nMotif '{motif_input}' Count: {len(motif_pos1)}")
    with col2:
        st.info(f"**{id2}**\n\nLength: {len(seq2)}\nGC%: {gc2:.2f}\nMotif '{motif_input}' Count: {len(motif_pos2)}")

    fig1, ax = plt.subplots(1, 2, figsize=(12, 5))
    codon_df1 = pd.DataFrame(codon1.items(), columns=["Codon", "Frequency"]).sort_values(by="Frequency", ascending=False)[:10]
    codon_df2 = pd.DataFrame(codon2.items(), columns=["Codon", "Frequency"]).sort_values(by="Frequency", ascending=False)[:10]
    ax[0].bar(codon_df1["Codon"], codon_df1["Frequency"], color="#2E86C1")
    ax[0].set_title(f"Top Codons - {id1}")
    ax[1].bar(codon_df2["Codon"], codon_df2["Frequency"], color="#28B463")
    ax[1].set_title(f"Top Codons - {id2}")
    st.pyplot(fig1)

    fig2, ax2 = plt.subplots(2, 1, figsize=(10, 5))
    ax2[0].stem(motif_pos1, [1]*len(motif_pos1), linefmt="C0-", markerfmt="C0o")
    ax2[0].set_title(f"Motif Positions - {id1}")
    ax2[1].stem(motif_pos2, [1]*len(motif_pos2), linefmt="C1-", markerfmt="C1o")
    ax2[1].set_title(f"Motif Positions - {id2}")
    st.pyplot(fig2)

    comparison_data = pd.DataFrame({
        "GC%": [gc1, gc2],
        "Motif Count": [len(motif_pos1), len(motif_pos2)],
        "Similarity %": [similarity, similarity]
    }, index=[id1, id2])

    fig3, ax3 = plt.subplots(figsize=(6, 4))
    sns.heatmap(comparison_data, annot=True, cmap="coolwarm", fmt=".2f", ax=ax3, cbar_kws={'label': 'Value'})
    ax3.set_title("DNA Comparison Heatmap")
    st.pyplot(fig3)

    st.markdown("## ⚖️ Comparison Summary")
    st.write(f"**GC% Difference:** {abs(gc1 - gc2):.2f}")
    st.write(f"**Motif Count Difference:** {abs(len(motif_pos1) - len(motif_pos2))}")
    st.write(f"**Sequence Similarity (first 2000 bases):** {similarity:.2f}%")

    st.markdown("## 🧪 Defect Scan")
    run_defect = st.checkbox("Run defect/quality scan")
    ref_file = st.file_uploader("Optional: Upload Reference FASTA for quick SNP check", type=["fasta","fa"])

    if run_defect:
        amb1 = find_ambiguous_bases(seq1)
        amb2 = find_ambiguous_bases(seq2)

        st.write("**Unknown/IUPAC bases**")
        colA, colB = st.columns(2)
        with colA:
            st.write(f"{id1}: total ambiguous = {amb1['ambiguous_total']}")
            st.write(f"IUPAC counts: {amb1['iupac_counts']}")
        with colB:
            st.write(f"{id2}: total ambiguous = {amb2['ambiguous_total']}")
            st.write(f"IUPAC counts: {amb2['iupac_counts']}")

        st.write("**GC window outliers (z≥2.5)**")
        gcw1 = sliding_gc(seq1, win=1000)
        gcw2 = sliding_gc(seq2, win=1000)
        out1 = gc_outliers(gcw1)
        out2 = gc_outliers(gcw2)
        st.write(f"{id1}: {len(out1)} outlier windows")
        st.write(f"{id2}: {len(out2)} outlier windows")

        st.write("**Premature stop flags (rough ORF check)**")
        st.write(f"{id1}: {premature_stop_flags(seq1)}")
        st.write(f"{id2}: {premature_stop_flags(seq2)}")

        if ref_file is not None:
            ref_id, ref_seq = read_fasta(ref_file)
            st.write(f"**Quick SNP diff vs reference ({ref_id})** (first 5000 bases)")
            diff1 = simple_snp_diff(seq1, ref_seq)
            st.write(f"{id1}: SNPs={diff1['snp_count']} (checked {diff1['checked_bases']} bases)")
            diff2 = simple_snp_diff(seq2, ref_seq)
            st.write(f"{id2}: SNPs={diff2['snp_count']} (checked {diff2['checked_bases']} bases)")

    # ----- PDF Report -----
    if download_report:
        report_path = "dna_report.pdf"
        codon_plot_img = "codon_plot.png"
        motif_plot_img = "motif_plot.png"
        heatmap_img = "heatmap.png"
        fig1.savefig(codon_plot_img)
        fig2.savefig(motif_plot_img)
        fig3.savefig(heatmap_img)

        pdf = FPDF()
        pdf.add_page()

        pdf.set_font("Arial", 'B', 16)
        pdf.cell(200, 10, "DNA Analysis Report", ln=True, align="C")
        pdf.ln(5)

        # Sequence details
        for name, length, gc_val, motif_count in [
            (id1, len(seq1), gc1, len(motif_pos1)),
            (id2, len(seq2), gc2, len(motif_pos2))
        ]:
            pdf.set_font("Arial", 'B', 12)
            pdf.cell(200, 8, f"Sequence: {name}", ln=True)
            pdf.set_font("Arial", '', 11)
            pdf.multi_cell(0, 6, f"Length: {length}\nGC%: {gc_val:.2f}\nMotif '{motif_input}' Count: {motif_count}")

        # Comparison
        pdf.set_font("Arial", 'B', 12)
        pdf.cell(200, 8, "Comparison", ln=True)
        pdf.set_font("Arial", '', 11)
        pdf.multi_cell(0, 6, f"GC% Difference: {abs(gc1 - gc2):.2f}\n"
                             f"Motif Count Difference: {abs(len(motif_pos1) - len(motif_pos2))}\n"
                             f"Similarity (first 2000 bases): {similarity:.2f}%")

        # Graphs
        pdf.ln(5)
        pdf.set_font("Arial", 'B', 12)
        pdf.cell(200, 8, "Codon Usage Comparison", ln=True)
        pdf.image(codon_plot_img, x=10, w=180)

        pdf.ln(5)
        pdf.cell(200, 8, f"Motif Positions for '{motif_input}'", ln=True)
        pdf.image(motif_plot_img, x=10, w=180)

        pdf.ln(5)
        pdf.cell(200, 8, "Comparison Heatmap", ln=True)
        pdf.image(heatmap_img, x=10, w=180)

        pdf.output(report_path)

        os.remove(codon_plot_img)
        os.remove(motif_plot_img)
        os.remove(heatmap_img)

        with open(report_path, "rb") as f:
            st.download_button("⬇️ Download Detailed Report", f, file_name="dna_report.pdf", mime="application/pdf")
