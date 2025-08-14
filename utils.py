import re
from io import StringIO
from collections import Counter
from statistics import mean, pstdev
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

IUPAC_VALID = set(list("ACGTNRYWSKMBDHV"))

def read_fasta(file):
    if hasattr(file, "getvalue"):
        text_stream = StringIO(file.getvalue().decode("utf-8"))
    else:
        text_stream = open(file, "r")
    record = next(SeqIO.parse(text_stream, "fasta"))
    return record.id, str(record.seq).upper()

def gc_content(seq):
    return (seq.count("G") + seq.count("C")) / len(seq) * 100

def codon_frequency(seq):
    counts = Counter()
    for i in range(0, len(seq)-2, 3):
        codon = seq[i:i+3]
        if re.fullmatch(r"[ACGT]{3}", codon):
            counts[codon] += 1
    total = sum(counts.values())
    return {cod: counts[cod] / total for cod in counts}

def motif_search(seq, motif="ATG"):
    return [m.start()+1 for m in re.finditer(f"(?={motif})", seq)]

def compare_sequences(seq1, seq2, max_len=2000):
    seq1 = seq1[:max_len]
    seq2 = seq2[:max_len]
    aligner = PairwiseAligner()
    aligner.mode = 'global'
    score = aligner.score(seq1, seq2)
    similarity = (score / max(len(seq1), len(seq2))) * 100
    return similarity

def find_ambiguous_bases(seq):
    bad_positions = [i+1 for i, b in enumerate(seq) if b not in "ACGT"]
    n_runs = []
    for m in re.finditer(r"N{5,}", seq):
        n_runs.append({"start": m.start()+1, "end": m.end(), "length": m.end()-m.start()})
    counts = {c: seq.count(c) for c in set(seq) if c not in "ACGT"}
    return {
        "ambiguous_positions": bad_positions[:1000],
        "ambiguous_total": len(bad_positions),
        "iupac_counts": counts,
        "n_runs": n_runs
    }

def sliding_gc(seq, win=500):
    vals = []
    for i in range(0, len(seq), win):
        chunk = seq[i:i+win]
        if not chunk:
            break
        gc = (chunk.count("G")+chunk.count("C"))/len(chunk)*100
        vals.append({"start": i+1, "end": i+len(chunk), "gc": gc})
    return vals

def gc_outliers(gc_windows, z=2.5):
    arr = [w["gc"] for w in gc_windows]
    if len(arr) < 3:
        return []
    mu, sd = mean(arr), pstdev(arr) or 1e-9
    out = []
    for w in gc_windows:
        zsc = (w["gc"]-mu)/sd
        if abs(zsc) >= z:
            w2 = w.copy()
            w2["z"] = zsc
            out.append(w2)
    return out

def premature_stop_flags(seq, min_orf=150):
    flags = []
    for frame in [0,1,2]:
        prot = Seq(seq[frame:]).translate(to_stop=False)
        stops = [i for i,aa in enumerate(prot) if aa=="*"]
        if stops and stops[0]*3 < min_orf:
            flags.append({"frame": frame, "first_stop_nt": stops[0]*3+frame+1})
    return flags

def simple_snp_diff(seq, ref, max_len=5000):
    L = min(len(seq), len(ref), max_len)
    snps = []
    for i in range(L):
        a, b = seq[i], ref[i]
        if a in IUPAC_VALID and b in IUPAC_VALID and a!=b and a in "ACGT" and b in "ACGT":
            snps.append({"pos": i+1, "ref": b, "alt": a})
    return {
        "checked_bases": L,
        "snp_count": len(snps),
        "snps_preview": snps[:1000]
    }
