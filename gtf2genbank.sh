#!/usr/bin/env bash
# SnapGene GenBank maker: GTF ⟶ (FASTA + minimal GFF3) ⟶ .gbk
# - Default: transcript + exon
# - Add UTR with -u/--utr
# - Input selector: exactly one of -g SYMBOL | --ensg ENSG | --enst ENST
# - ENSG/ENST: version suffix (.N) is ignored for matching and for filenames

set -Eeuo pipefail
trap 'echo "[ERR] line $LINENO: $BASH_COMMAND" >&2' ERR

usage() {
  cat <<'USAGE'
make_snapgene_views.sh — Subset GTF by gene symbol / ENSG / ENST, build locus FASTA, create minimal GFF3, and write GenBank (.gbk) for SnapGene.

Required (choose exactly one):
  -g,  --gene    SYMBOL     Gene symbol (e.g., TP53)
       --ensg    ENSG_ID    Ensembl gene ID (e.g., ENSG00000141510[.N])
       --enst    ENST_ID    Ensembl transcript ID (e.g., ENST00000413465[.N])

Inputs:
  -a,  --gtf     FILE       GTF file (GENCODE/Ensembl)
  -f,  --fasta   FILE       Genome FASTA (indexed with .fai)

Options:
  -o,  --outdir  DIR        Output dir (default: ./snapgene_<LABEL>)
  -u,  --utr                Include UTR features (default off)
  -h,  --help               Show this help

Outputs:
  <OUTDIR>/<LABEL>.locus.fa
  <OUTDIR>/<LABEL>.snapgene.gff3
  <OUTDIR>/<LABEL>.gbk
USAGE
}

# -------- args --------
GENE=""; ENSG=""; ENST=""
GTF=""; GENOME=""; OUTDIR=""
INCLUDE_UTR=0

[[ $# -gt 0 ]] || { usage; exit 1; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    -g|--gene)   GENE="${2:-}"; shift 2 ;;
    --ensg)      ENSG="${2:-}"; shift 2 ;;
    --enst)      ENST="${2:-}"; shift 2 ;;
    -a|--gtf)    GTF="${2:-}"; shift 2 ;;
    -f|--fasta)  GENOME="${2:-}"; shift 2 ;;
    -o|--outdir) OUTDIR="${2:-}"; shift 2 ;;
    -u|--utr)    INCLUDE_UTR=1; shift 1 ;;
    -h|--help)   usage; exit 0 ;;
    *) echo "[ERR] Unknown option: $1" >&2; usage; exit 1 ;;
  esac
done

# -------- validate selector / inputs --------
sel_count=0
[[ -n "$GENE" ]] && ((++sel_count))
[[ -n "$ENSG" ]] && ((++sel_count))
[[ -n "$ENST" ]] && ((++sel_count))
if [[ $sel_count -ne 1 ]]; then
  echo "[ERR] specify exactly one of: -g SYMBOL | --ensg ENSG | --enst ENST" >&2
  exit 1
fi

[[ -n "$GTF"    ]] || { echo "[ERR] --gtf is required" >&2; exit 1; }
[[ -n "$GENOME" ]] || { echo "[ERR] --fasta is required" >&2; exit 1; }

# -------- deps --------
command -v bedtools >/dev/null 2>&1 || { echo "[ERR] bedtools not found"; exit 1; }
python3 -c "import Bio" >/dev/null 2>&1 || { echo "[ERR] biopython not found (pip/mamba install biopython)" >&2; exit 1; }
[[ -f "$GTF"    ]] || { echo "[ERR] GTF not found: $GTF"; exit 1; }
[[ -f "$GENOME" ]] || { echo "[ERR] FASTA not found: $GENOME"; exit 1; }

# -------- subset from GTF (version-insensitive for ENSG/ENST) --------
TMPDIR="$(mktemp -d)"; trap 'rm -rf "$TMPDIR"' EXIT
GENE_GTF="${TMPDIR}/subset.gtf"

if [[ -n "$GENE" ]]; then
  # SYMBOL: fixed-string match
  awk -F'\t' -v g="$GENE"  'index($9, "gene_name \"" g "\"")>0' "$GTF" > "$GENE_GTF"
  LABEL="$GENE"
elif [[ -n "$ENSG" ]]; then
  # ENSG: strip version, match core id
  core="${ENSG%%.*}"
  awk -F'\t' -v g="$core" '
    match($9, /gene_id "([^"]+)"/, m) { split(m[1],a,"."); if (a[1]==g) print }
  ' "$GTF" > "$GENE_GTF"
  LABEL="$core"
else
  # ENST: strip version, match core id
  core="${ENST%%.*}"
  awk -F'\t' -v g="$core" '
    match($9, /transcript_id "([^"]+)"/, m) { split(m[1],a,"."); if (a[1]==g) print }
  ' "$GTF" > "$GENE_GTF"
  LABEL="$core"
fi

[[ -s "$GENE_GTF" ]] || { echo "[ERR] target not found in GTF (LABEL=${LABEL:-unknown})" >&2; exit 1; }

# -------- outdir / filenames --------
[[ -n "$OUTDIR" ]] || OUTDIR="./snapgene_${LABEL}"
mkdir -p "$OUTDIR"

GFF3="${OUTDIR}/${LABEL}.snapgene.gff3"
LOCUS_BED="${OUTDIR}/${LABEL}.locus.bed"
LOCUS_FA="${OUTDIR}/${LABEL}.locus.fa"
GBK="${OUTDIR}/${LABEL}.gbk"

# -------- locus: min exon start ~ max exon end --------
read -r CHR START END STRAND < <(
  awk -F'\t' '
    BEGIN{min=10^12; max=0; chr=""; strand="+"}
    $3=="exon"{
      if (chr=="") chr=$1
      if ($4<min) min=$4
      if ($5>max) max=$5
      if ($7=="+"||$7=="-") strand=$7
    }
    END{print chr, min, max, strand}
  ' "$GENE_GTF"
)
[[ -n "${CHR}" && "${START}" -gt 0 && "${END}" -ge "${START}" ]] || { echo "[ERR] exon locus not found (LABEL=$LABEL)" >&2; exit 1; }

# +strand 固定で getfasta（相対座標計算を単純化）
echo -e "${CHR}\t$((START-1))\t${END}\t${LABEL}\t.\t+" > "$LOCUS_BED"
bedtools getfasta -fi "$GENOME" -bed "$LOCUS_BED" -fo "$LOCUS_FA"

# -------- GFF3: default transcript+exon, add UTR if -u --------
awk -F'\t' -v OFS='\t' -v addutr="$INCLUDE_UTR" '
  BEGIN{ print "##gff-version 3" }
  ($3=="transcript" || $3=="mRNA" || $3=="exon" || (addutr==1 && $3=="UTR")) {
    if (match($9, /transcript_id "([^"]+)"/, m)) {
      tid=m[1]
      if (!(tid in chr))    chr[tid]=$1
      if (!(tid in strand)) strand[tid]=$7

      # mRNA span 更新
      if ($3=="exon") {
        if (!(tid in min)) { min[tid]=$4; max[tid]=$5 }
        else { if ($4<min[tid]) min[tid]=$4; if ($5>max[tid]) max[tid]=$5 }
      } else if (!(tid in min)) { min[tid]=$4; max[tid]=$5 }

      # exon
      if ($3=="exon") {
        exn=""
        if      (match($9, /exon_number[[:space:]]*"([^"]+)"/, e)) exn=e[1];
        else if (match($9, /exon_number[[:space:]]*([^;]+)/,   e)) exn=e[1];
        if (exn=="") { ex_i[tid]++; exn=ex_i[tid] }
        key=tid "|" $1 "|" $4 "|" $5
        if (!(key in ex_seen)) {
          ex_seen[key]=1
          name = tid "|exon" exn
          exon_lines[++nx] = sprintf("%s\tgencode\texon\t%d\t%d\t.\t%s\t.\tID=exon:%s:%s;Parent=transcript:%s;Name=%s",
                                     $1, $4, $5, $7, tid, exn, tid, name)
        }
      }

      # UTR（-u のときのみ）
      if (addutr==1 && $3=="UTR") {
        utr_lines[++nu] = sprintf("%s\tgencode\tUTR\t%d\t%d\t.\t%s\t.\tParent=transcript:%s;Name=UTR",
                                  $1, $4, $5, $7, tid)
      }
    }
  }
  END{
    for (t in min) {
      printf("%s\tgencode\tmRNA\t%d\t%d\t.\t%s\t.\tID=transcript:%s;Name=%s\n",
             chr[t], min[t], max[t], strand[t], t, t)
    }
    for (i=1;i<=nx;i++) print exon_lines[i]
    for (j=1;j<=nu;j++) print utr_lines[j]
  }
' "$GENE_GTF" > "$GFF3"

# -------- FASTA + GFF3 -> GenBank (.gbk) --------
python3 - "$LOCUS_FA" "$GFF3" "$GBK" "$LABEL" "$((START-1))" "$END" <<'PY'
import sys, re
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation

fa, gff, out_gbk, label, start0_str, end1_str = sys.argv[1:]
start0, end1 = int(start0_str), int(end1_str)

rec = next(SeqIO.parse(fa, "fasta"))
rec.id = label
rec.name = label[:10]
rec.description = f"{label} locus {start0+1}-{end1}"
rec.annotations["molecule_type"] = "DNA"

tid_re  = re.compile(r'(?:transcript_id=([^;]+))|(?:transcript_id "([^"]+)")')
name_re = re.compile(r'(?:Name=([^;]+))')

features = []
with open(gff) as fh:
  for line in fh:
    if line.startswith("#") or not line.strip():
      continue
    cols = line.rstrip("\n").split("\t")
    if len(cols) != 9:
      continue
    chrom, src, ftype, beg, end, score, strand, phase, attrs = cols
    beg_i, end_i = int(beg), int(end)     # GFF3: 1-based inclusive
    beg0 = beg_i - 1 - start0             # locus.fa は +鎖固定 → 0-based
    end0 = end_i - start0
    if beg0 < 0 or end0 > len(rec.seq):
      continue
    strand_val = +1 if strand == "+" else -1

    m = tid_re.search(attrs)
    tid = (m.group(1) or m.group(2)) if m else None
    nm = name_re.search(attrs)
    label_attr = (nm.group(1) if nm else (tid if tid else ftype))

    # feature key 正規化
    t = ftype.lower()
    if t in ("mrna","transcript"):
        ftype = "mRNA"
    elif t in ("utr",):
        ftype = "UTR"

    loc = FeatureLocation(beg0, end0, strand=strand_val)
    quals = {"label":[label_attr], "note":[label_attr]}
    if tid:
        quals["transcript_id"] = [tid]
    features.append(SeqFeature(loc, type=ftype, qualifiers=quals))

# 重複除去（座標×タイプ×ラベル）
seen=set(); uniq=[]
for f in features:
    key=(int(f.location.start), int(f.location.end), f.type, f.qualifiers.get("label",[""])[0])
    if key in seen:
        continue
    seen.add(key); uniq.append(f)
rec.features = uniq

SeqIO.write(rec, out_gbk, "genbank")
PY

echo "=== OUTPUT ==="
echo "FASTA : $LOCUS_FA"
echo "GFF3  : $GFF3"
echo "GBK   : $GBK"
echo "Open ${GBK} in SnapGene Viewer."
