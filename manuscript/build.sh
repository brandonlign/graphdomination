#!/bin/sh
# Build the article, then the two supplementary documents.
set -eu
cd "$(dirname "$0")"
engine=${LATEX_ENGINE:-pdflatex}
if ! command -v "$engine" >/dev/null 2>&1; then
  if [ -z "${LATEX_ENGINE:-}" ] && command -v tectonic >/dev/null 2>&1; then
    engine=tectonic
  else
    echo "Required LaTeX engine not found: $engine" >&2
    exit 1
  fi
fi
export SOURCE_DATE_EPOCH=1793448000 FORCE_SOURCE_DATE=1
build() {
  if [ "$engine" = tectonic ]; then
    tectonic --keep-logs --keep-intermediates "$1.tex"
  else
    for pass in 1 2 3; do
      "$engine" -interaction=nonstopmode -halt-on-error -file-line-error "$1.tex" > "$1.build.log" 2>&1 || {
        cat "$1.build.log" >&2
        exit 1
      }
    done
  fi
}
build manuscript
python3 export_article_references.py
build subdivision_supplement
build verification_guide
echo "Built all three PDFs with $engine."
