#!/usr/bin/env python3
"""Write the article label numbers used by the verification guide."""
from pathlib import Path
import re

root = Path(__file__).resolve().parent
aux = (root / 'manuscript.aux').read_text()
guide = (root / 'verification_guide.tex').read_text()
keys = sorted(set(re.findall(r'\\art(?:eq)?ref\{([^}]+)\}', guide)))
labels = dict(re.findall(r'\\newlabel\{([^}]+)\}\{\{([^{}]*)\}', aux))
missing = set(keys) - labels.keys()
if missing:
    raise SystemExit(f'Unresolved article labels: {sorted(missing)}')
lines = [
    '% Generated from manuscript.aux by export_article_references.py.',
        r'\newcommand{\artref}[1]{\ifcsname artref@#1\endcsname\csname artref@#1\endcsname\else\PackageError{article-references}{Unknown article label #1}{}\fi}',
    r'\newcommand{\arteqref}[1]{\textup{(\artref{#1})}}',
]
for key in keys:
    lines.append(r'\expandafter\def\csname artref@' + key + r'\endcsname{' + labels[key] + '}')
(root / 'article_references.tex').write_text('\n'.join(lines) + '\n')
print(f'Exported {len(keys)} article references for the guide.')
