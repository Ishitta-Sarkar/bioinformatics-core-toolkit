"""Write the HBB base-composition figure from the bundled FASTA file.

The counts come from ``DNATools``. The figure is an SVG so the repository
does not need a plotting library.
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FASTA = ROOT / "examples" / "sample_sequence.fasta"
OUT = ROOT / "docs" / "hbb_base_composition.svg"

COLORS = {
    "A": "#0072B2",
    "C": "#E69F00",
    "G": "#009E73",
    "T": "#D55E00",
}


def render_svg(sequence: str) -> str:
    """Return an SVG bar chart of A, C, G, and T counts."""
    from src.dna_tools import DNATools

    dna = DNATools(sequence)
    if not dna.is_valid_dna() or dna.length() == 0:
        raise ValueError("sequence must contain at least one A, C, G, or T")

    counts = {base: dna.sequence.count(base) for base in "ACGT"}
    length = dna.length()
    gc = dna.gc_content()
    percents = {base: round(counts[base] / length * 100, 2) for base in "ACGT"}

    width = 760
    height = 460
    left = 72
    right = 32
    top = 108
    bottom = 86
    plot_w = width - left - right
    plot_h = height - top - bottom
    y_max = 160
    baseline = top + plot_h
    bar_w = 78
    gap = 42
    total = 4 * bar_w + 3 * gap
    origin = left + (plot_w - total) / 2

    def y_of(value: float) -> float:
        return top + plot_h * (1 - value / y_max)

    parts = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
            f'height="{height}" viewBox="0 0 {width} {height}" '
            'role="img">'
        ),
        (
            "<title>Base composition of the human HBB coding sequence "
            "NM_000518.5:51-494</title>"
        ),
        (
            f"<desc>A {counts['A']} ({percents['A']:.2f}%), "
            f"C {counts['C']} ({percents['C']:.2f}%), "
            f"G {counts['G']} ({percents['G']:.2f}%), "
            f"T {counts['T']} ({percents['T']:.2f}%). "
            f"GC content {gc:.2f}% of {length} bp.</desc>"
        ),
        '<rect width="100%" height="100%" fill="#ffffff"/>',
        (
            '<text x="380" y="36" text-anchor="middle" '
            'font-family="Helvetica, Arial, sans-serif" font-size="20" '
            'fill="#1a1a1a">Base composition of the human HBB coding '
            "sequence</text>"
        ),
        (
            '<text x="380" y="62" text-anchor="middle" '
            'font-family="Helvetica, Arial, sans-serif" font-size="14" '
            f'fill="#333333">NM_000518.5:51-494 · {length} bp · '
            f"GC content {gc:.2f}%</text>"
        ),
    ]

    for tick in (0, 40, 80, 120, 160):
        y = y_of(tick)
        parts.append(
            f'<line x1="{left}" y1="{y:.1f}" x2="{width - right}" '
            f'y2="{y:.1f}" stroke="#e6e6e6" stroke-width="1"/>'
        )
        parts.append(
            f'<text x="{left - 10}" y="{y + 4:.1f}" text-anchor="end" '
            'font-family="Helvetica, Arial, sans-serif" font-size="12" '
            f'fill="#333333">{tick}</text>'
        )

    parts.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{baseline}" '
        'stroke="#1a1a1a" stroke-width="1.5"/>'
    )
    parts.append(
        f'<line x1="{left}" y1="{baseline}" x2="{width - right}" '
        f'y2="{baseline}" stroke="#1a1a1a" stroke-width="1.5"/>'
    )
    parts.append(
        f'<text x="22" y="{top + plot_h / 2:.1f}" '
        'transform="rotate(-90 22 '
        f'{top + plot_h / 2:.1f})" text-anchor="middle" '
        'font-family="Helvetica, Arial, sans-serif" font-size="13" '
        'fill="#1a1a1a">Count</text>'
    )

    for index, base in enumerate("ACGT"):
        x = origin + index * (bar_w + gap)
        y = y_of(counts[base])
        bar_h = baseline - y
        parts.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w}" '
            f'height="{bar_h:.1f}" fill="{COLORS[base]}" rx="3">'
            f"<title>{base}: {counts[base]} ({percents[base]:.2f}%)</title>"
            "</rect>"
        )
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{y - 22:.1f}" '
            'text-anchor="middle" '
            'font-family="Helvetica, Arial, sans-serif" font-size="13" '
            f'fill="#1a1a1a">{counts[base]}</text>'
        )
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{y - 6:.1f}" '
            'text-anchor="middle" '
            'font-family="Helvetica, Arial, sans-serif" font-size="12" '
            f'fill="#333333">{percents[base]:.2f}%</text>'
        )
        parts.append(
            f'<text x="{x + bar_w / 2:.1f}" y="{baseline + 24:.1f}" '
            'text-anchor="middle" '
            'font-family="Helvetica, Arial, sans-serif" font-size="16" '
            f'font-weight="bold" fill="#1a1a1a">{base}</text>'
        )

    parts.append(
        f'<text x="380" y="{height - 22}" text-anchor="middle" '
        'font-family="Helvetica, Arial, sans-serif" font-size="12" '
        'fill="#555555">Counts from DNATools on '
        "examples/sample_sequence.fasta</text>"
    )
    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def main() -> None:
    if str(ROOT) not in sys.path:
        sys.path.insert(0, str(ROOT))
    from src.dna_tools import DNATools

    records = DNATools.parse_fasta(str(FASTA))
    _header, sequence = next(iter(records.items()))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(render_svg(sequence), encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
