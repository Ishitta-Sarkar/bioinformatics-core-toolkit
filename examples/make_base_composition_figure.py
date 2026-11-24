"""Write the base-composition figure for the bundled HBB coding sequence.

The counts come from ``DNATools`` applied to ``examples/sample_sequence.fasta``.
Re-run this script from the repository root or from any working directory
after changing that file:

    python examples/make_base_composition_figure.py
"""

from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
FASTA_PATH = ROOT / "examples" / "sample_sequence.fasta"
SVG_PATH = ROOT / "examples" / "hbb_cds_base_composition.svg"

_COLORS = {
    "A": "#0072B2",
    "C": "#E69F00",
    "G": "#009E73",
    "T": "#D55E00",
}


def render_hbb_base_composition_svg(counts: dict[str, int], gc_percent: float) -> str:
    """Return an SVG bar chart of A/C/G/T counts.

    ``counts`` must contain A, C, G, and T. ``gc_percent`` is the value
    returned by ``DNATools.gc_content`` and is printed to two decimal places.
    Percentages on the bars are ``count / total * 100``, rounded to two
    decimal places. The title cites NM_000518.3 nucleotides 51-494 because
    that is the sequence this figure is built for.
    """
    order = ("A", "C", "G", "T")
    if set(counts) != set(order):
        raise ValueError("counts must contain exactly the keys A, C, G, and T")
    total = sum(counts[base] for base in order)
    if total <= 0:
        raise ValueError("total base count must be positive")
    if max(counts.values()) > 160:
        raise ValueError("a count above 160 does not fit the fixed axis")

    width = 720
    height = 460
    left = 76
    right = 32
    top = 84
    bottom = 110
    plot_w = width - left - right
    plot_h = height - top - bottom
    y_max = 160
    axis_bottom = top + plot_h
    gc_text = f"{gc_percent:.2f}"

    def y_of(count: int) -> int:
        return axis_bottom - round(count / y_max * plot_h)

    lines = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        "<!-- bases "
        + " ".join(f"{base}={counts[base]}" for base in order)
        + f" gc={gc_text} n={total} -->",
        (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" '
            f'height="{height}" viewBox="0 0 {width} {height}" role="img">'
        ),
        (
            "<title>Base composition of the HBB coding sequence "
            "(NM_000518.3:51-494)</title>"
        ),
        f'<rect x="0" y="0" width="{width}" height="{height}" fill="#ffffff"/>',
        (
            '<text x="360" y="32" text-anchor="middle" '
            'font-family="Arial, Helvetica, sans-serif" font-size="18" fill="#1a1a1a">'
            "Base composition of the HBB coding sequence</text>"
        ),
        (
            '<text x="360" y="56" text-anchor="middle" '
            'font-family="Arial, Helvetica, sans-serif" font-size="13" fill="#444444">'
            f"NM_000518.3 nucleotides 51-494, {total} bp. GC content {gc_text}%.</text>"
        ),
    ]

    for tick in (0, 40, 80, 120, 160):
        y = y_of(tick)
        lines.append(
            f'<line x1="{left}" y1="{y}" x2="{left + plot_w}" y2="{y}" '
            'stroke="#e6e6e6" stroke-width="1"/>'
        )
        lines.append(
            f'<text x="{left - 8}" y="{y + 4}" text-anchor="end" '
            'font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#333333">'
            f"{tick}</text>"
        )

    lines.append(
        f'<line x1="{left}" y1="{top}" x2="{left}" y2="{axis_bottom}" '
        'stroke="#222222" stroke-width="1.5"/>'
    )
    lines.append(
        f'<line x1="{left}" y1="{axis_bottom}" x2="{left + plot_w}" y2="{axis_bottom}" '
        'stroke="#222222" stroke-width="1.5"/>'
    )
    label_y = top + plot_h // 2
    lines.append(
        f'<text x="22" y="{label_y}" transform="rotate(-90 22 {label_y})" '
        'text-anchor="middle" font-family="Arial, Helvetica, sans-serif" '
        'font-size="13" fill="#222222">Count</text>'
    )

    slot = plot_w / 4
    bar_w = 78
    for index, base in enumerate(order):
        count = counts[base]
        percent = count / total * 100
        x = round(left + index * slot + (slot - bar_w) / 2)
        y = y_of(count)
        bar_height = axis_bottom - y
        center = x + bar_w // 2
        lines.append(
            f'<rect x="{x}" y="{y}" width="{bar_w}" height="{bar_height}" '
            f'fill="{_COLORS[base]}"/>'
        )
        lines.append(
            f'<text x="{center}" y="{y - 8}" text-anchor="middle" '
            'font-family="Arial, Helvetica, sans-serif" font-size="13" fill="#1a1a1a">'
            f"{count}</text>"
        )
        lines.append(
            f'<text x="{center}" y="{axis_bottom + 24}" text-anchor="middle" '
            'font-family="Arial, Helvetica, sans-serif" font-size="16" '
            'font-weight="bold" fill="#1a1a1a">'
            f"{base}</text>"
        )
        lines.append(
            f'<text x="{center}" y="{axis_bottom + 44}" text-anchor="middle" '
            'font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#444444">'
            f"{percent:.2f}%</text>"
        )

    lines.append(
        f'<text x="{left}" y="{height - 28}" '
        'font-family="Arial, Helvetica, sans-serif" font-size="12" fill="#555555">'
        "Counts from examples/sample_sequence.fasta using DNATools. "
        "Percentages are rounded to two decimal places.</text>"
    )
    lines.append("</svg>")
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    """Read the example FASTA and overwrite the SVG figure."""
    sys.path.insert(0, str(ROOT))
    from src.dna_tools import DNATools

    records = DNATools.parse_fasta(FASTA_PATH)
    if len(records) != 1:
        raise SystemExit(f"Expected one FASTA record, found {len(records)}")
    sequence = DNATools(next(iter(records.values()))).sequence
    counts = {base: sequence.count(base) for base in "ACGT"}
    svg = render_hbb_base_composition_svg(counts, DNATools(sequence).gc_content())
    SVG_PATH.write_text(svg, encoding="utf-8")


if __name__ == "__main__":
    main()
