"""Self-contained HTML report generation."""
from __future__ import annotations

import base64
from pathlib import Path


def _b64_image(path: str) -> str:
    ext = Path(path).suffix.lower().lstrip(".")
    mime = "image/svg+xml" if ext == "svg" else ("image/png" if ext == "png" else "image/jpeg")
    data = base64.b64encode(Path(path).read_bytes()).decode("ascii")
    return f"data:{mime};base64,{data}"


def render_html(title: str, sections: dict, output: str) -> str:
    """Render a simple HTML report.

    sections: dict {section_title: html_body}.
    """
    parts = [f"<html><head><meta charset='utf-8'><title>{title}</title>"
             "<style>body{font-family:Segoe UI,Arial,sans-serif;margin:2rem;color:#222}"
             "h1{border-bottom:3px solid #2c6fbb;padding-bottom:.5rem}"
             "h2{color:#2c6fbb;margin-top:2rem}"
             "table{border-collapse:collapse;width:100%;margin:1rem 0}"
             "th,td{border:1px solid #ccc;padding:.4rem .6rem;text-align:left;font-size:.9rem}"
             "img{max-width:100%;border:1px solid #eee;margin:.5rem 0}"
             ".note{color:#666;font-size:.85rem}</style></head><body>"]
    parts.append(f"<h1>{title}</h1>")
    for sec_title, body in sections.items():
        parts.append(f"<h2>{sec_title}</h2>{body}")
    parts.append("</body></html>")
    html = "".join(parts)
    Path(output).write_text(html, encoding="utf-8")
    return output


def figure_block(path: str, caption: str = "") -> str:
    img = _b64_image(path)
    cap = f"<div class='note'>{caption}</div>" if caption else ""
    return f"<div><img src='{img}'>{cap}</div>"


def table_block(df, max_rows: int = 50) -> str:
    return df.head(max_rows).to_html(index=False, border=0, classes="report-table")