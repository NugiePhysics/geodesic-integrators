"""MkDocs hooks: GitHub-flavoured math, figures from the repository, links that leave docs/.

The Markdown in docs/ is written to render on GitHub, which wants inline math as $`...`$ and
display math in ```math blocks. MathJax (through pymdownx.arithmatex) wants \\(...\\) and
\\[...\\]. The hook converts one into the other before the page is parsed, so that the same
files render in both places. Pipes inside math become \\vert, because Python-Markdown would
otherwise split table cells at them.

Tables are included with an HTML comment, `<!-- include: results/summary/t4_growth.md -->`,
which GitHub hides and the hook replaces by the file (paths relative to the repository root).

Figures live in figures/ at the repository root, outside docs/. They are added to the site
as extra files under figures/. Pages link to them as on GitHub (../figures/F07_*.png from
docs/), and the hook rewrites the link relative to the page. Other relative links that leave
docs/ (source files, the report) are rewritten to their GitHub URL.
"""

from __future__ import annotations

import posixpath
import re
from pathlib import Path

from mkdocs.structure.files import File

ROOT = Path(__file__).resolve().parents[1]
GITHUB = "https://github.com/NugiePhysics/geodesic-integrators/blob/main/"

MATH_BLOCK = re.compile(r"^```math\n(.*?)^```\n", re.M | re.S)
INLINE_MATH = re.compile(r"\$`(.+?)`\$")
LINK = re.compile(r"(\]\()([^)\s#]+)(#[^)\s]*)?(\))")
INCLUDE = re.compile(r"^<!-- include: (\S+) -->$", re.M)


def _protect(tex: str) -> str:
    return tex.replace(r"\|", r"\Vert ").replace("|", r"\vert ")


def _block(match: re.Match) -> str:
    body = match.group(1).strip()
    if body.startswith(r"\begin{"):
        return f"\n{body}\n\n"
    return f"\n\\[\n{body}\n\\]\n\n"


def on_page_markdown(markdown: str, page, config, files) -> str:
    markdown = INCLUDE.sub(lambda m: (ROOT / m.group(1)).read_text(), markdown)
    markdown = MATH_BLOCK.sub(_block, markdown)
    markdown = INLINE_MATH.sub(lambda m: rf"\({_protect(m.group(1))}\)", markdown)

    here = posixpath.dirname(page.file.src_uri)

    def relink(m: re.Match) -> str:
        target = m.group(2)
        if re.match(r"^[a-z]+:", target) or target.startswith("/"):
            return m.group(0)
        resolved = posixpath.normpath(posixpath.join("docs", here, target))
        if resolved.startswith("docs/"):
            return m.group(0)
        if resolved.startswith("figures/"):
            # Served next to the pages (on_files): link relative to the page's own folder.
            target = posixpath.relpath(resolved, here or ".")
            return f"{m.group(1)}{target}{m.group(3) or ''}{m.group(4)}"
        return f"{m.group(1)}{GITHUB}{resolved}{m.group(3) or ''}{m.group(4)}"

    return LINK.sub(relink, markdown)


def on_files(files, config):
    for path in sorted((ROOT / "figures").glob("*.png")):
        files.append(
            File(
                f"figures/{path.name}",
                src_dir=str(ROOT),
                dest_dir=config["site_dir"],
                use_directory_urls=config["use_directory_urls"],
            )
        )
    return files
