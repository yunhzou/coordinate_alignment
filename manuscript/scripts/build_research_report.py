"""Render the review's small Markdown subset to a linked, printable PDF.

Uses the same browser installation as validate_viewer.py. No network resources
are loaded. The Markdown evidence file remains the editable source.
"""
import html
import os
from pathlib import Path
import re

from playwright.sync_api import sync_playwright

MAN = Path(__file__).resolve().parents[1]
NAME = 'review-led-novelty-20260910'
SOURCE = MAN / 'evidence' / f'{NAME}.md'
OUT = MAN / 'build'


def inline(value):
    tokens = []

    def keep(value):
        tokens.append(value)
        return f'\x00{len(tokens)-1}\x00'

    value = re.sub(r'`([^`]+)`', lambda m: keep('<code>' + html.escape(m[1]) + '</code>'), value)
    value = re.sub(r'\[([^\]]+)\]\(([^)]+)\)',
                   lambda m: keep(f'<a href="{html.escape(m[2], quote=True)}">{html.escape(m[1])}</a>'), value)
    value = re.sub(r'\[\^(\d+)\]',
                   lambda m: keep(f'<sup><a href="#note-{m[1]}">{m[1]}</a></sup>'), value)
    value = html.escape(value)
    value = re.sub(r'\*([^*]+)\*', r'<em>\1</em>', value)
    return re.sub(r'\x00(\d+)\x00', lambda m: tokens[int(m[1])], value)


def render(markdown):
    blocks = re.split(r'\n\s*\n|\n(?=\[\^\d+\]:)', markdown.strip())
    result = []
    for block in blocks:
        if block.startswith('[^'):
            match = re.fullmatch(r'\[\^(\d+)\]:\s*(.*)', block, re.S)
            if not match:
                raise ValueError(f'Malformed note: {block[:80]}')
            result.append(f'<p class="note" id="note-{match[1]}"><b>{match[1]}.</b> {inline(match[2])}</p>')
        elif block.startswith('#'):
            match = re.fullmatch(r'(#{1,2})\s+(.*)', block)
            if not match:
                raise ValueError('Unsupported heading')
            level = len(match[1])
            result.append(f'<h{level}>{inline(match[2])}</h{level}>')
        elif block.startswith('|'):
            rows = [line.strip().strip('|').split('|') for line in block.splitlines()]
            assert all(len(row) == len(rows[0]) for row in rows)
            head = ''.join(f'<th>{inline(c.strip())}</th>' for c in rows[0])
            body = ''.join('<tr>' + ''.join(f'<td>{inline(c.strip())}</td>' for c in row) + '</tr>' for row in rows[2:])
            result.append(f'<table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table>')
        elif block.startswith('> '):
            result.append('<blockquote>' + inline(block[2:]) + '</blockquote>')
        else:
            result.append('<p>' + inline(block) + '</p>')
    return '\n'.join(result)


CSS = '''
@page { size: A4; margin: 19mm 20mm; }
body { font: 10.5pt/1.42 "Liberation Serif", "Times New Roman", serif; color: #171717; }
h1 { font: 23pt/1.13 "Liberation Sans", Arial, sans-serif; margin: 0 0 17pt; }
h2 { font: bold 13pt/1.2 "Liberation Sans", Arial, sans-serif; margin: 20pt 0 8pt; break-after: avoid; }
p { margin: 0 0 9pt; orphans: 3; widows: 3; }
a { color: #222; text-decoration: underline; }
sup { line-height: 0; font-size: 8pt; margin-left: 1pt; }
sup a { text-decoration: none; }
sup + sup::before { content: ","; margin-right: 1pt; }
code { font: 8.8pt "Liberation Mono", monospace; overflow-wrap: anywhere; }
table { width: 100%; border-collapse: collapse; font-size: 9.3pt; line-height: 1.32; margin: 10pt 0 12pt; }
th, td { border-bottom: 0.6pt solid #aaa; vertical-align: top; text-align: left; padding: 6pt 5pt; }
th { background: #eee; font-family: "Liberation Sans", Arial, sans-serif; }
th:first-child, td:first-child { width: 21%; }
tr { break-inside: avoid; }
thead { display: table-header-group; }
blockquote { margin: 10pt 14pt; font-style: italic; }
.note { font-size: 8.7pt; line-height: 1.3; margin-bottom: 6pt; break-inside: avoid; overflow-wrap: anywhere; }
'''


def main():
    OUT.mkdir(exist_ok=True)
    content = render(SOURCE.read_text())
    destination = OUT / f'{NAME}.html'
    destination.write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
                           '<title>Continuous fragment growth: literature and mechanism assessment</title>'
                           f'<style>{CSS}</style><body>{content}</body></html>')
    with sync_playwright() as pw:
        cached = Path('/h/399/yunhengzou/.cache/ms-playwright/chromium-1243/chrome-linux64/chrome')
        executable = os.environ.get('MANUSCRIPT_BROWSER', str(cached) if cached.exists() else None)
        environment = dict(os.environ)
        libs = OUT / 'browser-libs/root/usr/lib/x86_64-linux-gnu'
        if libs.exists():
            environment['LD_LIBRARY_PATH'] = str(libs) + ':' + environment.get('LD_LIBRARY_PATH', '')
        browser = pw.chromium.launch(executable_path=executable, env=environment,
                                    headless=True, args=['--no-sandbox'])
        page = browser.new_page()
        page.goto(destination.as_uri(), wait_until='load')
        missing = page.eval_on_selector_all('a[href^="#"]',
            '(links) => links.filter(a => !document.getElementById(a.hash.slice(1))).map(a => a.hash)')
        assert not missing, missing
        page.pdf(path=str(MAN / 'evidence' / f'{NAME}.pdf'), prefer_css_page_size=True,
                 print_background=True, display_header_footer=False)
        browser.close()
    print(f'Rendered evidence/{NAME}.pdf')


if __name__ == '__main__':
    main()
