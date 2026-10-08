# -*- coding: utf-8 -*-
"""
md_to_html_for_docx.py — 把装配母本 Markdown 转成 html-to-docx 可用的 HTML。

只做保真转换，不改任何文字。处理：标题(#~####)、段落、无序列表、
Markdown 表格、**粗体** / *斜体*、行内 `code`、> 引用、HTML 注释块(丢弃)、
水平线(丢弃)。输出 UTF-8 HTML，交给 html_to_docx convert。
"""
import os
import re
import sys
import html as _html

SRC = sys.argv[1]
DST = sys.argv[2]

CSS = """
body{font-family:"Times New Roman","SimSun",serif;font-size:11pt;line-height:1.6;color:#000;margin:28px 36px;}
h1{font-family:"Microsoft YaHei",sans-serif;font-size:15pt;color:#1F3864;margin:0 0 14px;}
h2{font-family:"Microsoft YaHei",sans-serif;font-size:12.5pt;color:#1F3864;margin:18px 0 6px;}
h3{font-family:"Microsoft YaHei",sans-serif;font-size:11pt;color:#2E5496;margin:14px 0 5px;}
h4{font-family:"Microsoft YaHei",sans-serif;font-size:10.5pt;color:#2E5496;margin:12px 0 4px;}
p{margin:5px 0;text-align:justify;}
ul{margin:6px 0 6px 20px;}
li{margin:3px 0;text-align:justify;}
blockquote{margin:8px 0 8px 14px;color:#444;font-size:10pt;}
table{border-collapse:collapse;width:100%;margin:6px 0 10px;font-size:9pt;}
th{background:#1F3864;color:#fff;border:0.5px solid #B8C4D9;padding:4px 6px;text-align:left;}
td{border:0.5px solid #B8C4D9;padding:4px 6px;vertical-align:top;}
tr:nth-child(even) td{background:#F7F9FC;}
code{font-family:Consolas,monospace;font-size:9pt;background:#F0F2F5;}
"""


def inline(t):
    t = _html.escape(t, quote=False)
    t = re.sub(r"`([^`]+)`", r"<code>\1</code>", t)
    t = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", t)
    t = re.sub(r"(?<!\*)\*([^*\n]+)\*(?!\*)", r"<em>\1</em>", t)
    return t


def split_row(line):
    line = line.strip()
    if line.startswith("|"):
        line = line[1:]
    if line.endswith("|"):
        line = line[:-1]
    return [c.strip() for c in line.split("|")]


def main():
    raw = open(SRC, encoding="utf-8").read()
    # 去掉装配注记 HTML 注释块
    raw = re.sub(r"<!--.*?-->", "", raw, flags=re.S)

    out = []
    lines = raw.split("\n")
    i = 0
    n = len(lines)
    while i < n:
        ln = lines[i]
        s = ln.strip()

        if not s:
            i += 1
            continue
        if s in ("---", "***", "___"):
            i += 1
            continue

        m = re.match(r"^(#{1,4})\s+(.*)$", s)
        if m:
            lvl = len(m.group(1))
            out.append(f"<h{lvl}>{inline(m.group(2))}</h{lvl}>")
            i += 1
            continue

        # 表格
        if s.startswith("|") and i + 1 < n and re.match(r"^\|[\s:|-]+\|$", lines[i + 1].strip()):
            header = split_row(lines[i])
            i += 2
            body = []
            while i < n and lines[i].strip().startswith("|"):
                body.append(split_row(lines[i]))
                i += 1
            t = ["<table><thead><tr>"]
            t += [f"<th>{inline(c)}</th>" for c in header]
            t.append("</tr></thead><tbody>")
            for r in body:
                t.append("<tr>" + "".join(f"<td>{inline(c)}</td>" for c in r) + "</tr>")
            t.append("</tbody></table>")
            out.append("".join(t))
            continue

        # 无序列表
        if re.match(r"^[-*+]\s+", s):
            items = []
            while i < n and re.match(r"^[-*+]\s+", lines[i].strip()):
                items.append(re.sub(r"^[-*+]\s+", "", lines[i].strip()))
                i += 1
            out.append("<ul>" + "".join(f"<li>{inline(x)}</li>" for x in items) + "</ul>")
            continue

        # 引用
        if s.startswith(">"):
            buf = []
            while i < n and lines[i].strip().startswith(">"):
                buf.append(lines[i].strip().lstrip(">").strip())
                i += 1
            out.append("<blockquote>" + " ".join(inline(x) for x in buf if x) + "</blockquote>")
            continue

        # 普通段落（含编号参考文献行——原样保留编号，不做成 list）
        buf = [s]
        i += 1
        while i < n:
            nx = lines[i].strip()
            if (not nx) or nx.startswith("#") or nx.startswith("|") or nx.startswith(">") \
               or re.match(r"^[-*+]\s+", nx) or nx in ("---", "***", "___"):
                break
            buf.append(nx)
            i += 1
        out.append("<p>" + " ".join(inline(x) for x in buf) + "</p>")

    html = ("<!DOCTYPE html><html lang=\"en\"><head><meta charset=\"utf-8\">"
            f"<style>{CSS}</style></head><body>" + "\n".join(out) + "</body></html>")
    open(DST, "w", encoding="utf-8").write(html)
    print("wrote", DST, "chars", len(html))


if __name__ == "__main__":
    main()
