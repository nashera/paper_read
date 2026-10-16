#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""PaperRenamer CLI —— 学术论文PDF智能重命名工具（命令行版）

基于GUI版本(paper_renamer_gui.py)重写，元数据提取与合并逻辑保持一致：
PDF内嵌元数据 / 首页文本 / BibTeX / Crossref(DOI) / arXiv 多源融合，
按自定义模板批量重命名。

用法示例:
  python paper_renamer_cli.py paper.pdf
  python paper_renamer_cli.py ~/papers -r
  python paper_renamer_cli.py ~/papers -t "{year}_{author}_{title}" -n
  python paper_renamer_cli.py a.pdf b.pdf -b ~/bibs -v
"""

import argparse
import os
import re
import sys
from string import Formatter

import requests
import bibtexparser
from PyPDF2 import PdfReader
from pdfminer.high_level import extract_text

try:
    import pymupdf as fitz  # PyMuPDF 新版导入方式
except ImportError:
    import fitz  # 兼容旧版

VERSION = "1.0"
HEADERS = {"User-Agent": f"PaperRenamer-CLI/{VERSION}"}

# 元数据提取（与GUI版一致） ===================================================

def extract_pdf_metadata(pdf_path):
    meta = {}
    try:
        r = PdfReader(pdf_path)
        info = r.metadata
        meta.update({k.lower(): v for k, v in info.items() if isinstance(v, str)})
    except Exception: pass
    try:
        doc = fitz.open(pdf_path)
        meta.update({k.lower(): v for k, v in doc.metadata.items() if isinstance(v, str)})
        doc.close()
    except Exception: pass
    return meta

def extract_text_info(pdf_path):
    try:
        doc = fitz.open(pdf_path)
        text = doc[0].get_text()
        doc.close()
        return text
    except Exception:
        try:
            return extract_text(pdf_path, maxpages=1)
        except Exception:
            return ""

def extract_doi(text):
    m = re.search(r'\b10\.\d{4,9}/[-._;()/:A-Z0-9]+\b', text, re.I)
    return m.group(0) if m else ""

def extract_arxiv_id(text):
    m = re.search(r'arXiv:\s?(\d{4}\.\d{4,5})', text)
    return m.group(1) if m else ""

def extract_bibtex_metadata(bibtex_path):
    """解析.bib文件首条目，返回小写字段字典。兼容bibtexparser v1/v2。

    注：GUI版使用 v1 的 bibtexparser.load()，在 bibtexparser>=2.0 下会静默失败。
    """
    try:
        with open(bibtex_path, encoding='utf-8') as f:
            content = f.read()
        if hasattr(bibtexparser, "parse_string"):  # bibtexparser v2
            entry = bibtexparser.parse_string(content).entries[0]
            return {k.lower(): v.value for k, v in entry.fields_dict.items()}
        bib = bibtexparser.loads(content)  # bibtexparser v1
        entry = bib.entries[0]
        return {k.lower(): v for k, v in entry.items()}
    except Exception:
        return {}

def query_crossref(doi):
    url = f"https://api.crossref.org/works/{doi}"
    try:
        r = requests.get(url, timeout=8, headers=HEADERS)
        if r.status_code == 200:
            d = r.json()["message"]
            authors = "; ".join([f"{a.get('family','')}, {a.get('given','')}" for a in d.get("author",[])])
            short_journal = (d.get("short-container-title") or [""])[0].replace(".", "").strip()
            return {
                "title": d.get("title", [""])[0],
                "author": authors,
                "journal": d.get("container-title", [""])[0],
                "journal_abbr": short_journal,
                "year": str(d.get("published-print", d.get("created", {})).get("date-parts", [[None]])[0][0]),
                "volume": d.get("volume", ""),
                "issue": d.get("issue", ""),
                "pages": d.get("page", ""),
                "doi": doi,
                "issn": d.get("ISSN", [""])[0] if d.get("ISSN") else "",
                "publisher": d.get("publisher", "")
            }
    except Exception:
        return {}
    return {}

def query_arxiv(arxiv_id):
    url = f"https://export.arxiv.org/api/query?id_list={arxiv_id}"
    try:
        r = requests.get(url, timeout=8, headers=HEADERS)
        if r.status_code == 200:
            import xml.etree.ElementTree as ET
            root = ET.fromstring(r.content)
            entry = root.find('{http://www.w3.org/2005/Atom}entry')
            title = entry.find('{http://www.w3.org/2005/Atom}title').text.strip()
            authors = "; ".join([a.find('{http://www.w3.org/2005/Atom}name').text for a in entry.findall('{http://www.w3.org/2005/Atom}author')])
            year = entry.find('{http://www.w3.org/2005/Atom}published').text[:4]
            return {
                "title": title, "author": authors, "journal": "arXiv", "year": year, "arxiv_id": arxiv_id
            }
    except Exception:
        return {}
    return {}

def extract_candidates(text):
    lines = [l.strip() for l in text.split('\n') if l.strip()]
    title, author, journal, year, volume, issue, pages, issn = "", "", "", "", "", "", "", ""
    # 优先检查前5行，排除摘要/引言等章节标题
    title_candidates = []
    for i, l in enumerate(lines[:5]):
        if 10 < len(l) < 120 and not re.search(r'Abstract|Introduction|References|Keywords', l, re.I):
            title_candidates.append(l)
    # 若前5行无候选，扩展至前10行
    if not title_candidates:
        title_candidates = [l for l in lines[:10] if 10 < len(l) < 120]
    if title_candidates:
        title = max(title_candidates, key=lambda l: len(l))  # 选择最长行作为标题
    author_candidates = [l for l in lines if (',' in l or 'and' in l or 'et al' in l) and len(l) < 70]
    if author_candidates:
        author = author_candidates[0]
    for l in lines:
        if re.match(r'\d{4}', l): year = l[:4]
        if 'Vol' in l or 'Volume' in l: volume = l
        if 'Issue' in l or 'No.' in l: issue = l
        if re.search(r'pp\.|pages', l): pages = l
        if 'ISSN' in l: issn = l
        if re.match(r'^[A-Z][A-Za-z0-9 \-]{8,}$', l) and not journal: journal = l
    return dict(title=title, author=author, journal=journal, year=year, volume=volume, issue=issue, pages=pages, issn=issn)

# 人名与期刊缩写 ================================================================

# ISO 4 风格常用期刊词缩写表（离线回退用；联网时优先使用Crossref的short-container-title）
JOURNAL_ABBR = {
    'nature': 'Nat', 'sciences': 'Sci', 'science': 'Sci', 'scientific': 'Sci',
    'communications': 'Commun', 'letters': 'Lett', 'medicine': 'Med', 'medical': 'Med',
    'biology': 'Biol', 'biological': 'Biol', 'biomedicine': 'Biomed',
    'oncology': 'Oncol', 'oncologia': 'Oncol', 'oncologist': 'Oncologist',
    'neurology': 'Neurol', 'neuroscience': 'Neurosci', 'neurosciences': 'Neurosci',
    'neurosurgery': 'Neurosurg', 'neurosurgical': 'Neurosurg', 'neuro': 'Neuro',
    'neuropathology': 'Neuropathol', 'neuropathologica': 'Neuropathol',
    'pathology': 'Pathol', 'pathogenesis': 'Pathog', 'molecular': 'Mol',
    'cellular': 'Cell', 'genetics': 'Genet', 'genomics': 'Genomics',
    'therapy': 'Ther', 'therapeutics': 'Ther', 'transduction': 'Transduct',
    'targeted': 'Target', 'research': 'Res', 'journal': 'J', 'international': 'Int',
    'annals': 'Ann', 'annual': 'Annu', 'proceedings': 'Proc',
    'national': 'Natl', 'academy': 'Acad', 'reviews': 'Rev', 'review': 'Rev',
    'clinical': 'Clin', 'clinics': 'Clin', 'chemistry': 'Chem', 'physics': 'Phys',
    'applied': 'Appl', 'advances': 'Adv', 'advanced': 'Adv',
    'materials': 'Mater', 'environmental': 'Environ', 'publications': 'Publ',
    'publication': 'Publ', 'surgery': 'Surg', 'surgical': 'Surg',
    'chinese': 'Chin', 'japanese': 'Jpn', 'european': 'Eur', 'american': 'Am',
    'british': 'Br', 'english': 'Engl', 'england': 'Engl', 'new': 'N',
    'accounts': 'Acc', 'reports': 'Rep', 'development': 'Dev', 'developmental': 'Dev',
    'immunology': 'Immunol', 'immunity': 'Immun', 'inflammation': 'Inflamm',
    'microbiology': 'Microbiol', 'virology': 'Virol', 'bacteriology': 'Bacteriol',
    'biochemistry': 'Biochem', 'biophysics': 'Biophys', 'biotechnology': 'Biotechnol',
    'pharmacology': 'Pharmacol', 'toxicology': 'Toxicol', 'therapeutic': 'Ther',
    'epidemiology': 'Epidemiol', 'health': 'Health', 'psychology': 'Psychol',
    'psychiatry': 'Psychiatry', 'cardiology': 'Cardiol', 'cardiovascular': 'Cardiovasc',
    'endocrinology': 'Endocrinol', 'metabolism': 'Metab', 'metabolic': 'Metab',
    'infectious': 'Infect', 'diseases': 'Dis', 'disease': 'Dis',
    'cells': 'Cells', 'cell': 'Cell', 'systems': 'Syst', 'system': 'Syst',
    'technology': 'Technol', 'technologies': 'Technol', 'engineering': 'Eng',
    'computer': 'Comput', 'computers': 'Comput', 'computing': 'Comput',
    'intelligence': 'Intell', 'statistical': 'Stat', 'statistics': 'Stat',
    'mathematical': 'Math', 'mathematics': 'Math', 'methods': 'Methods',
    'analysis': 'Anal', 'analyses': 'Anal', 'diagnostics': 'Diagn',
    'radiology': 'Radiol', 'radiation': 'Radiat', 'transplantation': 'Transplant',
    'translational': 'Transl', 'experimental': 'Exp', 'general': 'Gen',
}
JOURNAL_STOPWORDS = {'of', 'and', 'the', 'in', 'for', 'on', 'a', 'an', 'to',
                     'und', 'der', 'die', 'das', 'et', 'des', 'du', 'la', 'le'}

def abbreviate_person(name):
    """人名缩写：'Cui, Huanhuan' -> 'Cui H'；'Ashish Vaswani'（显示名）-> 'Vaswani A'。"""
    name = name.strip()
    if not name:
        return ""
    if ',' in name:  # "Family, Given" 形式（Crossref/BibTeX）
        family, given = (p.strip() for p in name.split(',', 1))
        initials = ''.join(seg[0].upper() for seg in re.split(r'[\s\-]+', given) if seg)
        return f"{family} {initials}".strip()
    words = name.split()
    if len(words) < 2:  # 单词名无法拆分
        return name
    return f"{words[-1]} {''.join(w[0].upper() for w in words[:-1])}"

def abbreviate_journal(name):
    """ISO 4 风格期刊名缩写（内置词表+停用词），未收录的词保持原样。"""
    if not name:
        return ""
    out = []
    for w in name.replace('.', ' ').split():
        wl = w.lower().strip(',;')
        if not wl:
            continue
        if wl in JOURNAL_STOPWORDS:
            continue
        # 连字符词（如 Neuro-Oncology）分段缩写
        if '-' in w:
            parts = [JOURNAL_ABBR.get(seg.lower(), seg) for seg in w.strip(',;').split('-')]
            out.append('-'.join(parts))
        else:
            out.append(JOURNAL_ABBR.get(wl, w.strip(',;')))
    return ' '.join(out)

def get_complete_metadata(pdf_path, bibtex_path=None, online=True):
    """多源融合元数据，优先级从高到低：Crossref > arXiv > BibTeX > PDF内嵌 > 文本候选。"""
    meta = extract_pdf_metadata(pdf_path)
    text = extract_text_info(pdf_path)
    doi = extract_doi(text)
    arxiv_id = extract_arxiv_id(text)
    bibtex_meta = extract_bibtex_metadata(bibtex_path) if bibtex_path else {}
    candidates = extract_candidates(text)
    crossref_info = query_crossref(doi) if (doi and online) else {}
    arxiv_info = query_arxiv(arxiv_id) if (arxiv_id and online) else {}

    result = {}
    # 优先使用Crossref数据（如果DOI存在且查询成功）
    if crossref_info:
        result.update({k: v for k, v in crossref_info.items() if v})
    # 补充arXiv信息（仅当必要字段缺失时）
    if arxiv_info and not result.get('title'):
        result.update({k: v for k, v in arxiv_info.items() if v and k not in result})
    # 合并BibTeX元数据
    result.update({k: v for k, v in bibtex_meta.items() if v and k not in result})
    # 合并PDF元数据
    result.update({k: v for k, v in meta.items() if v and k not in result})
    # 使用文本提取的候选信息作为最后补充
    result.update({k: v for k, v in candidates.items() if v and k not in result})
    # 确保所有字段都存在
    for k in ['title', 'author', 'journal', 'year', 'volume', 'issue', 'pages', 'doi', 'issn', 'publisher', 'arxiv_id']:
        if k not in result:
            result[k] = ""
    # 第一作者：按 "; "（Crossref/arXiv）或 " and "（BibTeX）分隔取首位
    result['first_author'] = re.split(r';|\band\b', result['author'])[0].strip().rstrip(',') if result['author'] else ""
    # 简写字段：期刊优先用Crossref短标题（ISO 4），否则按内置词表缩写；人名取姓+名首字母
    result['first_author_abbr'] = abbreviate_person(result['first_author'])
    result['author_abbr'] = "; ".join(filter(None, (
        abbreviate_person(a) for a in re.split(r';|\s+and\s+', result['author']))))
    jabbr = result.get('journal_abbr', '').replace('.', '').strip()
    # Crossref短标题缺失或与全名相同时，回退到内置词表缩写
    if not jabbr or jabbr.lower() == result['journal'].strip().lower():
        jabbr = abbreviate_journal(result['journal'])
    # 全小写的短标题（个别出版社数据不规范）规范为首字母大写
    if jabbr and jabbr == jabbr.lower():
        jabbr = jabbr.title()
    result['journal_abbr'] = jabbr
    return result

# 模板与文件名 ================================================================

KNOWN_FIELDS = ('title', 'author', 'first_author', 'first_author_abbr', 'author_abbr',
                'journal', 'journal_abbr', 'year', 'volume', 'issue',
                'pages', 'doi', 'issn', 'publisher', 'arxiv_id')

def sanitize_filename(s):
    """净化文件名：非法字符与逗号替换为 _，空白替换为 _，合并连续 _ 并去除首尾分隔符。"""
    s = re.sub(r'[\\/:*?"<>|,]+', '_', s)  # Windows非法字符及逗号
    s = re.sub(r'\s+', '_', s)             # 空格、换行等空白
    s = re.sub(r'_{2,}', '_', s)           # 合并连续下划线（如 ", " 折叠为一个 _）
    return s.strip('._')

def smart_truncate(s, limit=128):
    """超长时优先在词边界（_ 或 -）处截断，避免单词被切一半（GUI版为硬截断）。"""
    if len(s) <= limit:
        return s
    cut = s[:limit]
    idx = max(cut.rfind('_'), cut.rfind('-'))
    if idx >= limit - 32:  # 边界离截断点不远时才采用，避免丢弃过多内容
        cut = cut[:idx]
    return cut.rstrip(' .,;:-_')

def parse_template_fields(template):
    """解析模板中的字段占位符并校验格式，返回字段名列表。"""
    fields = []
    for _, field, _, _ in Formatter().parse(template):
        if field is None:
            continue
        if not re.fullmatch(r'[A-Za-z_][A-Za-z0-9_]*', field):
            raise ValueError(f"模板字段 {{{field}}} 无效（仅支持简单字段名，如 {{title}}）")
        fields.append(field)
    return fields

def make_filename(meta, template="{author}-{year}-{journal}-{title}"):
    """按模板生成文件名：未知字段替换为空串，逗号/空格替换为 _ 并在词边界截断到128字符。"""
    fields = {f: meta.get(f, "") for f in parse_template_fields(template)}
    return smart_truncate(sanitize_filename(template.format(**fields)), 128)

# 路径收集与冲突处理 ==========================================================

def collect_pdfs(paths, recursive=False):
    """收集输入路径下的所有PDF（按输入顺序、目录内排序、去重）。

    返回 (pdf列表, 警告消息列表)。
    """
    pdfs, seen, warnings = [], set(), []

    def add(path):
        key = os.path.normcase(os.path.abspath(path))
        if key not in seen:
            seen.add(key)
            pdfs.append(path)

    for p in paths:
        if os.path.isdir(p):
            if recursive:
                found = [os.path.join(root, f) for root, _, files in os.walk(p)
                         for f in files if f.lower().endswith('.pdf')]
            else:
                found = [os.path.join(p, f) for f in os.listdir(p)
                         if f.lower().endswith('.pdf')]
            for f in sorted(found):
                add(f)
            if not found:
                warnings.append(f"文件夹中未找到PDF：{p}")
        elif os.path.isfile(p):
            if p.lower().endswith('.pdf'):
                add(p)
            else:
                warnings.append(f"忽略非PDF文件：{p}")
        else:
            warnings.append(f"路径不存在：{p}")
    return pdfs, warnings

def find_bibtex(pdf_path, bibtex_dir):
    """在BibTeX文件夹中查找与PDF同名的 .bib 文件（与GUI版一致）。"""
    base = os.path.splitext(os.path.basename(pdf_path))[0]
    path = os.path.join(bibtex_dir, base + ".bib")
    return path if os.path.exists(path) else None

def resolve_target(new_path, used):
    """目标名冲突时追加 _1、_2…后缀，避免覆盖磁盘已有文件或本批次已占用的名字。"""
    stem, ext = os.path.splitext(new_path)
    cand, i = new_path, 1
    while os.path.exists(cand) or os.path.normcase(os.path.abspath(cand)) in used:
        cand = f"{stem}_{i}{ext}"
        i += 1
        if i > 999:
            return None
    used.add(os.path.normcase(os.path.abspath(cand)))
    return cand

# 主流程 ======================================================================

def build_parser():
    parser = argparse.ArgumentParser(
        prog="paper_renamer_cli.py",
        description="学术论文PDF智能重命名工具（命令行版）：提取元数据并按模板批量重命名。",
        epilog=(
            "可用字段: " + ", ".join(KNOWN_FIELDS) + "\n"
            "（BibTeX条目中的其他字段也可直接引用，如 {keyword}）\n"
            "人名/期刊简写示例: first_author_abbr -> \"Cui H\"，journal_abbr -> \"Nat Commun\"\n\n"
            "示例:\n"
            "  %(prog)s paper.pdf\n"
            "  %(prog)s ~/papers -r -n\n"
            "  %(prog)s ~/papers -t \"{first_author_abbr}-{year}-{journal_abbr}-{title}\"\n"
            "  %(prog)s *.pdf -b ~/bibs -v"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('paths', nargs='+', metavar='PATH',
                        help='PDF文件或文件夹，可同时给出多个')
    parser.add_argument('-t', '--template', default='{author}-{year}-{journal}-{title}',
                        help='命名模板（默认: %(default)s）；文件名中的非法字符、逗号与空格自动替换为 _')
    parser.add_argument('-b', '--bibtex-dir', metavar='DIR',
                        help='BibTeX文件夹，自动匹配与PDF同名的 .bib 文件')
    parser.add_argument('-r', '--recursive', action='store_true',
                        help='递归处理子文件夹中的PDF')
    parser.add_argument('-n', '--dry-run', action='store_true',
                        help='预览重命名结果，不实际改名')
    parser.add_argument('--no-network', dest='online', action='store_false',
                        help='跳过Crossref/arXiv联网查询（离线模式，速度更快）')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='显示每个文件提取到的元数据')
    parser.add_argument('--version', action='version',
                        version=f'PaperRenamer CLI {VERSION}')
    return parser

def run(args):
    pdfs, warnings = collect_pdfs(args.paths, args.recursive)
    for w in warnings:
        print(f"[警告] {w}", file=sys.stderr)
    if not pdfs:
        print("[错误] 未找到任何PDF文件。", file=sys.stderr)
        return 1

    total = len(pdfs)
    used = set()
    ok = fail = skip = 0
    for idx, pdf_path in enumerate(pdfs, 1):
        prefix = f"[{idx}/{total}]"
        try:
            bibtex_path = find_bibtex(pdf_path, args.bibtex_dir) if args.bibtex_dir else None
            meta = get_complete_metadata(pdf_path, bibtex_path, online=args.online)
            if args.verbose:
                shown = " | ".join(f"{k}={v}" for k, v in meta.items() if v)
                print(f"{prefix} 元数据: {shown or '(空)'}")
                if bibtex_path:
                    print(f"{prefix} BibTeX: {bibtex_path}")

            new_stem = make_filename(meta, args.template)
            old_name = os.path.basename(pdf_path)
            if not new_stem:
                skip += 1
                print(f"{prefix} {old_name} —— 跳过（模板字段全为空，无法生成文件名）")
                continue

            new_name = new_stem + ".pdf"
            new_path = os.path.join(os.path.dirname(pdf_path), new_name)
            if os.path.normcase(os.path.abspath(pdf_path)) == os.path.normcase(os.path.abspath(new_path)):
                skip += 1
                print(f"{prefix} {old_name} —— 跳过（文件名未变化）")
                continue

            # 当前文件名已带冲突后缀（如上次运行生成的 "xxx_1"）且目标名仍被占用：
            # 保持现状，避免重复运行时后缀来回跳动（_1→_2→_1…）
            suffixless = re.fullmatch(r'(.+)_\d+$', os.path.splitext(old_name)[0])
            target_taken = (os.path.exists(new_path)
                            or os.path.normcase(os.path.abspath(new_path)) in used)
            if suffixless and suffixless.group(1) == new_stem and target_taken:
                skip += 1
                print(f"{prefix} {old_name} —— 跳过（已带冲突后缀，目标名被占用）")
                continue

            target = resolve_target(new_path, used)
            if target is None:
                fail += 1
                print(f"{prefix} {old_name} -> {new_name} —— 失败（目标名冲突过多）")
                continue
            target_name = os.path.basename(target)

            if args.dry_run:
                ok += 1
                print(f"{prefix} {old_name} -> {target_name} —— 预览（未改名）")
            else:
                try:
                    os.rename(pdf_path, target)
                    ok += 1
                    print(f"{prefix} {old_name} -> {target_name} —— 成功")
                except OSError as e:
                    fail += 1
                    used.discard(os.path.normcase(os.path.abspath(target)))
                    print(f"{prefix} {old_name} -> {target_name} —— 失败（{e}）")
        except Exception as e:
            fail += 1
            print(f"{prefix} {os.path.basename(pdf_path)} —— 失败（{e}）")

    if args.dry_run:
        print(f"\n预览完成：可重命名 {ok}，失败 {fail}，跳过 {skip}（共 {total} 个文件，未实际改名）")
    else:
        print(f"\n处理完成：成功 {ok}，失败 {fail}，跳过 {skip}（共 {total} 个文件）")
    return 0 if fail == 0 else 1

def main(argv=None):
    args = build_parser().parse_args(argv)
    if not args.template.strip():
        print("[错误] 命名模板不能为空。", file=sys.stderr)
        return 2
    try:
        parse_template_fields(args.template)
    except ValueError as e:
        print(f"[错误] {e}", file=sys.stderr)
        return 2
    try:
        return run(args)
    except KeyboardInterrupt:
        print("\n已取消。", file=sys.stderr)
        return 130

if __name__ == "__main__":
    sys.exit(main())
