import os
import re
import docx

TEX_PATH = r"C:\ZZ - Antigravity\Mestrado\01 - Artigos\Artigo 02\main.tex"
BIB_PATH = r"C:\ZZ - Antigravity\Mestrado\01 - Artigos\Artigo 02\cas-refs.bib"
DOCX_PATH = r"C:\ZZ - Antigravity\Mestrado\01 - Artigos\Artigo 02\Envio do Artigo\REMUNOM\Artigo_02_REMUNOM_Versao_Completa.docx"
SCRIPT_PATH = r"C:\ZZ - Antigravity\Mestrado\01 - Artigos\Artigo 02\scripts\gerar_artigo_remunom.py"

def audit_docx():
    print("\n========================================")
    print("2. AUDITORIA FORENSE DOCX (REMUNOM)")
    print("========================================")
    if not os.path.exists(DOCX_PATH):
        print(f"Arquivo DOCX não encontrado: {DOCX_PATH}")
        return

    doc = docx.Document(DOCX_PATH)
    paragraphs = []
    full_text = []
    for p in doc.paragraphs:
        txt = p.text.strip()
        if txt:
            paragraphs.append(txt)
            full_text.append(txt)

    for table in doc.tables:
        for row in table.rows:
            for cell in row.cells:
                txt = cell.text.strip()
                if txt:
                    full_text.append(txt)

    combined_text = "\n".join(full_text)

    # Locate References section
    ref_start_idx = -1
    for i, p in enumerate(paragraphs):
        if re.match(r'^(REFERÊNCIAS|REFERENCES|REFERENCIAS)\b', p, re.IGNORECASE):
            ref_start_idx = i
            break

    if ref_start_idx != -1:
        body_paras = paragraphs[:ref_start_idx]
        ref_paras = paragraphs[ref_start_idx+1:]
    else:
        print("AVISO: Seção de Referências não identificada pelo cabeçalho padrão!")
        body_paras = paragraphs
        ref_paras = []

    body_text = "\n".join(body_paras)
    print(f"Total de parágrafos no corpo do texto: {len(body_paras)}")
    print(f"Total de entradas identificadas na seção de Referências: {len(ref_paras)}")

    # Extract all text inside parentheses
    raw_parens = re.findall(r'\(([^)]+)\)', body_text)
    year_cites = []
    for rp in raw_parens:
        if re.search(r'\b(19\d{2}|20\d{2})\b', rp):
            year_cites.append(rp)

    print(f"\nParenthetical citation candidates found: {len(year_cites)}")
    for yc in year_cites:
        print(f"  - ({yc})")

    # Narrative citations: Word (Year)
    narrative_cites = re.findall(r'\b([A-Z][a-zA-Z\-]+(?:\s+et\s+al\.)?)\s*\((19\d{2}|20\d{2}[a-z]?)\)', body_text)
    print(f"\nNarrative citations found: {len(narrative_cites)}")
    for nc in narrative_cites:
        print(f"  - {nc[0]} ({nc[1]})")

    # Numeric citations [X]
    num_cites = re.findall(r'\[(\d+(?:\s*,\s*\d+)*)\]', body_text)
    print(f"\nNumeric citations found: {len(num_cites)}")
    for nc in num_cites[:10]:
        print(f"  - [{nc}]")

    print("\n--- Lista de Referências no DOCX ---")
    for idx, rp in enumerate(ref_paras, 1):
        print(f"Ref {idx:02d}: {rp[:100]}...")

    print("\n--- Cruzamento Detalhado: Referência vs. Presença no Corpo do Texto ---")
    orphan_refs_docx = []
    
    # We will check each of the 26 references individually
    for idx, rp in enumerate(ref_paras, 1):
        # Extract surname
        m_author = re.match(r'^([A-ZÀ-Ú\s\-]+),', rp)
        m_year = re.search(r'\b(19\d{2}|20\d{2})\b', rp)
        
        surname = m_author.group(1).strip() if m_author else rp.split()[0].replace(',', '')
        year = m_year.group(1) if m_year else "????"

        # Search tokens
        if surname.startswith("THALES"):
            tokens = ["Thales", "ASCEND"]
        elif surname.startswith("VAN DER"):
            tokens = ["van der Merwe", "Merwe", "Van der Merwe"]
        elif surname.startswith("ORTIZ"):
            tokens = ["Ortiz", "del Castillo"]
        else:
            tokens = [surname.split()[0].title(), surname.split()[0].upper()]

        found_in_body = []
        for t in tokens:
            matches = re.findall(rf'([^.\n]*\b{re.escape(t)}\b[^.\n]*)', body_text, re.IGNORECASE)
            if matches:
                found_in_body.extend(matches)

        if not found_in_body:
            orphan_refs_docx.append((idx, surname, year, rp))
            print(f"[ÓRFÃ - NUNCA CITADA NO CORPO] Ref {idx:02d}: {surname} ({year})")
        else:
            print(f"[OK CITADA - {len(found_in_body)}x] Ref {idx:02d}: {surname} ({year})")
            for fb in found_in_body[:2]:
                print(f"     Trecho: \"...{fb.strip()}...\"")

    print(f"\nTotal de referências na lista final: {len(ref_paras)}")
    print(f"Total de referências ÓRFÃS (na lista mas NÃO citadas no corpo): {len(orphan_refs_docx)}")

    # Check the inverse: Are there citations in body that are NOT in ref_paras?
    print("\n--- Cruzamento Inverso: Citações no Corpo vs. Lista de Referências ---")
    orphan_cites_in_body = []
    
    # Check parenthetical candidates
    for yc in year_cites:
        # Ignore things like (1980), purely mathematical coordinates, dates, etc.
        # Check if any reference matches
        matched = False
        for rp in ref_paras:
            # check if year matches and some word matches
            m_year = re.search(r'\b(19\d{2}|20\d{2})\b', rp)
            if m_year and m_year.group(1) in yc:
                # check author word
                m_author = re.match(r'^([A-ZÀ-Ú\s\-]+),', rp)
                author_word = m_author.group(1).split()[0].lower() if m_author else ""
                if author_word and author_word in yc.lower():
                    matched = True
                    break
                elif "thales" in yc.lower() and "thales" in rp.lower():
                    matched = True
                    break
                elif "ascend" in yc.lower() and "ascend" in rp.lower():
                    matched = True
                    break
        if not matched:
            # Check if it looks like an academic citation: Author, Year
            if re.search(r'[A-Za-z]+.*,\s*\b(19\d{2}|20\d{2})\b', yc):
                orphan_cites_in_body.append(yc)

    for nc in narrative_cites:
        matched = False
        author_name = nc[0].split()[0].lower()
        year = nc[1][:4]
        for rp in ref_paras:
            if year in rp and author_name in rp.lower():
                matched = True
                break
        if not matched:
            orphan_cites_in_body.append(f"{nc[0]} ({nc[1]})")

    print(f"Citações no corpo sem referência correspondente ({len(orphan_cites_in_body)}):")
    for oc in sorted(set(orphan_cites_in_body)):
        print(f"  [CITAÇÃO ÓRFÃ NO CORPO] -> {oc}")

    # Check for Figure and Table cross-references in DOCX
    fig_mentions = set(re.findall(r'\b(?:Figura|Fig\.|Figure)\s+(\d+)\b', body_text, re.IGNORECASE))
    fig_captions = set(re.findall(r'\b(?:Figura|Figure)\s+(\d+)\s*[:\.\-–]', combined_text, re.IGNORECASE))
    print(f"\nFiguras citadas no texto: {sorted(fig_mentions)}")
    print(f"Figuras com legenda encontradas: {sorted(fig_captions)}")
    missing_figs = fig_mentions - fig_captions
    if missing_figs:
        print(f"[ALERTA] Figuras citadas sem legenda no texto: {missing_figs}")
    else:
        print("[OK] Todas as figuras citadas possuem legenda correspondente.")

    tab_mentions = set(re.findall(r'\b(?:Tabela|Table|Tab\.)\s+(\d+)\b', body_text, re.IGNORECASE))
    tab_captions = set(re.findall(r'\b(?:Tabela|Table)\s+(\d+)\s*[:\.\-–]', combined_text, re.IGNORECASE))
    print(f"\nTabelas citadas no texto: {sorted(tab_mentions)}")
    print(f"Tabelas com legenda encontradas: {sorted(tab_captions)}")
    missing_tabs = tab_mentions - tab_captions
    if missing_tabs:
        print(f"[ALERTA] Tabelas citadas sem legenda no texto: {missing_tabs}")
    else:
        print("[OK] Todas as tabelas citadas possuem legenda correspondente.")

if __name__ == '__main__':
    audit_docx()
