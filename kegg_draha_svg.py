#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
KEGG dráha -> vektorové SVG se strukturními vzorci
==================================================

Funguje pro libovolnou dráhu z KEGG: stačí zadat KGML (XML) soubor,
např. gga00360.xml stažený z https://www.kegg.jp/kegg-bin/show_pathway?gga00360
(tlačítko „Download KGML“).

Spuštění
--------
    pip install rdkit matplotlib
    python kegg_draha_svg.py                 # okno s výběrem metabolitů a reakcí
    python kegg_draha_svg.py gga00360.xml    # rovnou SVG z příkazové řádky
    python kegg_draha_svg.py gga00360.xml --vynechat C00423,C05629 --pdf

Co program dělá
---------------
1. Načte KGML. Je-li to mapa organismu (např. gga), stáhne z KEGG referenční
   mapu ko<číslo> (všechny reakce dráhy), rovnice/EC čísla reakcí, názvy
   a strukturní vzorce (MOL) metabolitů. Vše se ukládá do složky kegg_cache/
   vedle skriptu, takže podruhé už to funguje i offline.
   Pro map00360 (fenylalanin) jsou data vestavěná a pracuje i bez internetu.
2. V okně si vybereš metabolity a reakce. Když odškrtneš metabolit uprostřed
   dráhy, jeho sousedé se spojí tečkovanou fialovou šipkou „vynechány N kroky“.
   Odškrtnutá reakce se prostě nekreslí.
3. Rozložení: vrstvy podle vzdálenosti od výchozího metabolitu (Sugiyama),
   minimalizace křížení, pravoúhlé šipky s vlastní „kolejí“ pro každou hranu,
   popisky enzymů s kontrolou kolizí.
4. Výstup je čisté SVG: struktury, šipky i veškerý text jsou křivky (<path>).
"""

import argparse
import math
import os
import re
import shutil
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET
from collections import defaultdict, deque

from rdkit import Chem, RDLogger
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from matplotlib.textpath import TextPath
from matplotlib.font_manager import FontProperties

RDLogger.DisableLog("rdApp.*")
if getattr(sys, "frozen", False):          # zabalená aplikace (.exe z PyInstalleru)
    SLOZKA = os.path.dirname(sys.executable)
    ZDROJE = getattr(sys, "_MEIPASS", SLOZKA)
else:
    SLOZKA = os.path.dirname(os.path.abspath(__file__))
    ZDROJE = SLOZKA

# =============================================================================
# 1) VESTAVĚNÉ ZNALOSTI (české názvy, SMILES, kompletní map00360)
# =============================================================================
# České názvy podle KEGG C-čísel (jiné látky dostanou anglický název z KEGG).
CESKE_NAZVY = {
    "C00079": "L-fenylalanin", "C00082": "L-tyrosin", "C02265": "D-fenylalanin",
    "C03519": "N-acetyl-L-fenylalanin", "C00166": "fenylpyruvát", "C05332": "fenylethylamin",
    "C02505": "2-fenylacetamid", "C00423": "trans-cinnamát", "C00601": "fenylacetaldehyd",
    "C05607": "fenyllaktát", "C02763": "2-hydroxy-3-fenylpropenoát",
    "C23049": "kyselina polyporová", "C23116": "askokorynin", "C07086": "fenylacetát",
    "C05853": "2-fenylethanol", "C05852": "2-hydroxyfenylacetát",
    "C00642": "4-hydroxyfenylacetát", "C05593": "3-hydroxyfenylacetát",
    "C06207": "2,6-dihydroxyfenylacetát", "C00582": "fenylacetyl-CoA",
    "C05598": "fenylacetylglycin", "C04148": "fenylacetylglutamin",
    "C15524": "fenylglyoxylyl-CoA", "C02137": "fenylglyoxylát", "C00512": "benzoyl-CoA",
    "C20062": "2-(1,2-epoxy-1,2-dihydrofenyl)acetyl-CoA",
    "C19975": "2-oxepin-2(3H)-ylidenacetyl-CoA",
    "C19946": "3-oxo-5,6-dehydrosuberyl-CoA-semialdehyd",
    "C19945": "3-oxo-5,6-dehydrosuberyl-CoA", "C14144": "2,3-didehydroadipyl-CoA",
    "C14145": "(3S)-3-hydroxyadipyl-CoA", "C02232": "3-oxoadipyl-CoA",
    "C00091": "sukcinyl-CoA", "C05629": "fenylpropanoát", "C01772": "trans-2-hydroxycinnamát",
    "C12621": "trans-3-hydroxycinnamát", "C12622": "cinnamát-cis-dihydrodiol",
    "C11588": "fenylpropanoát-cis-dihydrodiol", "C01198": "3-(2-hydroxyfenyl)propanoát",
    "C11457": "3-(3-hydroxyfenyl)propanoát", "C12623": "trans-2,3-dihydroxycinnamát",
    "C04044": "3-(2,3-dihydroxyfenyl)propanoát",
    "C12624": "2-hydroxy-6-oxonona-2,4,7-trienedioát",
    "C04479": "2-hydroxy-6-oxonona-2,4-dien-1,9-dioát", "C00596": "2-hydroxypenta-2,4-dienoát",
    "C03589": "4-hydroxy-2-oxopentanoát", "C00084": "acetaldehyd", "C00024": "acetyl-CoA",
    "C00122": "fumarát", "C00042": "sukcinát", "C00022": "pyruvát",
}
# Strukturní vzorce (SMILES). "[*:1]" se vykreslí jako zkratka CoA.
SMILES = {
    "C00079": "N[C@@H](Cc1ccccc1)C(=O)O", "C00082": "N[C@@H](Cc1ccc(O)cc1)C(=O)O",
    "C02265": "N[C@H](Cc1ccccc1)C(=O)O", "C03519": "CC(=O)N[C@@H](Cc1ccccc1)C(=O)O",
    "C00166": "O=C(O)C(=O)Cc1ccccc1", "C05332": "NCCc1ccccc1", "C02505": "NC(=O)Cc1ccccc1",
    "C00423": "OC(=O)/C=C/c1ccccc1", "C00601": "O=CCc1ccccc1", "C05607": "OC(=O)C(O)Cc1ccccc1",
    "C02763": "OC(=O)/C(O)=C/c1ccccc1",
    "C23049": "O=C1C(O)=C(c2ccccc2)C(=O)C(O)=C1c1ccccc1",
    "C23116": "O=C1C(O)=C(c2ccc(O)cc2)C(=O)C(O)=C1c1ccccc1",
    "C07086": "OC(=O)Cc1ccccc1", "C05853": "OCCc1ccccc1", "C05852": "OC(=O)Cc1ccccc1O",
    "C00642": "OC(=O)Cc1ccc(O)cc1", "C05593": "OC(=O)Cc1cccc(O)c1",
    "C06207": "OC(=O)Cc1c(O)cccc1O", "C00582": "O=C(Cc1ccccc1)S[*:1]",
    "C05598": "OC(=O)CNC(=O)Cc1ccccc1", "C04148": "NC(=O)CC[C@H](NC(=O)Cc1ccccc1)C(=O)O",
    "C15524": "O=C(S[*:1])C(=O)c1ccccc1", "C02137": "OC(=O)C(=O)c1ccccc1",
    "C00512": "O=C(S[*:1])c1ccccc1", "C20062": "O=C(CC12OC1C=CC=C2)S[*:1]",
    "C19975": "O=C(/C=C1/CC=CC=CO1)S[*:1]", "C19946": "O=CC/C=C\\CC(=O)CC(=O)S[*:1]",
    "C19945": "OC(=O)C/C=C\\CC(=O)CC(=O)S[*:1]", "C14144": "OC(=O)CC/C=C/C(=O)S[*:1]",
    "C14145": "OC(=O)CC[C@H](O)CC(=O)S[*:1]", "C02232": "OC(=O)CCC(=O)CC(=O)S[*:1]",
    "C00091": "OC(=O)CCC(=O)S[*:1]", "C05629": "OC(=O)CCc1ccccc1",
    "C01772": "OC(=O)/C=C/c1ccccc1O", "C12621": "OC(=O)/C=C/c1cccc(O)c1",
    "C12622": "OC(=O)/C=C/C1=CC=C[C@H](O)[C@@H]1O", "C11588": "OC(=O)CCC1=CC=C[C@H](O)[C@@H]1O",
    "C01198": "OC(=O)CCc1ccccc1O", "C11457": "OC(=O)CCc1cccc(O)c1",
    "C12623": "OC(=O)/C=C/c1cccc(O)c1O", "C04044": "OC(=O)CCc1cccc(O)c1O",
    "C12624": "OC(=O)/C=C/C(=O)/C=C/C=C(\\O)C(=O)O", "C04479": "OC(=O)CCC(=O)/C=C/C=C(\\O)C(=O)O",
    "C00596": "C=C/C=C(\\O)C(=O)O", "C03589": "CC(O)CC(=O)C(=O)O", "C00084": "CC=O",
    "C00024": "CC(=O)S[*:1]", "C00122": "OC(=O)/C=C/C(=O)O", "C00042": "OC(=O)CCC(=O)O",
    "C00022": "CC(=O)C(=O)O",
}
ZKRATKY = {1: "CoA"}

# Kompletní map00360 (ověřeno proti KEGG): (reakce, z, do, EC, typ)
# typ "->" nevratná z→do, "<=>" vratná, "->+" vedlejší produkt téže reakce
# (plná šipka bez popisku). EC "" = enzym neznámý (čárkovaně).
VESTAVENE_REAKCE = {"00360": [
    ("R01795", "C00079", "C00082", "1.14.16.1", "->"),
    ("R00686", "C00079", "C02265", "5.1.1.11", "<=>"),
    ("R00693", "C00079", "C03519", "2.3.1.53", "->"),
    ("R00688 R00689 R00692 R00694", "C00079", "C00166",
     "1.4.1.20 · 1.4.3.2 · 2.6.1.1 · 2.6.1.5 · 2.6.1.9 · 2.6.1.57 · 2.6.1.58", "<=>"),
    ("R00699", "C00079", "C05332", "4.1.1.28 · 4.1.1.53", "->"),
    ("R11918", "C00079", "C00601", "4.1.1.109", "->"),
    ("R00690 R00698", "C00079", "C02505", "1.13.12.9 · 1.11.1.21", "->"),
    ("R00697", "C00079", "C00423", "4.3.1.24 · 4.3.1.25", "->"),
    ("R01374 R01582", "C02265", "C00166", "1.4.5.1 · 2.6.1.21", "->"),
    ("R01370 R01371", "C00166", "C05607", "1.1.1.110 · 1.1.1.237", "<=>"),
    ("R01378", "C00166", "C02763", "5.3.2.1", "<=>"),
    ("R01377", "C00166", "C00601", "4.1.1.43 · 4.1.1.-", "->"),
    ("R01372", "C00166", "C05852", "1.13.11.27", "->"),
    ("R13509", "C00166", "C23049", "6.4.1.11", "->"),
    ("R13613", "C23049", "C23116", "1.14.14.196", "->"),
    ("R02612 R02613", "C05332", "C00601", "1.4.3.4 · 1.4.3.21 · 1.4.9.2", "->"),
    ("R02611", "C00601", "C05853", "1.1.1.90", "<=>"),
    ("R02536 R02537", "C00601", "C07086", "1.2.1.5 · 1.2.1.39", "->"),
    ("R02540", "C02505", "C07086", "3.5.1.4", "->"),
    ("R05487", "C07086", "C05852", "1.14.14.54", "->"),
    ("R06790", "C07086", "C00642", "", "->"),
    ("R06791", "C07086", "C05593", "", "->"),
    ("R05001", "C05852", "C06207", "", "->"),
    ("R02539", "C07086", "C00582", "6.2.1.30", "->"),
    ("R09840", "C00582", "C07086", "3.1.2.-", "->"),
    ("R05841", "C00582", "C05598", "2.3.1.192", "->"),
    ("R02576", "C00582", "C04148", "2.3.1.14", "->"),
    ("R07222", "C00582", "C15524", "1.17.5.1", "->"),
    ("R07294", "C15524", "C02137", "3.1.2.25", "->"),
    ("R02450", "C02137", "C00512", "1.2.1.58", "->"),
    ("R09838", "C00582", "C20062", "1.14.13.149", "->"),
    ("R09837", "C20062", "C19975", "5.3.3.18", "->"),
    ("R09836", "C19975", "C19946", "3.3.2.12", "->"),
    ("R09820", "C19946", "C19945", "1.2.1.91", "->"),
    ("R09839", "C19945", "C14144", "2.3.1.223", "->"),
    ("R06942", "C14144", "C14145", "4.2.1.17", "<=>"),
    ("R06941", "C14145", "C02232", "1.1.1.157", "<=>"),
    ("R00829", "C02232", "C00091", "2.3.1.174", "<=>"),
    ("R02252", "C00423", "C05629", "1.3.1.31", "->"),
    ("R02254", "C00423", "C01772", "1.14.13.14", "->"),
    ("R06781", "C00423", "C12621", "", "->"),
    ("R06783", "C00423", "C12622", "1.14.12.19", "->"),
    ("R06785", "C12622", "C12623", "1.3.1.87", "->"),
    ("R06787", "C12621", "C12623", "1.14.13.127", "->"),
    ("R03709", "C01772", "C01198", "1.3.1.11", "->"),
    ("R04899", "C05629", "C01198", "", "->"),
    ("R06780", "C05629", "C11457", "", "->"),
    ("R06782", "C05629", "C11588", "1.14.12.19", "->"),
    ("R06784", "C11588", "C04044", "1.3.1.87", "->"),
    ("R06786", "C11457", "C04044", "1.14.13.127", "->"),
    ("R03369", "C01198", "C04044", "1.14.13.4", "->"),
    ("R06788", "C12623", "C12624", "1.13.11.16", "->"),
    ("R04376", "C04044", "C04479", "1.13.11.16", "->"),
    ("R06789", "C12624", "C00596", "3.7.1.14", "->"),
    ("R06789", "C12624", "C00122", "", "->+"),
    ("R02603", "C04479", "C00596", "3.7.1.14", "->"),
    ("R02603", "C04479", "C00042", "", "->+"),
    ("R02601", "C00596", "C03589", "4.2.1.80", "->"),
    ("R00750", "C03589", "C00084", "4.1.3.39", "->"),
    ("R00750", "C03589", "C00022", "", "->+"),
    ("R00228", "C00084", "C00024", "1.2.1.10", "->"),
]}
# Navazující dráhy pro vestavěné mapy: látka -> [(číslo mapy, směr)]
VESTAVENE_ODKAZY = {"00360": {
    "C00079": [("00400", "<-"), ("00960", "->")], "C00082": [("00350", "->")],
    "C05852": [("00643", "->")], "C00512": [("00362", "->")], "C00091": [("00020", "->")],
    "C00024": [("00020", "->")], "C00122": [("00020", "->")], "C00042": [("00020", "->")],
    "C00022": [("00010", "->")],
}}
VESTAVENY_START = {"00360": "C00079"}
# Poloha látek v mapě KEGG (y, x) – jen pomocné pořadí pro rozložení, když chybí KGML.
VESTAVENE_SOURADNICE = {"00360": {"C00022": (712, 898), "C00024": (703, 705), "C00042": (577, 898), "C00079": (166, 415), "C00082": (159, 67), "C00084": (736, 842), "C00091": (810, 456), "C00122": (636, 762), "C00166": (304, 415), "C00423": (219, 762), "C00512": (598, 587), "C00582": (598, 271), "C00596": (604, 842), "C00601": (398, 267), "C00642": (504, 193), "C01198": (284, 1004), "C01772": (218, 1004), "C02137": (598, 492), "C02232": (733, 456), "C02265": (166, 571), "C02505": (483, 143), "C02763": (421, 608), "C03519": (249, 334), "C03589": (669, 842), "C04044": (448, 842), "C04148": (717, 271), "C04479": (528, 842), "C05332": (286, 267), "C05593": (543, 193), "C05598": (717, 197), "C05607": (383, 478), "C05629": (288, 842), "C05852": (486, 415), "C05853": (445, 374), "C06207": (536, 415), "C07086": (483, 267), "C11457": (293, 936), "C11588": (370, 842), "C12621": (320, 709), "C12622": (370, 762), "C12623": (448, 762), "C12624": (528, 762), "C14144": (733, 645), "C14145": (733, 550), "C15524": (598, 392), "C19945": (655, 645), "C19946": (655, 551), "C19975": (655, 457), "C20062": (655, 363), "C23049": (403, 557), "C23116": (464, 557)}}

CESKE_MAPY = {
    "00350": "Metabolismus tyrosinu",
    "00400": "Biosyntéza fenylalaninu,\ntyrosinu a tryptofanu",
    "00960": "Biosyntéza tropanových,\npiperidinových a pyridinových alkaloidů",
    "00643": "Degradace styrenu", "00362": "Degradace benzoátu",
    "00020": "Citrátový cyklus", "00010": "Glykolýza / glukoneogeneze",
    "00360": "Metabolismus fenylalaninu", "00380": "Metabolismus tryptofanu",
    "00250": "Metabolismus Ala, Asp a Glu", "00260": "Metabolismus Gly, Ser a Thr",
    "00270": "Metabolismus Cys a Met", "00330": "Metabolismus Arg a Pro",
    "00220": "Biosyntéza argininu", "00620": "Metabolismus pyruvátu",
    "00630": "Metabolismus glyoxylátu a dikarboxylátů", "00071": "Degradace mastných kyselin",
    "00280": "Degradace Val, Leu a Ile", "00290": "Biosyntéza Val, Leu a Ile",
    "00310": "Degradace lysinu", "00340": "Metabolismus histidinu",
    "00230": "Metabolismus purinů", "00240": "Metabolismus pyrimidinů",
    "00030": "Pentózafosfátová dráha", "00052": "Metabolismus galaktózy",
    "00051": "Metabolismus fruktózy a manózy", "00500": "Metabolismus škrobu a sacharózy",
    "00190": "Oxidativní fosforylace", "00910": "Metabolismus dusíku",
    "00640": "Metabolismus propanoátu", "00650": "Metabolismus butanoátu",
}

# KEGG označuje tyto reakce jako vratné, fyziologicky jsou nevratné.
NEVRATNE = {"R02613", "R02612", "R02536", "R02537"}


# =============================================================================
# 2) KEGG REST (se složkou mezipaměti)
# =============================================================================
CACHE = os.path.join(SLOZKA, "kegg_cache")


def kegg_get(cesta, online=True):
    """Text z https://rest.kegg.jp/<cesta>; mezipaměť v kegg_cache/. Chyba -> None."""
    import urllib.request
    os.makedirs(CACHE, exist_ok=True)
    soubor = os.path.join(CACHE, re.sub(r"[^\w.+-]", "_", cesta) + ".txt")
    if os.path.exists(soubor):
        with open(soubor, encoding="utf-8") as f:
            txt = f.read()
        return txt or None
    if not online:
        return None
    try:
        with urllib.request.urlopen(f"https://rest.kegg.jp/{cesta}", timeout=20) as r:
            txt = r.read().decode("utf-8")
    except Exception:
        return None
    with open(soubor, "w", encoding="utf-8") as f:
        f.write(txt)
    return txt or None


def kegg_zaznamy(ids, online=True):
    """Stáhne ploché záznamy KEGG (po 10) a vrátí {id: text záznamu}."""
    out = {}
    ids = sorted(set(ids))
    for k in range(0, len(ids), 10):
        davka = ids[k:k + 10]
        txt = kegg_get("get/" + "+".join(davka), online)
        if not txt:
            continue
        for zaznam in txt.split("///"):
            m = re.search(r"^ENTRY\s+(\S+)", zaznam, re.M)
            if m:
                out[m.group(1)] = zaznam
    return out


def pole_zaznamu(zaznam, klic):
    """Obsah pole (i víceřádkového) z plochého záznamu KEGG."""
    m = re.search(rf"^{klic}\s+(.*(?:\n {{12}}.*)*)", zaznam, re.M)
    return m.group(1) if m else ""


# =============================================================================
# 3) NAČTENÍ DRÁHY
# =============================================================================
# Přenašeče a kofaktory: v obrázku se ve výchozím stavu nekreslí (KEGG je u některých
# drah, např. citrátového cyklu, uvádí jako samostatné látky).
PRENASECE = {
    "C00068",                       # thiamindifosfát (ThPP)
    "C15972", "C15973",             # enzym-N6-(lipoyl)lysin, enzym-N6-(dihydrolipoyl)lysin
    "C00248", "C00579",             # lipoamid, dihydrolipoamid
    "C00229",                       # acyl-carrier protein (ACP)
    "C00010", "C00003", "C00004", "C00005", "C00006", "C00016", "C01352",   # CoA, NAD(P)(H), FAD(H2)
    "C00002", "C00008", "C00020", "C00009", "C00013",                       # ATP, ADP, AMP, Pi, PPi
    "C00001", "C00011", "C00007", "C00080", "C00014",                       # H2O, CO2, O2, H+, NH3
    "C00138", "C00139", "C00342", "C00343",                                 # ferredoxin, thioredoxin
}
RE_PRENASEC = re.compile(
    r"^(thiamine? diphosphate|enzyme n6-\((dihydro)?lipoyl\)lysine|(dihydro)?lipoamide(-e)?|"
    r"acyl-carrier protein|(reduced |oxidized )?(ferredoxin|thioredoxin|glutaredoxin|flavodoxin))$",
    re.I)
# Meziprodukty vázané na enzymový komplex (PDH, OGDH, syntáza mastných kyselin …):
# když jsou vypnuté, jejich řetězec se sloučí do jedné běžné šipky.
RE_VAZANY = re.compile(r"thiamin|thpp|lipoyl|lipoamide|^\[|^enzyme |carrier protein|\[acp\]|-acp\b", re.I)


VAZANE = {
    "C05125", "C05381",             # 2-(α-hydroxyethyl)-ThPP, 3-karboxy-1-hydroxypropyl-ThPP
    "C16254", "C16255",             # S-acetyl- a S-sukcinyldihydrolipoyllysin (E2 PDH / OGDH)
    "C15974", "C16239", "C16240",   # 3-methyl-1-hydroxybutyl-ThPP, S-acyldihydrolipoamidy (BCKDH)
}


def role_latky(c, nazev_en):
    """'prenasec' (ThPP, lipoyllysin, NAD …), 'vazany' (meziprodukt na enzymu) nebo 'normal'."""
    if c in PRENASECE or RE_PRENASEC.match(nazev_en or ""):
        return "prenasec"
    if c in VAZANE or RE_VAZANY.search(nazev_en or ""):
        return "vazany"
    return "normal"


class Draha:
    """Načtená dráha: metabolity, reakce (dvojice látek), odkazy na další mapy."""

    def __init__(self):
        self.cislo = ""
        self.org = ""
        self.titul = ""
        self.metabolity = {}    # cid -> dict(nazev, zdroj, poradi)
        self.reakce = []        # dict(rn, a, b, ec, typ, geny)
        self.odkazy = []        # (cid, číslo mapy, název mapy, směr)
        self.start = None
        self.zdroj_dat = ""
        self.online = True


def _parse_kgml(koren):
    entries = {e.get("id"): e for e in koren.findall("entry")}
    cpd = {}
    souradnice = {}
    for eid, e in entries.items():
        if e.get("type") == "compound":
            c = re.sub(r"^\w+:", "", e.get("name").split()[0])   # cpd:C…, gl:G…, dr:D…
            cpd[eid] = c
            g = e.find("graphics")
            if g is not None and g.get("x"):
                xy = (float(g.get("y")), float(g.get("x")))
                souradnice[c] = min(souradnice.get(c, xy), xy)
    return entries, cpd, souradnice


def nacti_drahu(cesta, online=True, log=print):
    d = Draha()
    if cesta is None:               # bez souboru: vestavěná mapa 00360 (bez genů organismu)
        koren = ET.fromstring('<pathway org="map" number="00360" title="Phenylalanine metabolism"/>')
    else:
        koren = ET.parse(cesta).getroot()
    d.cislo = koren.get("number", "")
    d.org = koren.get("org", "")
    d.titul = koren.get("title", "")
    entries, cpd, souradnice = _parse_kgml(koren)
    log(f"KGML: {d.org}{d.cislo} – {d.titul}")

    # geny organismu: reakce -> názvy genů
    geny = defaultdict(list)
    for e in entries.values():
        if e.get("type") != "gene" or not e.get("reaction"):
            continue
        g = e.find("graphics").get("name")
        nazev = g.split(",")[0].replace("...", "").strip()
        if len(e.get("name").split()) > 1 or g.endswith("..."):
            nazev += "…"
        for rn in e.get("reaction").split():
            rn = rn.replace("rn:", "")
            if nazev not in geny[rn]:
                geny[rn].append(nazev)

    ec_rn = {}
    for c, xy in VESTAVENE_SOURADNICE.get(d.cislo, {}).items():
        souradnice.setdefault(c, xy)
    if d.cislo in VESTAVENE_REAKCE:
        d.zdroj_dat = "vestavěná kompletní mapa"
        for rn, a, b, ec, typ in VESTAVENE_REAKCE[d.cislo]:
            d.reakce.append(dict(rn=rn.split(), a=a, b=b, ec=[x for x in ec.split(" · ") if x], typ=typ))
        ref = None
    else:
        # referenční mapa (všechny reakce): buď je to už ko/map KGML, nebo ji stáhneme
        ma_reakce_ko = any(entries.get(r.get("id")) is not None and
                           entries[r.get("id")].get("type") == "ortholog"
                           for r in koren.findall("reaction"))
        ref = koren
        d.zdroj_dat = "KGML"
        if not ma_reakce_ko and d.org not in ("ko", "map", "ec", "rn"):
            log(f"Stahuji referenční mapu ko{d.cislo} z KEGG…")
            txt = kegg_get(f"get/ko{d.cislo}/kgml", online)
            if txt:
                try:
                    ref = ET.fromstring(txt)
                    d.zdroj_dat = f"referenční mapa ko{d.cislo} (KEGG)"
                except ET.ParseError:
                    ref = koren
            else:
                log("  → nedostupná, použiji jen reakce z KGML organismu.")
        r_entries, r_cpd, r_sour = _parse_kgml(ref)
        for c, xy in r_sour.items():
            souradnice.setdefault(c, xy)
        for r in ref.findall("reaction"):
            rn = [x.replace("rn:", "") for x in r.get("name").split()]
            typ = "<=>" if r.get("type") == "reversible" else "->"
            subs = [r_cpd.get(s.get("id")) for s in r.findall("substrate")]
            prod = [r_cpd.get(p.get("id")) for p in r.findall("product")]
            for a in subs:
                for b in prod:
                    if a and b and a != b:
                        d.reakce.append(dict(rn=rn, a=a, b=b, ec=[], typ=typ))
        # EC čísla z KEGG
        vsechny_rn = sorted({x for r in d.reakce for x in r["rn"]})
        log(f"Stahuji EC čísla pro {len(vsechny_rn)} reakcí…")
        for rid, zaz in kegg_zaznamy(vsechny_rn, online).items():
            ec_rn[rid] = pole_zaznamu(zaz, "ENZYME").split()
        for r in d.reakce:
            for x in r["rn"]:
                for ec in ec_rn.get(x, []):
                    if ec not in r["ec"]:
                        r["ec"].append(ec)
            # bez spojení s KEGG neznáme EC -> ukázat aspoň číslo reakce (ne „neznámý enzym“)
            if not r["ec"] and not any(x in ec_rn for x in r["rn"]):
                r["ec"] = list(r["rn"])

    if not d.reakce:
        typy = defaultdict(int)
        for rel in koren.findall("relation"):
            typy[rel.get("type")] += 1
        if typy.get("PPrel") or typy.get("GErel") or typy.get("PCrel"):
            duvod = ("Je to signální / regulační dráha – v KGML jsou jen vztahy mezi proteiny "
                     f"({', '.join(f'{k} {v}×' for k, v in typy.items())}), ne reakce mezi látkami.")
        elif d.cislo.startswith(("011", "012")):
            duvod = "Je to přehledová mapa KEGG, která nemá v KGML popsané jednotlivé reakce."
        else:
            duvod = "KGML neobsahuje žádné reakce mezi metabolity (a referenční mapu se nepodařilo stáhnout)."
        raise ValueError(f"Dráhu {d.org}{d.cislo} „{d.titul}“ nelze nakreslit.\n\n{duvod}\n\n"
                         "Program kreslí metabolické dráhy (mapy 00xxx), např. 00010, 00020, 00360.")
    for r in d.reakce:
        r["geny"] = [g for x in r["rn"] for g in geny.get(x, [])]
        if r["typ"] == "<=>" and any(x in NEVRATNE for x in r["rn"]):
            r["typ"] = "->"

    # metabolity
    cids = []
    for r in d.reakce:
        for c in (r["a"], r["b"]):
            if c not in cids:
                cids.append(c)
    glykany = [c for c in cids if c.startswith("G")]
    chybi_nazev = [c for c in cids if c not in CESKE_NAZVY]
    zaznamy = {}
    if chybi_nazev or glykany:
        log(f"Stahuji názvy {len(set(chybi_nazev + glykany))} metabolitů…")
        zaznamy = kegg_zaznamy(chybi_nazev + glykany, online)
    chybi_strukt = [c for c in cids if c not in SMILES and not c.startswith("G")]
    if chybi_strukt:
        log(f"Stahuji struktury {len(chybi_strukt)} metabolitů…")
    for c in cids:
        zaz = zaznamy.get(c, "")
        nazev = CESKE_NAZVY.get(c)
        nazev_en = pole_zaznamu(zaz, "NAME").split("\n")[0].strip().rstrip(";")
        if not nazev:
            n = nazev_en
            if not n and c.startswith("G"):                 # glykan bez názvu -> složení
                n = pole_zaznamu(zaz, "COMPOSITION").strip()
            nazev = n or c
        stejne = re.search(r"Same as:\s*(C\d{5})", zaz)
        if c.startswith("G"):
            zdroj = ("KCF\n" + zaz) if "NODE" in zaz else None   # strom cukrů (SNFG)
        else:
            zdroj = SMILES.get(c) or kegg_get(f"get/{c}/mol", online)
        d.metabolity[c] = dict(nazev=nazev, zdroj=zdroj, poradi=souradnice.get(c, (9e9, 9e9)),
                               stejne_jako=stejne.group(1) if stejne else None,
                               nazev_en=nazev_en, role=role_latky(c, nazev_en or nazev))
    d.online = online

    # odkazy na navazující dráhy
    if d.cislo in VESTAVENE_ODKAZY:
        for c, lst in VESTAVENE_ODKAZY[d.cislo].items():
            for m, smer in lst:
                d.odkazy.append((c, m, CESKE_MAPY.get(m, m), smer))
    else:
        videno = set()
        for koren_ in ([koren] + ([ref] if ref is not None and ref is not koren else [])):
            ent, cp, _ = _parse_kgml(koren_)
            for rel in koren_.findall("relation"):
                if rel.get("type") != "maplink":
                    continue
                for st in rel.findall("subtype"):
                    c = cp.get(st.get("value"))
                    if c not in d.metabolity:
                        continue
                    for eid in (rel.get("entry1"), rel.get("entry2")):
                        e = ent.get(eid)
                        if e is None or e.get("type") != "map":
                            continue
                        m = re.sub(r"\D", "", e.get("name"))[-5:]
                        if m == d.cislo or (c, m) in videno:
                            continue
                        videno.add((c, m))
                        nazev = CESKE_MAPY.get(m) or e.find("graphics").get("name")
                        d.odkazy.append((c, m, nazev, "->"))

    # výchozí metabolit
    if d.cislo in VESTAVENY_START and VESTAVENY_START[d.cislo] in d.metabolity:
        d.start = VESTAVENY_START[d.cislo]
    else:
        # výchozí = látka, ze které je po směru reakcí dosažitelné nejvíc dalších
        nasl = defaultdict(set)
        bez_nosicu = [r for r in d.reakce
                      if d.metabolity.get(r["a"], {}).get("role", "normal") == "normal"
                      and d.metabolity.get(r["b"], {}).get("role", "normal") == "normal"]
        for r in bez_nosicu:
            nasl[r["a"]].add(r["b"])
            if r["typ"] == "<=>":
                nasl[r["b"]].add(r["a"])

        def dosah(c):
            vid, q = {c}, [c]
            while q:
                for v in nasl[q.pop()]:
                    if v not in vid:
                        vid.add(v)
                        q.append(v)
            return len(vid)
        vstup = defaultdict(int)
        for r in bez_nosicu:
            vstup[r["b"]] += 1
            if r["typ"] == "<=>":
                vstup[r["a"]] += 1
        s_odkazem = {c for c, *_ in d.odkazy}
        # přednost: látka, do které nic nevede a která je napojená na jinou dráhu
        bezne = [c for c, m in d.metabolity.items() if m.get("role", "normal") == "normal"]
        d.start = (max(bezne or list(d.metabolity), key=lambda c: (vstup[c] == 0 and c in s_odkazem,
                                                    vstup[c] == 0, dosah(c),
                                                    -d.metabolity[c]["poradi"][0]))
                   if d.metabolity else None)
    log(f"Načteno: {len(d.metabolity)} metabolitů, {len(d.reakce)} reakčních kroků, "
        f"{len(d.odkazy)} odkazů na další dráhy ({d.zdroj_dat}).")
    if len(d.metabolity) > VELKA_MAPA:
        log(f"Pozor: velká mapa ({len(d.metabolity)} metabolitů) – doporučuji vybrat jen část.")
    return d


# =============================================================================
# 4) VÝBĚR -> GRAF PRO KRESLENÍ (včetně šipek přes vynechané kroky)
# =============================================================================
PODPOROVANE_DRAHY = """Které dráhy program umí

✔ Metabolické dráhy KEGG – mapy 00xxx (např. 00010 glykolýza, 00020 citrátový cyklus,
   00230 puriny, 00360 fenylalanin, 00380 tryptofan, 00600 sfingolipidy …),
   pro libovolný organismus (gga, hsa, mmu, eco …) i referenční mapy (map, ko).
✔ Glykanové dráhy (např. 00510, 00604) – glykany se kreslí symboly SNFG,
   volitelně jako atomové vzorce.
✔ Cyklické dráhy – kruhové rozložení (cyklus po kružnici, větve ven).
✔ Látky s koenzymem A se kreslí se zkratkou „CoA“ (lze vypnout).

✘ Signální a regulační dráhy (04xxx), dráhy nemocí (05xxx) – v KGML nejsou reakce
   mezi látkami, jen vztahy mezi proteiny.
✘ Přehledové mapy (01100, 01200 …) – příliš velké, bez popsaných reakcí.

Na co si dát pozor
• Cyklické dráhy (citrátový cyklus, močovinový cyklus): zvol Rozložení → kruhové.
  Program vezme nejdelší cyklus ve vybrané části dráhy a nakreslí ho po kružnici.
• Přenašeče (ThPP, lipoyllysin, NAD …) a meziprodukty vázané na enzymový komplex
  (PDH, OGDH) jsou po načtení vypnuté – komplex je pak jedna šipka. Lze je zapnout.
• Látky bez struktury v KEGG (protein, akceptor, polymery) mají jen rámeček s názvem.
• KGML někdy uvádí nevratnou reakci obráceně – karta Reakce → „Otočit směr“.
• U map nad 120 látek doporučuji vypnout části, které nepotřebuješ.
• Data z KEGG se ukládají do složky kegg_cache/ (podruhé to jde i offline).
  KEGG je pro akademické použití zdarma – viz https://www.kegg.jp/kegg/legal.html
"""

# Tyto globální proměnné čte kreslicí jádro (sekce 7).
METABOLITY = {}
MAPY = {}
REAKCE = []
VYCHOZI = None
NADPIS = ""
PODNADPIS = ""
ORGANISMUS = ""
ZVYRAZNIT = False
POSLEDNI_SVG = ""  # text posledního vytvořeného SVG (pro náhled)
MEZERY = {}       # id rámečku vynechaných kroků -> seznam skrytých metabolitů


def priprav_graf(d, vybrane_met, vybrane_rx, start=None, otocene=(), odkazy=True,
                 geny=True, zkratky=True):
    """Z vybraných metabolitů/reakcí sestaví seznam uzlů a REAKCE pro kreslení."""
    global METABOLITY, MAPY, VYCHOZI
    start_rozl = start if start in d.metabolity else d.start   # i když je skrytý
    start = start_rozl if start_rozl in vybrane_met else None   # zvýraznění jen viditelného
    rx = []
    def role_(c):
        return d.metabolity.get(c, {}).get("role", "normal")
    for i in sorted(vybrane_rx):
        r = dict(d.reakce[i])
        if i in otocene:
            r["a"], r["b"] = r["b"], r["a"]
        if any(c not in vybrane_met and role_(c) == "prenasec" for c in (r["a"], r["b"])):
            continue                        # vypnutý přenašeč (ThPP, lipoyllysin …) se vynechá úplně
        rx.append(r)

    # pořadí pro rozložení: vzdálenost od výchozího metabolitu přes VŠECHNY vybrané reakce
    sous = defaultdict(set)
    for r in rx:
        sous[r["a"]].add(r["b"])
        sous[r["b"]].add(r["a"])
    poradi = {c: d.metabolity[c]["poradi"] for c in d.metabolity}
    dist = {}
    zdroje = [start_rozl] if start_rozl in sous else []
    zdroje += sorted(sous, key=lambda c: poradi[c])
    for z in zdroje:                       # i nesouvislé části grafu
        if z in dist:
            continue
        dist[z] = 0
        q = deque([z])
        while q:
            u = q.popleft()
            for v in sorted(sous[u], key=lambda c: poradi[c]):
                if v not in dist:
                    dist[v] = dist[u] + 1
                    q.append(v)
    klic = lambda c: (dist.get(c, 99), poradi[c])

    # přímé hrany (oba konce vybrané), sloučené podle dvojice
    hrany = {}
    for r in rx:
        a, b = r["a"], r["b"]
        if a not in vybrane_met or b not in vybrane_met:
            continue
        k = frozenset((a, b))
        h = hrany.setdefault(k, dict(smery=set(), ec=[], geny=[], vratna=False, plna=False, rn=set()))
        if r["typ"] == "->+":
            h["plna"] = True
        if r["typ"] == "<=>":
            h["vratna"] = True
        h["smery"].add((a, b))
        h["rn"].update(r["rn"])
        for x in r["ec"]:
            if x not in h["ec"]:
                h["ec"].append(x)
        for x in r["geny"]:
            if x not in h["geny"]:
                h["geny"].append(x)

    REAKCE.clear()
    for k, h in hrany.items():
        a, b = sorted(k, key=klic)
        vratna = h["vratna"] or len(h["smery"]) > 1
        zpet = (not vratna) and (b, a) in h["smery"]
        carkovana = not h["ec"] and not h["plna"]
        REAKCE.append((a, b, " · ".join(h["geny"]) if geny else "", " · ".join(h["ec"]),
                       vratna, zpet, bool(h["geny"]), "", carkovana, 0))

    # vynechané metabolity: každá souvislá skupina skrytých látek, která leží MEZI
    # vybranými látkami, se nahradí jedním rámečkem „⋯ vynechány N metabolity“
    # a tečkovanými šipkami. Skrytá slepá větev se nekreslí vůbec.
    MEZERY.clear()
    mezery_uzly = []
    if True:
        skryte = {c for r in rx for c in (r["a"], r["b"]) if c not in vybrane_met}
        sous_s = defaultdict(set)
        for r in rx:
            if r["a"] in skryte and r["b"] in skryte:
                sous_s[r["a"]].add(r["b"])
                sous_s[r["b"]].add(r["a"])
        videno = set()
        for c0 in sorted(skryte, key=klic):
            if c0 in videno:
                continue
            komp, q = [], deque([c0])
            videno.add(c0)
            while q:
                u = q.popleft()
                komp.append(u)
                for v in sous_s[u]:
                    if v not in videno:
                        videno.add(v)
                        q.append(v)
            kset = set(komp)
            k_gap = min(klic(c) for c in komp)      # rámeček stojí na místě 1. skryté látky
            if all(role_(c) == "vazany" for c in komp):
                # řetězec meziproduktů vázaných na enzymový komplex -> jedna běžná šipka
                pred_, nasl_, ecs, gs = set(), set(), [], []
                for r in rx:
                    a, b = r["a"], r["b"]
                    if a not in kset and b not in kset:
                        continue
                    ecs += [x for x in r["ec"] if x not in ecs]
                    gs += [x for x in r["geny"] if x not in gs]
                    # směr podle zápisu reakce v KEGG (substrát -> produkt)
                    if a in vybrane_met and b in kset:
                        pred_.add(a)
                    if b in vybrane_met and a in kset:
                        nasl_.add(b)
                for u in sorted(pred_, key=klic):
                    for v in sorted(nasl_, key=klic):
                        if u == v or frozenset((u, v)) in hrany:
                            continue
                        hrany[frozenset((u, v))] = None
                        a_, b_ = sorted((u, v), key=klic)
                        REAKCE.append((a_, b_, " · ".join(gs) if geny else "", " · ".join(ecs),
                                       False, a_ != u, bool(gs), "", False, 0))
                continue
            if not zkratky:
                continue
            pred_, nasl_ = set(), set()
            for r in rx:
                for s_, w in ((r["a"], r["b"]), (r["b"], r["a"])):
                    if s_ in vybrane_met and w in kset:
                        (pred_ if klic(s_) < k_gap else nasl_).add(s_)
            if len(pred_ | nasl_) < 2:
                continue                    # slepá větev – vynechat úplně
            gid = f"gap:{len(MEZERY)}"
            MEZERY[gid] = sorted(komp, key=klic)
            mezery_uzly.append((gid, k_gap))
            for u in sorted(pred_, key=klic):
                REAKCE.append((u, gid, "", "", False, False, False, "", False, -1))
            for v in sorted(nasl_, key=klic):
                REAKCE.append((gid, v, "", "", False, False, False, "", False, -1))

    klic_gap = {g: k for g, k in mezery_uzly}
    klic2 = lambda c: klic_gap[c] if c in klic_gap else klic(c)
    uzly = sorted(list(vybrane_met) + list(klic_gap), key=klic2)
    MAPY = {}
    if odkazy:
        for c, m, nazev, smer in d.odkazy:
            if c not in vybrane_met:
                continue
            mid = f"map:{m}:{c}"
            if mid in uzly:
                continue
            MAPY[m] = nazev
            uzly.append(mid)
            if smer == "->":
                REAKCE.append((c, mid, "", "", False, False, False, "", True, 0))
            else:
                REAKCE.append((mid, c, "", "", False, False, False, "", True, 0))
    REAKCE.sort(key=lambda r: (klic2(r[0]) if not r[0].startswith("map:") else (-1, (0, 0)),
                               klic2(r[1]) if not r[1].startswith("map:") else (99, (0, 0))))
    METABOLITY = {c: (d.metabolity[c]["nazev"], d.metabolity[c]["zdroj"], None) for c in d.metabolity}
    VYCHOZI = start
    return uzly


def vytvor_svg(d, vybrane_met, vybrane_rx, vystup, start=None, otocene=(), popisky="pod",
               odkazy=True, geny=True, zvyraznit=False, zkratky=True, nadpis=None,
               glykany_atomy=False, zkratit_coa=True, delka_vazby=None, velikost_atomu=None,
               rozlozeni="radky", vzhled=None):
    """Hlavní funkce: z výběru vytvoří SVG soubor. Vrací (počet uzlů, počet křížení)."""
    global NADPIS, PODNADPIS, ORGANISMUS, ZVYRAZNIT, ZKRACOVAT_COA, DELKA_VAZBY, FONT_ATOMU
    pouzij_vzhled(vzhled)
    ZVYRAZNIT = zvyraznit
    ZKRACOVAT_COA = zkratit_coa
    if delka_vazby:
        DELKA_VAZBY = float(delka_vazby)
    if velikost_atomu:
        FONT_ATOMU = float(velikost_atomu)
    ORGANISMUS = d.org
    uzly = priprav_graf(d, set(vybrane_met), set(vybrane_rx), start, set(otocene),
                        odkazy, geny, zkratky)
    if not uzly:
        raise ValueError("Není vybrán žádný metabolit.")
    if glykany_atomy:               # glykany jako plné atomové vzorce (přes „Same as: C…“)
        for c, m in d.metabolity.items():
            if c.startswith("G") and m.get("stejne_jako"):
                mol = SMILES.get(m["stejne_jako"]) or kegg_get(f'get/{m["stejne_jako"]}/mol', d.online)
                if mol:
                    METABOLITY[c] = (m["nazev"], mol, None)
    POUZITE_CUKRY.clear()
    NADPIS = (nadpis or "").strip() or CESKE_MAPY.get(d.cislo, d.titul) or "Metabolická dráha"
    n_met = len(vybrane_met)
    n_rx = sum(1 for r in REAKCE if not r[9] and ":" not in r[0] + r[1])
    n_zk = len(MEZERY)
    PODNADPIS = (f"KEGG {d.org}{d.cislo} · {n_met} metabolitů, {n_rx} přeměn"
                 + (f", {n_zk}× vynechané kroky" if n_zk else "")
                 + " · struktury RDKit · celý obrázek vektorově")
    global POSLEDNI_SVG
    svg, krizeni = (sestav_kruh if rozlozeni == "kruh" else sestav)(uzly, popisky)
    POSLEDNI_SVG = svg
    if vystup:                      # vystup=None -> jen náhled, nic se neukládá
        with open(vystup, "w", encoding="utf-8") as f:
            f.write(svg)
    return len(uzly), krizeni




# =============================================================================
# 5) VZHLED
# =============================================================================
DELKA_VAZBY = 20.0        # px – stejné měřítko pro všechny struktury
FONT_ATOMU = 0.50         # relativní velikost atomových symbolů (RDKit default 0.6)
TLOUSTKA_VAZBY = 1.3
FONT = "DejaVu Sans"      # dodává se s matplotlibem -> reprodukovatelné
VEL_NAZEV = 12.5          # velikost názvu metabolitu
VEL_GEN = 10.5
VEL_EC = 8.5
VEL_MAPA = 10.5           # text v rámečku navazující dráhy
VEL_NADPIS = 22
VEL_PODNADPIS = 10

MEZERA_SVISLE = 34        # min. svislá mezera mezi uzly ve vrstvě
MIN_MEZERA_VRSTEV = 120   # min. vodorovná mezera mezi vrstvami
OKRAJ_POPISKU = 22        # rezerva kolem popisku enzymu v mezeře
OKRAJ = 40                # okraj obrázku
ODSTUP_SIPKY = 7          # vzdálenost hrotu šipky od struktury
HROT_DELKA = 9
HROT_SIRKA = 3.6
ROZESTUP_VRATNE = 2.6     # polovina vzdálenosti dvou čar u vratné reakce
PAD_KOLEJ = 18            # odstup prvni svisle koleje od sloupce
ROZESTUP_KOLEJI = 11      # vzdálenost sousedních svislých kolejí
LIMIT_RAZENI_S = 20      # max. čas na minimalizaci křížení (s)
VELKA_MAPA = 120          # od kolika metabolitů varovat
POLOMER_ROHU = 4          # malé zaoblení pravoúhlých ohybů (0 = ostré rohy)
BARVA_SIPKY = "#3a3a3a"     # všechny šipky stejně
BARVA_GEN = "#222222"
BARVA_EC = "#555555"
BARVA_NAZEV = "#111111"
BARVA_ORG = "#1b7a3a"     # jen při volbě „zvýraznit enzymy organismu“
BARVA_ZKRATKA = "#7b3fa0" # šipka přes vynechané metabolity
BARVA_POZADI = "#ffffff"
BARVA_NADPIS = "#111111"
BARVA_POPIS = "#555555"   # podnadpis, poznámky, legenda
BARVA_VAZEB = "#000000"   # vazby a uhlíky ve strukturách
BAREVNE_ATOMY = True      # O červeně, N modře … (jinak vše barvou vazeb)
BARVA_VYCHOZI = "#fff4dc"
BARVA_VYCHOZI_OBRYS = "#e0a43a"
BARVA_MAPY = "#f4f7fb"
BARVA_MAPY_OBRYS = "#4a6fa5"
BARVA_MAPY_TEXT = "#1f3a60"
BARVA_MEZERY = "#faf6fd"  # výplň rámečku vynechaných kroků
TVAR_HROTU = "plny"       # plny | trojuhelnik | otevreny
TLOUSTKA_SIPEK = 1.3
ZAOBLENI_RAMECKU = 10


# --- Předvolby vzhledu a jejich použití ---------------------------------------
VYCHOZI_VZHLED = dict(
    font="DejaVu Sans", barva_pozadi="#ffffff", barva_nadpisu="#111111", barva_popisu="#555555",
    barva_sipek="#3a3a3a", barva_nazvu="#111111", barva_genu="#222222", barva_ec="#555555",
    barva_vazeb="#000000", barevne_atomy=True, barva_vychozi="#fff4dc", barva_vychozi_obrys="#e0a43a",
    barva_mapy="#f4f7fb", barva_mapy_obrys="#4a6fa5", barva_mapy_text="#1f3a60",
    barva_zkratky="#7b3fa0", barva_mezery="#faf6fd", barva_org="#1b7a3a",
    hrot="plny", tloustka=1.3, velikost_hrotu=9.0, polomer_rohu=4.0, zaobleni_ramecku=10.0,
    delka_vazby=20.0, velikost_atomu=0.5, pismo_nazvu=12.5, pismo_enzymu=10.5, rozestup=1.0,
)
PREDVOLBY = {
    "Klasický": {},
    "Barevný": dict(barva_sipek="#2f5d8a", barva_genu="#1f4e8c", barva_ec="#5b6f86",
                    barva_vychozi="#e6f2ff", barva_vychozi_obrys="#3b82c4",
                    barva_mapy="#fff4e0", barva_mapy_obrys="#e08a1e", barva_mapy_text="#8a4b00",
                    barva_zkratky="#c2185b", barva_mezery="#fdf0f5", hrot="trojuhelnik"),
    "Černobílý (tisk)": dict(barva_sipek="#000000", barva_genu="#000000", barva_ec="#333333",
                             barva_popisu="#333333", barevne_atomy=False, barva_vychozi="#f0f0f0",
                             barva_vychozi_obrys="#000000", barva_mapy="#ffffff", barva_mapy_obrys="#000000",
                             barva_mapy_text="#000000", barva_zkratky="#555555", barva_mezery="#ffffff",
                             barva_org="#000000", tloustka=1.1),
    "Tmavý": dict(barva_pozadi="#1e2026", barva_nadpisu="#f2f2f2", barva_popisu="#a9adb5",
                  barva_sipek="#c9ccd2", barva_nazvu="#f0f0f0", barva_genu="#e6e6e6", barva_ec="#a9adb5",
                  barva_vazeb="#e8e8e8", barva_vychozi="#3a3324", barva_vychozi_obrys="#e0a43a",
                  barva_mapy="#2a3140", barva_mapy_obrys="#7d9cd1", barva_mapy_text="#c9d7f0",
                  barva_zkratky="#c792ea", barva_mezery="#2e2536", barva_org="#6fd38f"),
}


def pouzij_vzhled(v=None):
    """Nastaví globální proměnné vzhledu ze slovníku (chybějící klíče = výchozí)."""
    global FONT, BARVA_POZADI, BARVA_NADPIS, BARVA_POPIS, BARVA_SIPKY, BARVA_NAZEV, BARVA_GEN, BARVA_EC
    global BARVA_VAZEB, BAREVNE_ATOMY, BARVA_VYCHOZI, BARVA_VYCHOZI_OBRYS, BARVA_MAPY, BARVA_MAPY_OBRYS
    global BARVA_MAPY_TEXT, BARVA_ZKRATKA, BARVA_MEZERY, BARVA_ORG, TVAR_HROTU, TLOUSTKA_SIPEK
    global HROT_DELKA, HROT_SIRKA, POLOMER_ROHU, ZAOBLENI_RAMECKU, DELKA_VAZBY, FONT_ATOMU
    global VEL_NAZEV, VEL_GEN, VEL_EC, VEL_MAPA, MEZERA_SVISLE, MIN_MEZERA_VRSTEV, ROZESTUP_KOLEJI
    z = dict(VYCHOZI_VZHLED)
    z.update({k: val for k, val in (v or {}).items() if k in VYCHOZI_VZHLED})
    FONT = z["font"]
    BARVA_POZADI, BARVA_NADPIS, BARVA_POPIS = z["barva_pozadi"], z["barva_nadpisu"], z["barva_popisu"]
    BARVA_SIPKY, BARVA_NAZEV, BARVA_GEN, BARVA_EC = z["barva_sipek"], z["barva_nazvu"], z["barva_genu"], z["barva_ec"]
    BARVA_VAZEB, BAREVNE_ATOMY = z["barva_vazeb"], bool(z["barevne_atomy"])
    BARVA_VYCHOZI, BARVA_VYCHOZI_OBRYS = z["barva_vychozi"], z["barva_vychozi_obrys"]
    BARVA_MAPY, BARVA_MAPY_OBRYS, BARVA_MAPY_TEXT = z["barva_mapy"], z["barva_mapy_obrys"], z["barva_mapy_text"]
    BARVA_ZKRATKA, BARVA_MEZERY, BARVA_ORG = z["barva_zkratky"], z["barva_mezery"], z["barva_org"]
    TVAR_HROTU = z["hrot"]
    TLOUSTKA_SIPEK = float(z["tloustka"])
    HROT_DELKA = float(z["velikost_hrotu"])
    HROT_SIRKA = HROT_DELKA * 0.4
    POLOMER_ROHU = float(z["polomer_rohu"])
    ZAOBLENI_RAMECKU = float(z["zaobleni_ramecku"])
    DELKA_VAZBY = float(z["delka_vazby"])
    FONT_ATOMU = float(z["velikost_atomu"])
    VEL_NAZEV = float(z["pismo_nazvu"])
    VEL_GEN = float(z["pismo_enzymu"])
    VEL_EC = VEL_GEN * 0.81
    VEL_MAPA = VEL_GEN
    MEZERA_SVISLE = 34 * float(z["rozestup"])
    MIN_MEZERA_VRSTEV = 120 * float(z["rozestup"])
    ROZESTUP_KOLEJI = max(7.0, 11 * min(1.3, float(z["rozestup"])), TLOUSTKA_SIPEK * 2 + ROZESTUP_VRATNE * 2 + 3)


def nacti_vzhled_ze_souboru(cesta):
    import json
    with open(cesta, encoding="utf-8") as f:
        return {k: v for k, v in json.load(f).items() if k in VYCHOZI_VZHLED}


def tmave_pozadi():
    from matplotlib.colors import to_rgb
    r, g, b = to_rgb(BARVA_POZADI)
    return 0.299 * r + 0.587 * g + 0.114 * b < 0.5


# =============================================================================
# 5b) TEXT -> KŘIVKY
# =============================================================================
def _fmt(v):
    s = f"{v:.2f}".rstrip("0").rstrip(".")
    return "0" if s in ("-0", "") else s


class Text:
    """Text převedený na obrysové křivky. Rozměry v px, počátek = levý horní roh."""

    def __init__(self, s, size, weight="normal", style="normal"):
        prop = FontProperties(family=FONT, weight=weight, style=style)
        self.tp = TextPath((0, 0), s, size=size, prop=prop)
        # stabilní výška řádku nezávislá na konkrétních písmenech
        ref = TextPath((0, 0), "ÁHgjy", size=size, prop=prop).get_extents()
        ext = self.tp.get_extents()
        self.ascent, self.descent = ref.y1, -ref.y0
        self.w = ext.x1 if s else 0.0
        self.h = self.ascent + self.descent

    def path(self, x, y, fill, halo=None):
        """SVG <path>; (x, y) = levý horní roh textového boxu."""
        base = y + self.ascent
        d = []
        for verts, code in self.tp.iter_segments(curves=True, simplify=False):
            pts = [(_fmt(x + verts[i]), _fmt(base - verts[i + 1]))
                   for i in range(0, len(verts), 2)]
            if code == 1:
                d.append("M%s,%s" % pts[0])
            elif code == 2:
                d.append("L%s,%s" % pts[0])
            elif code == 3:
                d.append("Q%s,%s %s,%s" % (*pts[0], *pts[1]))
            elif code == 4:
                d.append("C%s,%s %s,%s %s,%s" % (*pts[0], *pts[1], *pts[2]))
            elif code == 79:
                d.append("Z")
        d = "".join(d)
        out = ""
        if halo:  # bílý obrys pod textem – čitelnost přes čáry
            out += (f'<path d="{d}" fill="{halo}" stroke="{halo}" '
                    f'stroke-width="3.2" stroke-linejoin="round"/>')
        return out + f'<path d="{d}" fill="{fill}"/>'


# =============================================================================
# 5c) STRUKTURY (RDKit -> SVG s cestami)
# =============================================================================
def _orientuj(mol):
    """Natočí molekulu tak, aby byla protáhlá vodorovně a fenyl byl vlevo
    (všechny struktury pak ‚tečou‘ stejným směrem jako dráha)."""
    conf = mol.GetConformer()
    pts = [conf.GetAtomPosition(i) for i in range(mol.GetNumAtoms())]
    n = len(pts)
    cx = sum(p.x for p in pts) / n
    cy = sum(p.y for p in pts) / n
    sxx = sum((p.x - cx) ** 2 for p in pts)
    syy = sum((p.y - cy) ** 2 for p in pts)
    sxy = sum((p.x - cx) * (p.y - cy) for p in pts)
    ang = 0.5 * math.atan2(2 * sxy, sxx - syy)   # hlavní osa
    c, s = math.cos(-ang), math.sin(-ang)
    new = []
    for p in pts:
        x, y = p.x - cx, p.y - cy
        new.append([x * c - y * s, x * s + y * c])
    # fenylový kruh (aromatické atomy) doleva
    arom = [i for i, a in enumerate(mol.GetAtoms()) if a.GetIsAromatic()]
    if arom:
        ax = sum(new[i][0] for i in arom) / len(arom)
        if ax > 0:
            for q in new:
                q[0] = -q[0]
    for i, (x, y) in enumerate(new):
        conf.SetAtomPosition(i, (x, y, 0.0))


# --- glykany: symbolická nomenklatura SNFG (Varki et al., Glycobiology 2015) ---
SNFG = {
    "Glc": ("kruh", "#0072BC"), "Gal": ("kruh", "#FFD400"), "Man": ("kruh", "#00A651"),
    "GlcNAc": ("ctverec", "#0072BC"), "GalNAc": ("ctverec", "#FFD400"), "ManNAc": ("ctverec", "#00A651"),
    "GlcN": ("ctverec_pul", "#0072BC"), "GalN": ("ctverec_pul", "#FFD400"),
    "GlcA": ("kosoctverec_pul", "#0072BC"), "GalA": ("kosoctverec_pul", "#FFD400"),
    "IdoA": ("kosoctverec_pul", "#A17A4D"),
    "Neu5Ac": ("kosoctverec", "#A54399"), "Neu5Gc": ("kosoctverec", "#8FCCE9"),
    "Kdn": ("kosoctverec", "#00A651"), "Fuc": ("trojuhelnik", "#ED1C24"),
    "Rha": ("trojuhelnik", "#00A651"), "Xyl": ("hvezda", "#F47920"), "Ara": ("hvezda", "#00A651"),
}
POUZITE_CUKRY = set()
GLYKAN_MERITKO = 4.2      # px na jednotku souřadnic KCF
GLYKAN_SYMBOL = 7.0       # poloměr symbolu cukru


def parsuj_kcf(txt):
    """NODE/EDGE ze záznamu KEGG GLYCAN -> ({id: (zbytek, x, y)}, [(dítě, rodič, popisek)])."""
    uzly_, hrany_, sekce = {}, [], None
    for radek in txt.splitlines():
        if radek[:1].strip():
            sekce = radek.split()[0]
            radek = radek[12:]
        if sekce == "NODE":
            c = radek.split()
            if len(c) >= 4 and c[0].isdigit():
                uzly_[int(c[0])] = (c[1], float(c[2]), float(c[3]))
        elif sekce == "EDGE":
            c = radek.split()
            if len(c) >= 3 and c[0].isdigit():
                dite, _, a = c[1].partition(":")
                rodic, _, b = c[2].partition(":")
                a = a.replace("a", "α").replace("b", "β")
                popis = f"{a}-{b}" if (a or b) else ""
                hrany_.append((int(re.sub(r"\D", "", dite)), int(re.sub(r"\D", "", rodic)), popis))
    return uzly_, hrany_


def symbol_cukru(tvar, barva, x, y, r):
    """SVG cesty jednoho symbolu SNFG se středem (x, y)."""
    f = _fmt
    obrys = f'stroke="{BARVA_VAZEB}" stroke-width="0.9" stroke-linejoin="round"'
    if tvar == "kruh":
        d = f"M{f(x - r)},{f(y)} a{f(r)},{f(r)} 0 1,0 {f(2 * r)},0 a{f(r)},{f(r)} 0 1,0 {f(-2 * r)},0 Z"
        return f'<path d="{d}" fill="{barva}" {obrys}/>'
    if tvar in ("ctverec", "ctverec_pul"):
        q = r * 0.9
        d = f"M{f(x - q)},{f(y - q)} H{f(x + q)} V{f(y + q)} H{f(x - q)} Z"
        if tvar == "ctverec":
            return f'<path d="{d}" fill="{barva}" {obrys}/>'
        t = f"M{f(x - q)},{f(y - q)} H{f(x + q)} L{f(x - q)},{f(y + q)} Z"
        return f'<path d="{d}" fill="#ffffff" {obrys}/><path d="{t}" fill="{barva}"/><path d="{d}" fill="none" {obrys}/>'
    if tvar in ("kosoctverec", "kosoctverec_pul"):
        q = r * 1.15
        d = f"M{f(x)},{f(y - q)} L{f(x + q)},{f(y)} L{f(x)},{f(y + q)} L{f(x - q)},{f(y)} Z"
        if tvar == "kosoctverec":
            return f'<path d="{d}" fill="{barva}" {obrys}/>'
        t = f"M{f(x)},{f(y - q)} L{f(x + q)},{f(y)} L{f(x - q)},{f(y)} Z"
        return f'<path d="{d}" fill="#ffffff" {obrys}/><path d="{t}" fill="{barva}"/><path d="{d}" fill="none" {obrys}/>'
    if tvar == "trojuhelnik":
        q = r * 1.1
        d = f"M{f(x)},{f(y - q)} L{f(x + q)},{f(y + q * 0.8)} L{f(x - q)},{f(y + q * 0.8)} Z"
        return f'<path d="{d}" fill="{barva}" {obrys}/>'
    if tvar == "hvezda":
        body_ = []
        for k in range(10):
            rr = r * 1.15 if k % 2 == 0 else r * 0.5
            u = -math.pi / 2 + k * math.pi / 5
            body_.append(f"{f(x + rr * math.cos(u))},{f(y + rr * math.sin(u))}")
        return f'<path d="M{" L".join(body_)} Z" fill="{barva}" {obrys}/>'
    return ""


def nakresli_glykan(kcf):
    """Glykan jako symboly SNFG (vše vektorově). Vrací (šířka, výška, SVG) nebo None."""
    uzly_, hrany_ = parsuj_kcf(kcf)
    if not uzly_:
        return None
    S, r, pad = GLYKAN_MERITKO, GLYKAN_SYMBOL, 10
    xs = [u[1] for u in uzly_.values()]
    ys = [u[2] for u in uzly_.values()]
    texty = {}
    for i, (zb, x, y) in uzly_.items():
        if zb not in SNFG:
            texty[i] = Text(zb, VEL_EC + 0.5, weight="bold")
    extra = max([t.w / 2 for t in texty.values()], default=0)
    x0, y0 = min(xs), min(ys)
    pos = {i: ((x - x0) * S + pad + extra, (y - y0) * S + pad + 4) for i, (_, x, y) in uzly_.items()}
    w = (max(xs) - x0) * S + 2 * pad + 2 * extra
    h = (max(ys) - y0) * S + 2 * pad + 8
    out = []
    for dite, rodic, popis in hrany_:          # vazby
        if dite in pos and rodic in pos:
            (ax, ay), (bx, by) = pos[dite], pos[rodic]
            out.append(f'<path d="M{_fmt(ax)},{_fmt(ay)} L{_fmt(bx)},{_fmt(by)}" '
                       f'stroke="{BARVA_VAZEB}" stroke-width="1"/>')
    for dite, rodic, popis in hrany_:          # popisky vazeb (α2-3 …)
        if popis and dite in pos and rodic in pos:
            (ax, ay), (bx, by) = pos[dite], pos[rodic]
            t = Text(popis, 6.3)
            mx, my = ax + (bx - ax) * 0.42, ay + (by - ay) * 0.42
            out.append(t.path(mx - t.w / 2, my - t.h - 2, BARVA_EC, halo=BARVA_POZADI))
    for i, (zb, _, _) in uzly_.items():        # symboly / text
        x, y = pos[i]
        if zb in SNFG:
            POUZITE_CUKRY.add(zb)
            out.append(symbol_cukru(*SNFG[zb], x, y, r))
        else:
            t = texty[i]
            out.append(t.path(x - t.w / 2, y - t.h / 2, BARVA_NAZEV, halo=BARVA_POZADI))
    return w, h, "".join(out)


# Velké kofaktory se ve stažených strukturách nahradí zkratkou (vazebný atom zůstane).
# (zkratka, SMILES části, která se nahradí – první atom je ten navázaný na zbytek molekuly)
ZKRACOVANE_SKUPINY = [
    ("CoA", "CCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)(O)O"),
    ("CoA", "CCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)([O-])OP(=O)([O-])OCC1OC(n2cnc3c(N)ncnc32)C(O)C1OP(=O)([O-])[O-]"),
]
_VZORY_ZKRATEK = None


def zkrat_skupiny(mol):
    """Nahradí koenzym A (a podobné velké skupiny) jedním atomem s popiskem „CoA“."""
    global _VZORY_ZKRATEK
    if _VZORY_ZKRATEK is None:
        _VZORY_ZKRATEK = []
        for nazev, smi in ZKRACOVANE_SKUPINY:
            v = Chem.MolFromSmiles(smi)
            if v is not None:
                _VZORY_ZKRATEK.append((nazev, v))
    for nazev, vzor in _VZORY_ZKRATEK:
        while True:
            shoda = mol.GetSubstructMatch(vzor)
            if not shoda:
                break
            sada = set(shoda)
            kotva = None            # atom mimo skupinu navázaný na její první atom
            for nb in mol.GetAtomWithIdx(shoda[0]).GetNeighbors():
                if nb.GetIdx() not in sada:
                    kotva = nb.GetIdx()
            if kotva is None:
                break
            rw = Chem.RWMol(mol)
            novy = rw.AddAtom(Chem.Atom(0))
            rw.GetAtomWithIdx(novy).SetProp("atomLabel", nazev)
            rw.AddBond(kotva, novy, Chem.BondType.SINGLE)
            for idx in sorted(sada, reverse=True):
                rw.RemoveAtom(idx)
            mol = rw.GetMol()
            try:
                Chem.SanitizeMol(mol)
            except Exception:
                break
    return mol


ZKRACOVAT_COA = True       # False = koenzym A se kreslí celý
CELY_COA = ("CCNC(=O)CCNC(=O)C(O)C(C)(C)COP(=O)(O)OP(=O)(O)OCC%91OC(n%92cnc%93c(N)ncnc%93%92)"
            "C(O)C%91OP(=O)(O)O")


def nakresli_strukturu(zdroj):
    """zdroj = SMILES, MOL blok nebo KCF glykanu; vrací (šířka, výška, SVG obsah) nebo None."""
    if not zdroj:
        return None
    if zdroj.startswith("KCF"):
        return nakresli_glykan(zdroj)
    if not ZKRACOVAT_COA and "M  END" not in zdroj:
        zdroj = zdroj.replace("[*:1]", CELY_COA)          # vestavěné SMILES: rozepsat CoA
    mol = Chem.MolFromMolBlock(zdroj) if "M  END" in zdroj else Chem.MolFromSmiles(zdroj)
    if mol is None:
        return None
    if ZKRACOVAT_COA and mol.GetNumAtoms() > 30:
        mol = zkrat_skupiny(mol)
    for a in mol.GetAtoms():
        if a.GetAtomicNum() == 0 and a.GetAtomMapNum() in ZKRATKY:
            a.SetProp("atomLabel", ZKRATKY[a.GetAtomMapNum()])
            a.SetAtomMapNum(0)
    rdDepictor.SetPreferCoordGen(True)
    rdDepictor.Compute2DCoords(mol)
    _orientuj(mol)
    mol = rdMolDraw2D.PrepareMolForDrawing(mol)

    d = rdMolDraw2D.MolDraw2DSVG(-1, -1)       # plátno dle velikosti molekuly
    o = d.drawOptions()
    o.fixedBondLength = DELKA_VAZBY
    from matplotlib.colors import to_rgb
    cv = to_rgb(BARVA_VAZEB)
    if not BAREVNE_ATOMY:
        o.useBWAtomPalette()
        o.updateAtomPalette({k: cv for k in (-1, 1, 5, 6, 7, 8, 9, 11, 12, 14, 15, 16, 17, 19, 20,
                                               26, 29, 30, 34, 35, 53)})
    else:
        paleta = {-1: cv, 1: cv, 6: cv}
        if tmave_pozadi():                 # na tmavém pozadí světlejší heteroatomy
            paleta.update({7: (0.45, 0.66, 1.0), 8: (1.0, 0.45, 0.45), 16: (0.95, 0.85, 0.3),
                           15: (1.0, 0.65, 0.3)})
        o.updateAtomPalette(paleta)
    o.baseFontSize = FONT_ATOMU
    o.bondLineWidth = TLOUSTKA_VAZBY
    o.scaleBondWidth = False
    o.clearBackground = False
    o.padding = 0.04
    o.additionalAtomLabelPadding = 0.08
    d.DrawMolecule(mol)
    d.FinishDrawing()
    svg = d.GetDrawingText()
    w = float(re.search(r"width='([\d.]+)px'", svg).group(1))
    h = float(re.search(r"height='([\d.]+)px'", svg).group(1))
    body = svg.split("<!-- END OF HEADER -->", 1)[1].rsplit("</svg>", 1)[0]
    body = re.sub(r"<rect[^>]*/>", "", body)            # pro jistotu bez pozadí
    body = re.sub(r" class='[^']*'", "", body)          # kratší soubor
    return w, h, body.strip()


# =============================================================================
# 6) ROZLOŽENÍ GRAFU (Sugiyama)
# =============================================================================
class Uzel:
    def __init__(self, uid, w, h, dummy=False):
        self.id, self.w, self.h, self.dummy = uid, w, h, dummy
        self.layer = 0
        self.x = self.y = 0.0   # střed


def prirad_vrstvy(uzly, hrany):
    """Nejdelší cesta od zdrojů (topologicky)."""
    vstup = defaultdict(int)
    nasl = defaultdict(list)
    for u, v in hrany:
        vstup[v] += 1
        nasl[u].append(v)
    q = deque([u for u in uzly if vstup[u] == 0])
    vrstva = {u: 0 for u in uzly}
    seen = 0
    while q:
        u = q.popleft()
        seen += 1
        for v in nasl[u]:
            vrstva[v] = max(vrstva[v], vrstva[u] + 1)
            vstup[v] -= 1
            if vstup[v] == 0:
                q.append(v)
    if seen != len(uzly):
        raise ValueError("Graf reakcí obsahuje cyklus – zkuste otočit směr některé reakce.")
    return vrstva


def pocet_krizeni(horni, dolni, hrany_mezi):
    pos_h = {u: i for i, u in enumerate(horni)}
    pos_d = {u: i for i, u in enumerate(dolni)}
    e = sorted((pos_h[u], pos_d[v]) for u, v in hrany_mezi)
    c = 0
    for i in range(len(e)):
        for j in range(i + 1, len(e)):
            if e[i][0] < e[j][0] and e[i][1] > e[j][1]:
                c += 1
    return c


def celkove_krizeni(vrstvy, hrany_seg):
    return sum(pocet_krizeni(vrstvy[i], vrstvy[i + 1], hrany_seg[i])
               for i in range(len(vrstvy) - 1))


def serad_vrstvy(vrstvy, pred, nasl, hrany_seg, iterace=30):
    """Barycentrická heuristika + lokální výměny; vrací nejlepší pořadí."""
    def bary(layer, sousedi, ref):
        pos = {u: i for i, u in enumerate(ref)}
        klice = {}
        for i, u in enumerate(layer):
            s = [pos[v] for v in sousedi[u] if v in pos]
            klice[u] = (sum(s) / len(s)) if s else i
        return sorted(layer, key=lambda u: (klice[u], layer.index(u)))

    nej = [list(l) for l in vrstvy]
    nej_c = celkove_krizeni(nej, hrany_seg)
    cur = [list(l) for l in vrstvy]
    konec = time.time() + LIMIT_RAZENI_S       # u obřích map nečekat donekonečna
    for it in range(iterace):
        if time.time() > konec:
            break
        if it % 2 == 0:
            for i in range(1, len(cur)):
                cur[i] = bary(cur[i], pred, cur[i - 1])
        else:
            for i in range(len(cur) - 2, -1, -1):
                cur[i] = bary(cur[i], nasl, cur[i + 1])
        # transpozice sousedních uzlů (počítá se jen křížení sousedních vrstev)
        def mistni(i):
            c = 0
            if i > 0:
                c += pocet_krizeni(cur[i - 1], cur[i], hrany_seg[i - 1])
            if i < len(cur) - 1:
                c += pocet_krizeni(cur[i], cur[i + 1], hrany_seg[i])
            return c
        zlepseno = True
        while zlepseno and time.time() < konec:
            zlepseno = False
            for i, layer in enumerate(cur):
                for j in range(len(layer) - 1):
                    pred_c = mistni(i)
                    layer[j], layer[j + 1] = layer[j + 1], layer[j]
                    if mistni(i) < pred_c:
                        zlepseno = True
                    else:
                        layer[j], layer[j + 1] = layer[j + 1], layer[j]
        c = celkove_krizeni(cur, hrany_seg)
        if c < nej_c:
            nej, nej_c = [list(l) for l in cur], c
    return nej, nej_c


def pav_s_mezerami(zadane, vysky, mezera):
    """Izotonická regrese: najde y co nejblíže 'zadane' (středy) tak,
    aby se sousední uzly nepřekrývaly (y[i+1]-y[i] >= (h_i+h_{i+1})/2 + mezera)."""
    n = len(zadane)
    off = [0.0] * n
    for i in range(1, n):
        off[i] = off[i - 1] + (vysky[i - 1] + vysky[i]) / 2 + mezera
    z = [zadane[i] - off[i] for i in range(n)]
    bloky = []  # [součet, počet]
    for val in z:
        bloky.append([val, 1])
        while len(bloky) > 1 and bloky[-2][0] / bloky[-2][1] > bloky[-1][0] / bloky[-1][1]:
            s, c = bloky.pop()
            bloky[-1][0] += s
            bloky[-1][1] += c
    res = []
    for s, c in bloky:
        res += [s / c] * c
    return [res[i] + off[i] for i in range(n)]


# =============================================================================
# 6b) GEOMETRIE ŠIPEK
# =============================================================================
def bez(p0, p1, p2, p3, t):
    m = 1 - t
    return (m**3 * p0[0] + 3 * m * m * t * p1[0] + 3 * m * t * t * p2[0] + t**3 * p3[0],
            m**3 * p0[1] + 3 * m * m * t * p1[1] + 3 * m * t * t * p2[1] + t**3 * p3[1])


def segmenty(body):
    """Hladká křivka přes body s vodorovnými tečnami v každém bodě."""
    segs = []
    for (x0, y0), (x1, y1) in zip(body, body[1:]):
        k = (x1 - x0) * 0.5
        segs.append(((x0, y0), (x0 + k, y0), (x1 - k, y1), (x1, y1)))
    return segs


def d_krivky(segs):
    d = "M%s,%s" % (_fmt(segs[0][0][0]), _fmt(segs[0][0][1]))
    for _, c1, c2, p in segs:
        d += " C%s,%s %s,%s %s,%s" % tuple(_fmt(v) for v in (*c1, *c2, *p))
    return d


def d_lomene(pl):
    """SVG cesta pravoúhlé lomené čáry s malým zaoblením rohů."""
    d = "M%s,%s" % (_fmt(pl[0][0]), _fmt(pl[0][1]))
    for k in range(1, len(pl) - 1):
        (x0, y0), (x1, y1), (x2, y2) = pl[k - 1], pl[k], pl[k + 1]
        l1 = math.hypot(x1 - x0, y1 - y0)
        l2 = math.hypot(x2 - x1, y2 - y1)
        r = min(POLOMER_ROHU, l1 / 2, l2 / 2)
        if r < 0.5:
            d += " L%s,%s" % (_fmt(x1), _fmt(y1))
            continue
        ax, ay = x1 - (x1 - x0) / l1 * r, y1 - (y1 - y0) / l1 * r
        bx, by = x1 + (x2 - x1) / l2 * r, y1 + (y2 - y1) / l2 * r
        d += " L%s,%s Q%s,%s %s,%s" % tuple(_fmt(v) for v in (ax, ay, x1, y1, bx, by))
    d += " L%s,%s" % (_fmt(pl[-1][0]), _fmt(pl[-1][1]))
    return d


def posun_lomene(pl, o):
    """Rovnoběžná kopie pravoúhlé lomené čáry posunutá o 'o' doleva (ve směru jízdy)."""
    def normala(p, q):
        dx, dy = q[0] - p[0], q[1] - p[1]
        L = math.hypot(dx, dy) or 1
        return dy / L * o, -dx / L * o
    n = [normala(p, q) for p, q in zip(pl, pl[1:])]
    out = [(pl[0][0] + n[0][0], pl[0][1] + n[0][1])]
    for k in range(1, len(pl) - 1):
        out.append((pl[k][0] + n[k - 1][0] + n[k][0], pl[k][1] + n[k - 1][1] + n[k][1]))
    out.append((pl[-1][0] + n[-1][0], pl[-1][1] + n[-1][1]))
    return out


def hrot(tip, smer, pul=None, barva=None):
    """Hrot vodorovné šipky. smer = +1 (doprava) / -1 (doleva);
    pul = None (celý), 'nahore', 'dole' (poloviční hrot)."""
    strana = None
    if pul == "nahore":
        strana = "vlevo" if smer > 0 else "vpravo"
    elif pul == "dole":
        strana = "vpravo" if smer > 0 else "vlevo"
    return hrot_smer(tip, float(smer), 0.0, strana, barva)


def y_na_krivce(segs, x):
    """Najde y křivky v dané x-souřadnici (bisekce na monotónním segmentu)."""
    for s in segs:
        if min(s[0][0], s[3][0]) - 1e-6 <= x <= max(s[0][0], s[3][0]) + 1e-6:
            lo, hi = 0.0, 1.0
            for _ in range(50):
                mid = (lo + hi) / 2
                if bez(*s, mid)[0] < x:
                    lo = mid
                else:
                    hi = mid
            return bez(*s, (lo + hi) / 2)[1]
    return (segs[0][0][1] + segs[-1][3][1]) / 2


# =============================================================================
# 7) SESTAVENÍ SVG
# =============================================================================
def cesky_vynechano(n):
    if n == 1:
        return "vynechán 1 metabolit"
    if 2 <= n <= 4:
        return f"vynechány {n} metabolity"
    return f"vynecháno {n} metabolitů"


def zalom(slova, max_znaku=34):
    radky, akt = [], ""
    for s_ in slova:
        kandidat = (akt + ", " + s_) if akt else s_
        if len(kandidat) > max_znaku and akt:
            radky.append(akt + ",")
            akt = s_
        else:
            akt = kandidat
    if akt:
        radky.append(akt)
    return radky


def popisek_reakce(geny, ec, pritomna, vedl, zkratka=0):
    """Řádky popisku: gen, EC čísla (zalomená po 4), vedlejší produkty;
    u zkratky přes vynechané metabolity jen počet vynechaných kroků."""
    if zkratka:
        return []                   # šipky k rámečku vynechaných kroků jsou bez popisku
    zelena = pritomna and ZVYRAZNIT
    radky = []
    if geny:
        radky.append((Text(geny, VEL_GEN, weight="bold"), BARVA_ORG if zelena else BARVA_GEN))
    pol = [x for x in ec.split(" · ") if x]
    for j in range(0, len(pol), 4):
        radky.append((Text(" · ".join(pol[j:j + 4]), VEL_EC), BARVA_ORG if zelena else BARVA_EC))
    if vedl:
        radky.append((Text(vedl, VEL_EC, style="italic"), BARVA_POPIS))
    return radky


def priprav_data_uzlu(uzly_drahy, popisky="pod"):
    """Nakreslí struktury / rámečky všech uzlů a spočítá jejich rozměry."""
    data = {}
    for mid in uzly_drahy:
        if mid.startswith("map:"):          # rámeček navazující dráhy
            cislo = mid.split(":")[1]
            radky = [Text(x, VEL_MAPA) for x in MAPY.get(cislo, cislo).split("\n")]
            radky.append(Text(f"map{cislo}", VEL_EC, style="italic"))
            bw = max(t.w for t in radky) + 24
            bh = sum(t.h + 2 for t in radky) + 14
            body = [f'<rect x="0.6" y="0.6" width="{_fmt(bw - 1.2)}" height="{_fmt(bh - 1.2)}" '
                    f'rx="{_fmt(min(ZAOBLENI_RAMECKU * 1.4, bh / 2))}" fill="{BARVA_MAPY}" '
                    f'stroke="{BARVA_MAPY_OBRYS}" stroke-width="1.2"/>']
            yy = 7
            for j, t in enumerate(radky):
                barva = BARVA_MAPY_OBRYS if j == len(radky) - 1 else BARVA_MAPY_TEXT
                body.append(t.path((bw - t.w) / 2, yy, barva))
                yy += t.h + 2
            data[mid] = dict(sw=bw, sh=bh, body="".join(body), t_n=None, t_p=[], w=bw, h=bh)
            continue
        if mid.startswith("gap:"):          # rámeček vynechaných kroků
            jmena = [METABOLITY[c][0] for c in MEZERY[mid]]
            if len(jmena) > 6:
                jmena = jmena[:6] + [f"… (+{len(jmena) - 6})"]
            radky = [(Text("⋯ " + cesky_vynechano(len(MEZERY[mid])), VEL_EC + 1.5,
                           weight="bold", style="italic"), BARVA_ZKRATKA)]
            radky += [(Text(r_, VEL_EC), BARVA_POPIS) for r_ in zalom(jmena)]
            bw = max(t.w for t, _ in radky) + 22
            bh = sum(t.h + 2 for t, _ in radky) + 12
            body = [f'<rect x="0.6" y="0.6" width="{_fmt(bw - 1.2)}" height="{_fmt(bh - 1.2)}" '
                    f'rx="{_fmt(ZAOBLENI_RAMECKU * 0.8)}" fill="{BARVA_MEZERY}" stroke="{BARVA_ZKRATKA}" stroke-width="1.1" '
                    f'stroke-dasharray="4 3"/>']
            yy = 6
            for t, barva in radky:
                body.append(t.path((bw - t.w) / 2, yy, barva))
                yy += t.h + 2
            data[mid] = dict(sw=bw, sh=bh, body="".join(body), t_n=None, t_p=[], w=bw, h=bh)
            continue
        nazev, smi, pozn = METABOLITY[mid]
        st = nakresli_strukturu(smi)
        if st is None:      # struktura není k dispozici -> jen rámeček
            sw, sh = 70.0, 34.0
            body = ('<rect x="1" y="1" width="68" height="32" rx="6" fill="none" '
                    f'stroke="{BARVA_POPIS}" stroke-dasharray="3 3"/>')
        else:
            sw, sh, body = st
        t_n = Text(nazev, VEL_NAZEV, weight="bold")
        t_p = [Text(x.strip(), VEL_EC + 0.5, style="italic") for x in pozn.split(";")] if pozn else []
        extra_h = sum(t.h + 2 for t in t_p) + (2 if t_p else 0)
        pw = max([t.w for t in t_p], default=0)
        if popisky == "vpravo":
            w = sw + 8 + max(t_n.w, pw)
            h = max(sh, t_n.h + extra_h)
        else:
            w = max(sw, t_n.w, pw)
            h = sh + 4 + t_n.h + extra_h
        data[mid] = dict(sw=sw, sh=sh, body=body, t_n=t_n, t_p=t_p, w=w, h=h)
    return data


def sestav(uzly_drahy, popisky="pod"):
    """Řádkové rozložení (zleva doprava) s pravoúhlými šipkami."""
    data = priprav_data_uzlu(uzly_drahy, popisky)

    # --- graf ---
    hrany = [(r[0], r[1]) for r in REAKCE]
    vrstva = prirad_vrstvy(list(uzly_drahy), hrany)
    uzly = {m: Uzel(m, data[m]["w"], data[m]["h"]) for m in uzly_drahy}
    for m, u in uzly.items():
        u.layer = vrstva[m]

    retezy, seg_hrany, k = {}, [], 0
    for i, (a, b) in enumerate(hrany):
        cesta = [a]
        for L in range(vrstva[a] + 1, vrstva[b]):
            did = f"_d{k}"
            k += 1
            uzly[did] = Uzel(did, 0, 16, dummy=True)
            uzly[did].layer = L
            cesta.append(did)
        cesta.append(b)
        retezy[i] = cesta
        seg_hrany += list(zip(cesta, cesta[1:]))

    nL = max(u.layer for u in uzly.values()) + 1
    pred, nasl = defaultdict(list), defaultdict(list)
    for u, v in seg_hrany:
        nasl[u].append(v)
        pred[v].append(u)

    vrstvy = [[] for _ in range(nL)]
    videno = set()
    q = deque([m for m in uzly_drahy if vrstva[m] == 0])
    while q:
        u = q.popleft()
        if u in videno:
            continue
        videno.add(u)
        vrstvy[uzly[u].layer].append(u)
        q.extend(nasl[u])
    for m in uzly:                       # pro jistotu (nedosažitelné uzly)
        if m not in videno:
            vrstvy[uzly[m].layer].append(m)
    hrany_seg = [[(u, v) for u, v in seg_hrany if uzly[u].layer == i] for i in range(nL - 1)]
    vrstvy, _ = serad_vrstvy(vrstvy, pred, nasl, hrany_seg)

    # --- šířky vrstev a mezer ---
    popisky_reakci = [popisek_reakce(r[2], r[3], r[6], r[7], r[9]) for r in REAKCE]
    sirka_vrstvy = [max([uzly[u].w for u in l], default=0) for l in vrstvy]
    mezery = [MIN_MEZERA_VRSTEV] * (nL - 1)
    max_lw = defaultdict(float)
    for i, (a, b) in enumerate(hrany):
        lw = max([t.w for t, _ in popisky_reakci[i]], default=0)
        max_lw[vrstva[a]] = max(max_lw[vrstva[a]], lw)
    # mezera = pás kolejí (svislé úseky) + místo pro popisky na vodorovných úsecích
    for g in range(nL - 1):
        n = len(hrany_seg[g])
        pas = PAD_KOLEJ + max(n - 1, 0) * ROZESTUP_KOLEJI
        mezery[g] = max(MIN_MEZERA_VRSTEV, pas + max_lw[g] + 2 * OKRAJ_POPISKU + HROT_DELKA)

    top0 = OKRAJ + 80
    x = OKRAJ
    sloupec_x = []
    for i in range(nL):
        sloupec_x.append(x)
        for u in vrstvy[i]:
            uzly[u].x = x + sirka_vrstvy[i] / 2
        x += sirka_vrstvy[i] + (mezery[i] if i < nL - 1 else 0)
    celkova_sirka = x + OKRAJ

    # --- svislé souřadnice (PAV) ---
    for l in vrstvy:
        y = 0
        for u in l:
            uzly[u].y = y + uzly[u].h / 2
            y += uzly[u].h + MEZERA_SVISLE
    for it in range(120):
        for i in (range(nL) if it % 2 == 0 else range(nL - 1, -1, -1)):
            l = vrstvy[i]
            if not l:
                continue
            zadane = []
            for u in l:
                sous = pred[u] + nasl[u]
                if sous:
                    wsum = ysum = 0.0
                    for v in sous:
                        wt = 2.0 if (uzly[u].dummy or uzly[v].dummy) else 1.0
                        wsum += wt
                        ysum += wt * uzly[v].y
                    zadane.append(ysum / wsum)
                else:
                    zadane.append(uzly[u].y)
            nove = pav_s_mezerami(zadane, [uzly[u].h for u in l], MEZERA_SVISLE)
            for u, y in zip(l, nove):
                uzly[u].y = y
    # nesouvislé části grafu naskládat těsně pod sebe (bez prázdných mezer)
    sous_k = defaultdict(set)
    for u, v in seg_hrany:
        sous_k[u].add(v)
        sous_k[v].add(u)
    komponenty, vid = [], set()
    for u0 in uzly:
        if u0 in vid:
            continue
        komp, q = [], [u0]
        vid.add(u0)
        while q:
            u = q.pop()
            komp.append(u)
            for v in sous_k[u]:
                if v not in vid:
                    vid.add(v)
                    q.append(v)
        komponenty.append(komp)
    if len(komponenty) > 1:
        komponenty.sort(key=lambda k: (-len(k), min(uzly[u].y - uzly[u].h / 2 for u in k)))
        spodek = None
        for komp in komponenty:
            horni = min(uzly[u].y - uzly[u].h / 2 for u in komp)
            if spodek is not None:
                posun = spodek + 2 * MEZERA_SVISLE - horni
                for u in komp:
                    uzly[u].y += posun
            spodek = max(uzly[u].y + uzly[u].h / 2 for u in komp)
    ymin = min(u.y - u.h / 2 for u in uzly.values())
    for u in uzly.values():
        u.y += top0 - ymin
    ymax = max(u.y + u.h / 2 for u in uzly.values())

    def box_struktury(m):
        u, dd = uzly[m], data[m]
        if popisky == "vpravo":
            return u.x - u.w / 2, u.y - dd["sh"] / 2, dd["sw"], dd["sh"]
        return u.x - dd["sw"] / 2, u.y - u.h / 2, dd["sw"], dd["sh"]

    # porty – více šipek z/do jednoho uzlu se rozprostře svisle
    out_port, in_port = {}, {}
    for m in uzly_drahy:
        sx, sy, sw, sh = box_struktury(m)
        cy = sy + sh / 2
        for smer, slovnik in (("out", out_port), ("in", in_port)):
            idx = [i for i, (a, b) in enumerate(hrany) if (a if smer == "out" else b) == m]
            idx.sort(key=lambda i: uzly[retezy[i][1] if smer == "out" else retezy[i][-2]].y)
            n = len(idx)
            krok = min(14.0, max(9.0, sh * 0.9 / max(n, 1)))
            for j, i in enumerate(idx):
                slovnik[i] = cy + (j - (n - 1) / 2) * krok

    # --- šipky: pravoúhlé vedení s vlastní „kolejí“ pro každou hranu ---
    # Každý úsek hrany mezi dvěma sousedními vrstvami = vodorovně -> svisle -> vodorovně.
    # Svislá část leží v mezeře na vlastní koleji (x). Pořadí kolejí se volí tak,
    # aby bylo co nejméně průsečíků (lokální prohledávání přesunů).
    def koncove_body(i):
        a, b = hrany[i]
        sx, sy, sw, sh = box_struktury(a)
        tx, ty, tw, th = box_struktury(b)
        x_start = (uzly[a].x + uzly[a].w / 2) if popisky == "vpravo" else sx + sw
        return (x_start + ODSTUP_SIPKY, out_port[i]), (tx - ODSTUP_SIPKY, in_port[i])

    body_hran = {}
    hopy = defaultdict(list)          # mezera -> [(index hrany, pořadí úseku, ys, yt)]
    porty_cile = defaultdict(list)
    for i, (a, b) in enumerate(hrany):
        porty_cile[b].append(i)
    for i in range(len(hrany)):
        p0, p1 = koncove_body(i)
        cesta = retezy[i]
        pts = [p0] + [(uzly[d].x, uzly[d].y) for d in cesta[1:-1]] + [p1]
        # malé „schůdky“ (< 12 px) srovnat do rovné čáry, pokud tím cílový port
        # nesplyne s jiným portem téhož uzlu
        for k in range(len(pts) - 1):
            dy = pts[k + 1][1] - pts[k][1]
            if 0.5 < abs(dy) < 12:
                y_new = pts[k][1]
                if k + 1 == len(pts) - 1:
                    b = hrany[i][1]
                    sx_, sy_, sw_, sh_ = box_struktury(b)
                    volno = all(abs(in_port[j] - y_new) > 7 for j in porty_cile[b] if j != i)
                    if volno and sy_ + 4 < y_new < sy_ + sh_ - 4:
                        in_port[i] = y_new
                        pts[k + 1] = (pts[k + 1][0], y_new)
                else:
                    pts[k + 1] = (pts[k + 1][0], y_new)
        body_hran[i] = pts
        for k in range(len(pts) - 1):
            g = vrstva[hrany[i][0]] + k
            hopy[g].append((i, k, pts[k][1], pts[k + 1][1]))

    def cena(e, f):                   # počet křížení, když je e vlevo od f
        def uvnitr(y, a, b):
            lo, hi = min(a, b), max(a, b)
            return lo + 0.5 < y < hi - 0.5
        c = 0
        if uvnitr(f[2], e[2], e[3]):  # zdrojový úsek f protne svislici e
            c += 1
        if uvnitr(e[3], f[2], f[3]):  # cílový úsek e protne svislici f
            c += 1
        return c

    def usporadej(h):
        if len(h) < 2:
            return h
        o = sorted(range(len(h)), key=lambda i: (min(h[i][2], h[i][3]), h[i][3] - h[i][2]))
        n = len(h)
        C = [[cena(h[a], h[b]) for b in range(n)] for a in range(n)]   # C[a][b]: a vlevo od b
        konec_ = time.time() + 5
        zlepseno = True
        while zlepseno and time.time() < konec_:
            zlepseno = False
            for x in range(n):              # přesun prvku na lepší místo (počítá se jen změna)
                e = o[x]
                nej_d, nej_y, d = 0, x, 0
                for y in range(x - 1, -1, -1):          # doleva
                    f = o[y]
                    d += C[e][f] - C[f][e]
                    if d < nej_d:
                        nej_d, nej_y = d, y
                d = 0
                for y in range(x + 1, n):                # doprava
                    f = o[y]
                    d += C[f][e] - C[e][f]
                    if d < nej_d:
                        nej_d, nej_y = d, y
                if nej_y != x:
                    o.pop(x)
                    o.insert(nej_y, e)
                    zlepseno = True
        return [h[i] for i in o]

    kolej_x = {}
    for g, h in hopy.items():
        svisle = [e for e in h if abs(e[2] - e[3]) > 0.5]
        g_l = sloupec_x[g] + sirka_vrstvy[g]
        for j, e in enumerate(usporadej(svisle)):
            kolej_x[(e[0], e[1])] = g_l + PAD_KOLEJ + j * ROZESTUP_KOLEJI

    polyline = {}
    for i, pts in body_hran.items():
        pl = [pts[0]]
        for k in range(len(pts) - 1):
            (x0, y0), (x1, y1) = pts[k], pts[k + 1]
            if (i, k) in kolej_x:
                kx = kolej_x[(i, k)]
                pl += [(kx, y0), (kx, y1)]
            pl.append((x1, y1))
        # sloučit body ležící na jedné přímce
        cist = [pl[0]]
        for p in pl[1:]:
            if abs(p[0] - cist[-1][0]) < 1e-6 and abs(p[1] - cist[-1][1]) < 1e-6:
                continue
            if len(cist) >= 2:
                a_, b_ = cist[-2], cist[-1]
                if (abs(a_[0] - b_[0]) < 1e-6 and abs(b_[0] - p[0]) < 1e-6) or \
                   (abs(a_[1] - b_[1]) < 1e-6 and abs(b_[1] - p[1]) < 1e-6):
                    cist[-1] = p
                    continue
            cist.append(p)
        polyline[i] = cist

    sipky = []
    for i, (a, b, geny, ec, vratna, zpet, pritomna, vedl, carkovana, zkratka) in enumerate(REAKCE):
        if zkratka:                 # vynechané kroky: tečkovaná fialová šipka
            barva, sirka = BARVA_ZKRATKA, TLOUSTKA_SIPEK * 1.3
            carky = ' stroke-dasharray="0.1 4.5" stroke-linecap="round"'
        else:
            zelena = pritomna and ZVYRAZNIT
            barva = BARVA_ORG if zelena else BARVA_SIPKY
            sirka = TLOUSTKA_SIPEK * 1.4 if zelena else TLOUSTKA_SIPEK
            carky = ' stroke-dasharray="5 4"' if carkovana else ''

        def kresli(pl, dopredu, pul):
            pl = list(pl)
            if dopredu:
                tip = pl[-1]
                pl[-1] = (tip[0] - HROT_DELKA + 1.5, tip[1])
            else:
                tip = pl[0]
                pl[0] = (tip[0] + HROT_DELKA - 1.5, tip[1])
            sipky.append(f'<path d="{d_lomene(pl)}" fill="none" stroke="{barva}" '
                         f'stroke-width="{sirka}" stroke-linejoin="round"{carky}/>')
            sipky.append(hrot(tip, 1 if dopredu else -1, pul, barva))

        if vratna:
            kresli(posun_lomene(polyline[i], ROZESTUP_VRATNE), True, "nahore")
            kresli(posun_lomene(polyline[i], -ROZESTUP_VRATNE), False, "dole")
        else:
            kresli(polyline[i], not zpet, None)

    # skutečný počet průsečíků v hotovém obrázku
    def useky(pl):
        return list(zip(pl, pl[1:]))
    vsechny = [(i, s) for i, pl in polyline.items() for s in useky(pl)]
    krizeni = 0
    for x in range(len(vsechny)):
        i, (p, q) = vsechny[x]
        for y in range(x + 1, len(vsechny)):
            j, (r, t) = vsechny[y]
            if i == j:
                continue
            for (h1, h2), (v1, v2) in (((p, q), (r, t)), ((r, t), (p, q))):
                if abs(h1[1] - h2[1]) < 1e-6 and abs(v1[0] - v2[0]) < 1e-6:
                    if (min(h1[0], h2[0]) + 1 < v1[0] < max(h1[0], h2[0]) - 1 and
                            min(v1[1], v2[1]) + 1 < h1[1] < max(v1[1], v2[1]) - 1):
                        krizeni += 1

    # --- popisky enzymů s kontrolou kolizí ---
    def vzorkuj(pl):              # bod zhruba každé 4 px
        pts = []
        for (x0, y0), (x1, y1) in zip(pl, pl[1:]):
            n = max(2, int((abs(x1 - x0) + abs(y1 - y0)) / 4))
            pts += [(x0 + (x1 - x0) * t / n, y0 + (y1 - y0) * t / n) for t in range(n + 1)]
        return pts
    mrizka = defaultdict(list)       # body šipek v mřížce 16 px (rychlé hledání kolizí)
    for pl in polyline.values():
        for px, py in vzorkuj(pl):
            mrizka[(int(px // 16), int(py // 16))].append((px, py))

    def body_v_boxu(box):
        n_ = 0
        for gx_ in range(int(box[0] // 16), int(box[2] // 16) + 1):
            for gy_ in range(int(box[1] // 16), int(box[3] // 16) + 1):
                for px, py in mrizka.get((gx_, gy_), ()):
                    if box[0] <= px <= box[2] and box[1] <= py <= box[3]:
                        n_ += 1
        return n_
    boxy_uzlu = [(u.x - u.w / 2 - 4, u.y - u.h / 2 - 4, u.x + u.w / 2 + 4, u.y + u.h / 2 + 4)
                 for u in uzly.values() if not u.dummy]
    umistene, popisky_svg = [], []

    def prekryv(b1, b2):
        return not (b1[2] < b2[0] or b2[2] < b1[0] or b1[3] < b2[1] or b2[3] < b1[1])

    for i in sorted(range(len(REAKCE)), key=lambda i: vrstva[REAKCE[i][0]]):
        radky = popisky_reakci[i]
        if not radky:
            continue
        lw = max(t.w for t, _ in radky)
        lh = sum(t.h + 1 for t, _ in radky)
        pl = polyline[i]
        kandidati = []            # (x středu, y horního okraje, přirážka)
        n_us = len(pl) - 1
        for k, ((x0, y0), (x1, y1)) in enumerate(zip(pl, pl[1:])):
            if abs(y0 - y1) < 1e-6 and abs(x1 - x0) >= lw + 8:      # vodorovný úsek
                pref = 0 if k == n_us - 1 else 6                     # nejradši u cíle
                lo, hi = min(x0, x1) + lw / 2 + 4, max(x0, x1) - lw / 2 - 4 - HROT_DELKA
                for f in (0.5, 0.0, 1.0, 0.25, 0.75):
                    cx = lo + (hi - lo) * f
                    kandidati.append((cx, y0 - lh - 4, pref))
                    kandidati.append((cx, y0 + 5, pref + 3))
            elif abs(x0 - x1) < 1e-6 and abs(y1 - y0) >= lh:         # svislý úsek
                cy = (y0 + y1) / 2 - lh / 2
                kandidati.append((x0 + lw / 2 + 7, cy, 10))
                kandidati.append((x0 - lw / 2 - 7, cy, 12))
        if not kandidati:
            (x0, y0), (x1, y1) = pl[-2], pl[-1]
            kandidati.append(((x0 + x1) / 2, y1 - lh - 4, 0))
        nej = None
        for cx, ytop, pref in kandidati:
            box = (cx - lw / 2 - 3, ytop - 2, cx + lw / 2 + 3, ytop + lh + 2)
            skore = pref + 15 * body_v_boxu(box)
            skore += 500 * sum(prekryv(box, b2) for b2 in umistene)
            skore += 500 * sum(prekryv(box, b2) for b2 in boxy_uzlu)
            if nej is None or skore < nej[0]:
                nej = (skore, cx, ytop, box)
        _, cx, ytop, box = nej
        umistene.append(box)
        for t, barva in radky:
            popisky_svg.append(t.path(cx - t.w / 2, ytop, barva, halo=BARVA_POZADI))
            ytop += t.h + 1


    return dokonci_svg(uzly_drahy, uzly, data, popisky, sipky, popisky_svg, krizeni,
                       ymax, celkova_sirka)


def dokonci_svg(uzly_drahy, uzly, data, popisky, sipky, popisky_svg, krizeni, ymax, celkova_sirka):
    """Společné pro všechna rozložení: pozadí uzlů, struktury, názvy, nadpis, legenda."""
    def box_struktury(m):
        u, dd = uzly[m], data[m]
        if popisky == "vpravo":
            return u.x - u.w / 2, u.y - dd["sh"] / 2, dd["sw"], dd["sh"]
        return u.x - dd["sw"] / 2, u.y - u.h / 2, dd["sw"], dd["sh"]

    # --- pozadí uzlů ---
    pozadi = []
    for m in uzly_drahy:
        u = uzly[m]
        pad = 9
        rect = (f'x="{_fmt(u.x - u.w/2 - pad)}" y="{_fmt(u.y - u.h/2 - pad)}" '
                f'width="{_fmt(u.w + 2*pad)}" height="{_fmt(u.h + 2*pad)}" rx="{_fmt(ZAOBLENI_RAMECKU)}"')
        if ":" in m:
            continue
        if m == VYCHOZI:
            pozadi.append(f'<rect {rect} fill="{BARVA_VYCHOZI}" stroke="{BARVA_VYCHOZI_OBRYS}" stroke-width="1.2"/>')
        elif METABOLITY[m][2]:
            pozadi.append(f'<rect {rect} fill="none" stroke="{BARVA_POPIS}" stroke-width="1" '
                          f'stroke-dasharray="5 4"/>')

    # --- struktury a názvy ---
    struktury = []
    for m in uzly_drahy:
        dd, u = data[m], uzly[m]
        sx, sy, sw, sh = box_struktury(m)
        struktury.append(f'<g transform="translate({_fmt(sx)},{_fmt(sy)})">{dd["body"]}</g>')
        t_n, t_p = dd["t_n"], dd["t_p"]
        if t_n is None:
            continue
        if popisky == "vpravo":
            lx = sx + sw + 8
            ly = u.y - (t_n.h + sum(t.h + 2 for t in t_p)) / 2
        else:
            lx = u.x - t_n.w / 2
            ly = sy + sh + 4
        struktury.append(t_n.path(lx, ly, BARVA_NAZEV))
        ly += t_n.h + 3
        for t in t_p:
            px = lx if popisky == "vpravo" else u.x - t.w / 2
            struktury.append(t.path(px, ly, BARVA_POPIS))
            ly += t.h + 2

    # --- nadpis a legenda ---
    t1 = Text(NADPIS, VEL_NADPIS, weight="bold")
    t2 = Text(PODNADPIS, VEL_PODNADPIS)
    hlavicka = [t1.path(OKRAJ, OKRAJ - 8, BARVA_NADPIS), t2.path(OKRAJ, OKRAJ - 8 + t1.h + 4, BARVA_POPIS)]

    ly = ymax + 40
    leg, lx = [], OKRAJ

    def leg_text(t_, x):
        t = Text(t_, VEL_EC + 1.5)
        leg.append(t.path(x, ly - t.h / 2, BARVA_POPIS))
        return x + t.w + 34

    def leg_cara(x, barva, sirka=None, styl=""):
        sirka = sirka or TLOUSTKA_SIPEK
        leg.append(f'<path d="M{x},{ly} L{x+40},{ly}" fill="none" stroke="{barva}" '
                   f'stroke-width="{sirka}"{styl}/>')
        leg.append(hrot((x + 48, ly), 1, None, barva))
        return x + 56

    lx = leg_text("nevratná reakce", leg_cara(lx, BARVA_SIPKY))
    o = ROZESTUP_VRATNE
    leg.append(f'<path d="M{lx},{_fmt(ly-o)} L{lx+40},{_fmt(ly-o)}" stroke="{BARVA_SIPKY}" stroke-width="{TLOUSTKA_SIPEK}"/>')
    leg.append(hrot((lx + 48, ly - o), 1, "nahore", BARVA_SIPKY))
    leg.append(f'<path d="M{lx+8},{_fmt(ly+o)} L{lx+48},{_fmt(ly+o)}" stroke="{BARVA_SIPKY}" stroke-width="{TLOUSTKA_SIPEK}"/>')
    leg.append(hrot((lx, ly + o), -1, "dole", BARVA_SIPKY))
    lx = leg_text("vratná reakce", lx + 56)
    if any(r[8] and ":" not in r[0] + r[1] for r in REAKCE):
        lx = leg_text("enzym neznámý / spontánní",
                      leg_cara(lx, BARVA_SIPKY, TLOUSTKA_SIPEK, ' stroke-dasharray="5 4"'))
    if any(r[9] for r in REAKCE):
        lx = leg_text("vynechané kroky (rámeček uvádí skryté metabolity)",
                      leg_cara(lx, BARVA_ZKRATKA, TLOUSTKA_SIPEK * 1.3,
                               ' stroke-dasharray="0.1 4.5" stroke-linecap="round"'))
    if ZVYRAZNIT and any(r[6] for r in REAKCE):
        lx = leg_text(f"enzym přítomný v organismu {ORGANISMUS}", leg_cara(lx, BARVA_ORG, TLOUSTKA_SIPEK * 1.4))
    if any(u.startswith("map:") for u in uzly_drahy):
        tm = Text("Dráha", VEL_MAPA)
        leg.append(f'<rect x="{lx}" y="{_fmt(ly - 10)}" width="{_fmt(tm.w + 20)}" height="20" rx="10" '
                   f'fill="{BARVA_MAPY}" stroke="{BARVA_MAPY_OBRYS}" stroke-width="1.2"/>')
        leg.append(tm.path(lx + 10, ly - tm.h / 2, BARVA_MAPY_TEXT))
        lx = leg_text("navazující metabolická dráha", lx + tm.w + 30)
    for cukr in sorted(POUZITE_CUKRY):          # legenda symbolů SNFG
        leg.append(symbol_cukru(*SNFG[cukr], lx + 7, ly, GLYKAN_SYMBOL))
        lx = leg_text(cukr, lx + 20) - 18
    if POUZITE_CUKRY:
        lx += 16
    tg = Text("GEN", VEL_GEN, weight="bold")
    te = Text("1.2.3.4", VEL_EC)
    leg.append(tg.path(lx, ly - (tg.h + te.h) / 2, BARVA_GEN))
    leg.append(te.path(lx, ly - (tg.h + te.h) / 2 + tg.h + 1, BARVA_EC))
    konec_legendy = leg_text("gen (organismus) a EC číslo enzymu", lx + max(tg.w, te.w) + 10)

    W = math.ceil(max(celkova_sirka, OKRAJ * 2 + t1.w, OKRAJ * 2 + t2.w, konec_legendy - 34 + OKRAJ))
    H = math.ceil(ly + 30)
    svg = [
        '<?xml version="1.0" encoding="UTF-8"?>',
        f'<svg xmlns="http://www.w3.org/2000/svg" version="1.1" '
        f'width="{W}" height="{H}" viewBox="0 0 {W} {H}">',
        f'<title>{NADPIS}</title>',
        f'<rect x="0" y="0" width="{W}" height="{H}" fill="{BARVA_POZADI}"/>',
        '<g id="pozadi">', *pozadi, '</g>',
        '<g id="sipky">', *sipky, '</g>',
        '<g id="struktury">', *struktury, '</g>',
        '<g id="enzymy">', *popisky_svg, '</g>',
        '<g id="nadpis">', *hlavicka, '</g>',
        '<g id="legenda">', *leg, '</g>',
        '</svg>',
    ]
    return "\n".join(svg), krizeni


# =============================================================================
# 7b) KRUHOVÉ ROZLOŽENÍ PRO CYKLICKÉ DRÁHY (citrátový cyklus, močovinový cyklus …)
# =============================================================================
def smery_hrany(r):
    """Skutečné směry reakce z položky REAKCE: množina dvojic (z, do)."""
    a, b = r[0], r[1]
    if r[4]:
        return {(a, b), (b, a)}
    return {(b, a)} if r[5] else {(a, b)}


def najdi_cyklus(uzly_, hrany_smery, limit_s=2.0, min_delka=3):
    """Nejdelší jednoduchý orientovaný cyklus (prohledávání do hloubky s časovým limitem).
    hrany_smery: iterovatelné dvojice (z, do) ve směru reakcí."""
    nasl = defaultdict(set)
    for a, b in hrany_smery:
        if a != b:
            nasl[a].add(b)
    uzly_ = sorted(u for u in uzly_ if nasl[u])
    poradi_u = {u: i for i, u in enumerate(uzly_)}
    nej = []
    konec_ = time.time() + limit_s
    for s0 in uzly_:                    # cyklus začínající nejmenším uzlem (bez duplicit)
        if time.time() > konec_:
            break
        zasobnik = [(s0, iter(sorted(nasl[s0])))]
        cesta, na_ceste = [s0], {s0}
        while zasobnik and time.time() < konec_:
            u, it = zasobnik[-1]
            v = next(it, None)
            if v is None:
                zasobnik.pop()
                na_ceste.discard(cesta.pop())
                continue
            if v == s0 and len(cesta) >= min_delka and len(cesta) > len(nej):
                nej = list(cesta)
            elif v not in na_ceste and poradi_u.get(v, -1) > poradi_u[s0]:
                cesta.append(v)
                na_ceste.add(v)
                zasobnik.append((v, iter(sorted(nasl[v]))))
    return nej


def cyklus_v_draze(d):
    """Nejdelší cyklus v celé načtené dráze (pro nápovědu v okně)."""
    # přenašeče (ThPP, NAD …) se nepočítají; vázané meziprodukty ano (spojují kroky komplexu)
    uzly_ = {c for c, m in d.metabolity.items() if m.get("role", "normal") != "prenasec"}
    smery = []
    for r in d.reakce:
        if r["a"] in uzly_ and r["b"] in uzly_:
            smery.append((r["a"], r["b"]))
            if r["typ"] == "<=>":
                smery.append((r["b"], r["a"]))
    cyk = najdi_cyklus(sorted(uzly_), smery, limit_s=1.5, min_delka=4)
    return [c for c in cyk if d.metabolity[c].get("role", "normal") == "normal"]


def oblouk_d(cx, cy, r, a0, a1, zacatek=True):
    """Kruhový oblouk z úhlu a0 do a1 jako Bézierovy křivky (vektorově, bez příkazu A)."""
    n = max(1, int(math.ceil(abs(a1 - a0) / (math.pi / 2))))
    d = []
    if zacatek:
        d.append("M%s,%s" % (_fmt(cx + r * math.cos(a0)), _fmt(cy + r * math.sin(a0))))
    for k in range(n):
        t0 = a0 + (a1 - a0) * k / n
        t1 = a0 + (a1 - a0) * (k + 1) / n
        h = 4 / 3 * math.tan((t1 - t0) / 4)
        p0 = (cx + r * math.cos(t0), cy + r * math.sin(t0))
        p3 = (cx + r * math.cos(t1), cy + r * math.sin(t1))
        p1 = (p0[0] - h * r * math.sin(t0), p0[1] + h * r * math.cos(t0))
        p2 = (p3[0] + h * r * math.sin(t1), p3[1] - h * r * math.cos(t1))
        d.append(" C%s,%s %s,%s %s,%s" % tuple(_fmt(v) for v in (*p1, *p2, *p3)))
    return "".join(d)


def hrot_smer(tip, ux, uy, pul=None, barva=None):
    """Hrot šipky v libovolném směru (ux, uy = jednotkový směr jízdy).
    pul = None (celý), "vlevo" / "vpravo" (poloviční hrot na dané straně vzhledem ke směru).
    Tvar podle TVAR_HROTU: plny (se zářezem), trojuhelnik, otevreny (čárový)."""
    barva = barva or BARVA_SIPKY
    x, y = tip
    nx, ny = uy, -ux                      # levá normála
    if pul == "vpravo":
        nx, ny = -nx, -ny
    bx, by = x - ux * HROT_DELKA, y - uy * HROT_DELKA
    f = _fmt
    if TVAR_HROTU == "otevreny":
        if pul is None:
            d = (f"M{f(bx + nx * HROT_SIRKA)},{f(by + ny * HROT_SIRKA)} L{f(x)},{f(y)} "
                 f"L{f(bx - nx * HROT_SIRKA)},{f(by - ny * HROT_SIRKA)} M{f(bx)},{f(by)} L{f(x)},{f(y)}")
        else:
            d = (f"M{f(bx + nx * HROT_SIRKA * 1.4)},{f(by + ny * HROT_SIRKA * 1.4)} L{f(x)},{f(y)} "
                 f"M{f(bx)},{f(by)} L{f(x)},{f(y)}")
        return (f'<path d="{d}" fill="none" stroke="{barva}" stroke-width="{f(TLOUSTKA_SIPEK)}" '
                f'stroke-linecap="round" stroke-linejoin="round"/>')
    if pul is None:
        if TVAR_HROTU == "trojuhelnik":
            pts = [(x, y), (bx + nx * HROT_SIRKA, by + ny * HROT_SIRKA), (bx - nx * HROT_SIRKA, by - ny * HROT_SIRKA)]
        else:
            pts = [(x, y), (bx + nx * HROT_SIRKA, by + ny * HROT_SIRKA),
                   (bx + ux * HROT_DELKA * 0.22, by + uy * HROT_DELKA * 0.22),
                   (bx - nx * HROT_SIRKA, by - ny * HROT_SIRKA)]
    else:
        pts = [(x, y), (bx + nx * HROT_SIRKA * 1.4, by + ny * HROT_SIRKA * 1.4), (bx, by)]
    d = "M" + " L".join("%s,%s" % (f(a), f(b)) for a, b in pts) + " Z"
    return f'<path d="{d}" fill="{barva}"/>'


def sestav_kruh(uzly_drahy, popisky="pod"):
    """Kruhové rozložení: nejdelší cyklus po kružnici (po směru hodinových ručiček),
    ostatní látky a rámečky paprskovitě ven od kružnice."""
    data = priprav_data_uzlu(uzly_drahy, popisky)
    uzly = {m: Uzel(m, data[m]["w"], data[m]["h"]) for m in uzly_drahy}

    smery = set()
    for r in REAKCE:
        if ":" in r[0] + r[1] and not (r[0].startswith("gap:") or r[1].startswith("gap:")):
            continue                              # rámečky drah do cyklu nepatří
        smery |= smery_hrany(r)
    cyklus = najdi_cyklus([u for u in uzly_drahy if not u.startswith("map:")], smery)
    if len(cyklus) < 3:
        raise ValueError("Ve vybrané části dráhy nebyl nalezen žádný cyklus (alespoň 3 látky "
                         "spojené reakcemi do kruhu).\n\nPoužij řádkové rozložení, nebo zkontroluj "
                         "směry reakcí (karta Reakce → Otočit směr).")
    # nahoru dát výchozí látku, nebo látku cyklu, která je k ní nejblíž
    if VYCHOZI and VYCHOZI in uzly_drahy:
        sous_v = defaultdict(set)
        for r in REAKCE:
            sous_v[r[0]].add(r[1])
            sous_v[r[1]].add(r[0])
        vzd, q = {VYCHOZI: 0}, deque([VYCHOZI])
        while q:
            u = q.popleft()
            for v in sous_v[u]:
                if v not in vzd:
                    vzd[v] = vzd[u] + 1
                    q.append(v)
        nejbliz = min(cyklus, key=lambda c: (vzd.get(c, 999), cyklus.index(c)))
        k = cyklus.index(nejbliz)
        cyklus = cyklus[k:] + cyklus[:k]
    n = len(cyklus)
    ncyk = set(cyklus)

    # --- poloměr kružnice podle velikosti uzlů ---
    obvod = sum(math.hypot(uzly[c].w, uzly[c].h) * 0.85 + 125 for c in cyklus)   # + místo na šipku a popisek
    R = max(obvod / (2 * math.pi), 170.0)
    uhel = {c: -math.pi / 2 + 2 * math.pi * i / n for i, c in enumerate(cyklus)}
    cx0 = cy0 = 0.0
    for c in cyklus:
        uzly[c].x, uzly[c].y = cx0 + R * math.cos(uhel[c]), cy0 + R * math.sin(uhel[c])

    # --- ostatní uzly: strom od kružnice ven ---
    sous = defaultdict(set)
    for r in REAKCE:
        sous[r[0]].add(r[1])
        sous[r[1]].add(r[0])
    rodic, hloubka, deti = {}, {c: 0 for c in cyklus}, defaultdict(list)
    q = deque(cyklus)
    while q:
        u = q.popleft()
        for v in sorted(sous[u], key=lambda m: uzly_drahy.index(m)):
            if v not in hloubka:
                hloubka[v] = hloubka[u] + 1
                rodic[v] = u
                deti[u].append(v)
                q.append(v)

    def listy(u):
        return max(1, sum(listy(v) for v in deti[u]))

    krok = max([math.hypot(uzly[u].w, uzly[u].h) for u in uzly_drahy], default=150) * 0.75 + 90
    sektor = {}
    for c in cyklus:
        sektor[c] = (uhel[c] - math.pi / n * 0.92, uhel[c] + math.pi / n * 0.92)

    def rozmisti(u):
        a0, a1 = sektor[u]
        ch = deti[u]
        if not ch:
            return
        celkem = sum(listy(v) for v in ch)
        akt = a0
        for v in ch:
            sirka = (a1 - a0) * listy(v) / celkem
            sektor[v] = (akt, akt + sirka)
            stred = akt + sirka / 2
            rr = R + hloubka[v] * krok
            uzly[v].x, uzly[v].y = cx0 + rr * math.cos(stred), cy0 + rr * math.sin(stred)
            akt += sirka
            rozmisti(v)
    for c in cyklus:
        rozmisti(c)

    # látka napojená na víc blízkých látek cyklu (např. acetyl-CoA k oxaloacetátu i citrátu)
    # se natočí doprostřed mezi ně – i s celou svou větví
    def potomci(u):
        out = [u]
        for v in deti[u]:
            out += potomci(v)
        return out
    for c in cyklus:
        for v in list(deti[c]):
            sousede_c = [x for x in sous[v] if x in ncyk]
            if len(sousede_c) < 2:
                continue
            sx_ = sum(math.cos(uhel[x]) for x in sousede_c)
            sy_ = sum(math.sin(uhel[x]) for x in sousede_c)
            if math.hypot(sx_, sy_) < 0.5 * len(sousede_c):
                continue                    # sousedé jsou na opačných stranách kruhu
            cil = math.atan2(sy_, sx_)
            akt = math.atan2(uzly[v].y - cy0, uzly[v].x - cx0)
            dl = cil - akt
            cs, sn = math.cos(dl), math.sin(dl)
            for u in potomci(v):
                x_, y_ = uzly[u].x - cx0, uzly[u].y - cy0
                uzly[u].x, uzly[u].y = cx0 + x_ * cs - y_ * sn, cy0 + x_ * sn + y_ * cs

    # uzly nenapojené na cyklus: do řady pod kružnici
    volne = [u for u in uzly_drahy if u not in hloubka]
    if volne:
        y_rad = max(uzly[u].y + uzly[u].h / 2 for u in hloubka) + 120
        x_akt = -R
        for u in volne:
            uzly[u].x, uzly[u].y = x_akt + uzly[u].w / 2, y_rad
            x_akt += uzly[u].w + 60

    # odstranit překryvy: necyklové uzly se odsouvají po svém paprsku ven
    def prekryv_uzlu(a, b, okraj=26):
        return (abs(a.x - b.x) < (a.w + b.w) / 2 + okraj and
                abs(a.y - b.y) < (a.h + b.h) / 2 + okraj)
    seznam = list(uzly_drahy)
    for _ in range(300):
        pohyb = False
        for i, a in enumerate(seznam):
            for b in seznam[i + 1:]:
                ua, ub = uzly[a], uzly[b]
                if not prekryv_uzlu(ua, ub):
                    continue
                kdo = b if (b not in ncyk and (a in ncyk or hloubka.get(b, 99) >= hloubka.get(a, 99))) else a
                if kdo in ncyk:
                    continue
                u = uzly[kdo]
                dx, dy = u.x - cx0, u.y - cy0
                L = math.hypot(dx, dy) or 1
                u.x += dx / L * 18
                u.y += dy / L * 18
                pohyb = True
        if not pohyb:
            break

    # posun celého obrázku do kladných souřadnic
    top0 = OKRAJ + 80
    xmin = min(u.x - u.w / 2 for u in uzly.values()) - 60
    ymin = min(u.y - u.h / 2 for u in uzly.values()) - 60
    for u in uzly.values():
        u.x += OKRAJ - xmin
        u.y += top0 - ymin
    cx0 += OKRAJ - xmin
    cy0 += top0 - ymin

    def box_struktury(m):
        u, dd = uzly[m], data[m]
        if popisky == "vpravo":
            return u.x - u.w / 2, u.y - dd["sh"] / 2, dd["sw"], dd["sh"]
        return u.x - dd["sw"] / 2, u.y - u.h / 2, dd["sw"], dd["sh"]

    def stred_struktury(m):
        sx, sy, sw, sh = box_struktury(m)
        return sx + sw / 2, sy + sh / 2

    def uvnitr_uzlu(m, x, y, pad=8):
        u = uzly[m]
        return abs(x - u.x) <= u.w / 2 + pad and abs(y - u.y) <= u.h / 2 + pad

    def okraj_bod(m, x0, y0, x1, y1):
        """Bod, kde úsečka (x0,y0)->(x1,y1) opouští uzel m."""
        lo, hi = 0.0, 1.0
        for _ in range(40):
            mid = (lo + hi) / 2
            if uvnitr_uzlu(m, x0 + (x1 - x0) * mid, y0 + (y1 - y0) * mid):
                lo = mid
            else:
                hi = mid
        return x0 + (x1 - x0) * hi, y0 + (y1 - y0) * hi

    # --- šipky ---
    sipky, cary = [], {}          # cary: index -> body pro kolize a popisky
    nakreslene = set()

    def styl(r):
        a, b, geny, ec, vratna, zpet, pritomna, vedl, carkovana, zkratka = r
        if zkratka:
            return BARVA_ZKRATKA, TLOUSTKA_SIPEK * 1.3, ' stroke-dasharray="0.1 4.5" stroke-linecap="round"'
        zel = pritomna and ZVYRAZNIT
        return (BARVA_ORG if zel else BARVA_SIPKY, TLOUSTKA_SIPEK * 1.4 if zel else TLOUSTKA_SIPEK,
                ' stroke-dasharray="5 4"' if carkovana else '')

    # a) hrany po kružnici mezi sousedy v cyklu
    def konec_oblouku(c, a_start, smer):
        """Úhel, kde oblouk opouští uzel c (binární hledání)."""
        lo, hi = 0.0, math.pi / n
        for _ in range(40):
            mid = (lo + hi) / 2
            t = a_start + smer * mid
            if uvnitr_uzlu(c, cx0 + R * math.cos(t), cy0 + R * math.sin(t)):
                lo = mid
            else:
                hi = mid
        return a_start + smer * hi

    for i, c in enumerate(cyklus):
        c2 = cyklus[(i + 1) % n]
        idx = [j for j, r in enumerate(REAKCE) if {r[0], r[1]} == {c, c2}]
        if not idx:
            continue
        j = idx[0]
        nakreslene.add(j)
        r = REAKCE[j]
        barva, sirka, carky = styl(r)
        a0 = uhel[c]
        a1 = uhel[c2] if uhel[c2] > a0 else uhel[c2] + 2 * math.pi
        t0 = konec_oblouku(c, a0, +1) + 0.01
        t1 = konec_oblouku(c2, a1, -1) - 0.01
        dozadu_delka = HROT_DELKA / R
        smery_r = smery_hrany(r)
        vratna = (c, c2) in smery_r and (c2, c) in smery_r
        if vratna:
            o = ROZESTUP_VRATNE
            for rr, od, do_, smer in ((R + o, t0, t1 - dozadu_delka, +1), (R - o, t1, t0 + dozadu_delka, -1)):
                sipky.append(f'<path d="{oblouk_d(cx0, cy0, rr, od, do_)}" fill="none" stroke="{barva}" '
                             f'stroke-width="{sirka}"{carky}/>')
                tip_u = do_ + smer * dozadu_delka
                tip = (cx0 + rr * math.cos(tip_u), cy0 + rr * math.sin(tip_u))
                ux, uy = -math.sin(tip_u) * smer, math.cos(tip_u) * smer
                sipky.append(hrot_smer(tip, ux, uy, "vlevo", barva))
        else:
            if (c, c2) in smery_r:
                od, do_, smer = t0, t1, +1
            else:
                od, do_, smer = t1, t0, -1
            sipky.append(f'<path d="{oblouk_d(cx0, cy0, R, od, do_ - smer * dozadu_delka)}" fill="none" '
                         f'stroke="{barva}" stroke-width="{sirka}"{carky}/>')
            tip = (cx0 + R * math.cos(do_), cy0 + R * math.sin(do_))
            ux, uy = -math.sin(do_) * smer, math.cos(do_) * smer
            sipky.append(hrot_smer(tip, ux, uy, None, barva))
        pts = [(cx0 + R * math.cos(t0 + (t1 - t0) * k / 40), cy0 + R * math.sin(t0 + (t1 - t0) * k / 40))
               for k in range(41)]
        cary[j] = ("oblouk", pts, (t0 + t1) / 2)

    # b) spojnice mezi blízkými (nesousedními) látkami cyklu: oblouk po vnitřní straně kruhu
    poz_c = {c: i for i, c in enumerate(cyklus)}
    R_in = R - max(math.hypot(uzly[c].w, uzly[c].h) / 2 for c in cyklus) - 22
    for j, r in enumerate(REAKCE):
        a, b = r[0], r[1]
        if j in nakreslene or a not in ncyk or b not in ncyk or R_in < 60:
            continue
        krok_c = (poz_c[b] - poz_c[a]) % n
        if min(krok_c, n - krok_c) > max(2, n // 3):
            continue                        # vzdálené látky -> rovná čára (níže)
        smery_r = smery_hrany(r)
        u_, v_ = (a, b) if (a, b) in smery_r else (b, a)
        ta = uhel[u_]
        tb = uhel[v_]
        # kratší směr po kruhu
        dt = (tb - ta) % (2 * math.pi)
        if dt > math.pi:
            dt -= 2 * math.pi
        tb = ta + dt
        barva, sirka, carky = styl(r)
        (xa, ya), (xb, yb) = stred_struktury(u_), stred_struktury(v_)
        p_in_a = (cx0 + R_in * math.cos(ta), cy0 + R_in * math.sin(ta))
        p_in_b = (cx0 + R_in * math.cos(tb), cy0 + R_in * math.sin(tb))
        pa = okraj_bod(u_, xa, ya, *p_in_a)
        pb = okraj_bod(v_, xb, yb, *p_in_b)
        L = math.hypot(pb[0] - p_in_b[0], pb[1] - p_in_b[1]) or 1
        ux, uy = (pb[0] - p_in_b[0]) / L, (pb[1] - p_in_b[1]) / L
        konec = (pb[0] - ux * (HROT_DELKA - 1.5), pb[1] - uy * (HROT_DELKA - 1.5))
        d_ = (f"M{_fmt(pa[0])},{_fmt(pa[1])} L{_fmt(p_in_a[0])},{_fmt(p_in_a[1])}"
              + oblouk_d(cx0, cy0, R_in, ta, tb, zacatek=False)
              + f" L{_fmt(konec[0])},{_fmt(konec[1])}")
        sipky.append(f'<path d="{d_}" fill="none" stroke="{barva}" stroke-width="{sirka}" '
                     f'stroke-linejoin="round"{carky}/>')
        sipky.append(hrot_smer(pb, ux, uy, None, barva))
        if (v_, u_) in smery_r:             # vratná: hrot i na začátku
            L2 = math.hypot(pa[0] - p_in_a[0], pa[1] - p_in_a[1]) or 1
            sipky.append(hrot_smer(pa, (pa[0] - p_in_a[0]) / L2, (pa[1] - p_in_a[1]) / L2, None, barva))
        nakreslene.add(j)
        pts = [(cx0 + R_in * math.cos(ta + (tb - ta) * k / 30), cy0 + R_in * math.sin(ta + (tb - ta) * k / 30))
               for k in range(31)]
        cary[j] = ("oblouk_in", pts, (ta + tb) / 2)

    # c) ostatní hrany: rovné čáry mezi okraji uzlů
    def okraj(m, x0, y0, x1, y1):
        """Bod, kde úsečka (x0,y0)->(x1,y1) opouští uzel m."""
        lo, hi = 0.0, 1.0
        for _ in range(40):
            mid = (lo + hi) / 2
            if uvnitr_uzlu(m, x0 + (x1 - x0) * mid, y0 + (y1 - y0) * mid):
                lo = mid
            else:
                hi = mid
        return x0 + (x1 - x0) * hi, y0 + (y1 - y0) * hi

    for j, r in enumerate(REAKCE):
        if j in nakreslene:
            continue
        a, b = r[0], r[1]
        barva, sirka, carky = styl(r)
        (xa, ya), (xb, yb) = stred_struktury(a), stred_struktury(b)
        pa = okraj(a, xa, ya, xb, yb)
        pb = okraj(b, xb, yb, xa, ya)
        L = math.hypot(pb[0] - pa[0], pb[1] - pa[1])
        if L < HROT_DELKA * 2:
            continue
        ux, uy = (pb[0] - pa[0]) / L, (pb[1] - pa[1]) / L
        smery_r = smery_hrany(r)
        vratna = (a, b) in smery_r and (b, a) in smery_r

        def usecka(p, q, ux_, uy_, pul):
            konec = (q[0] - ux_ * (HROT_DELKA - 1.5), q[1] - uy_ * (HROT_DELKA - 1.5))
            sipky.append(f'<path d="M{_fmt(p[0])},{_fmt(p[1])} L{_fmt(konec[0])},{_fmt(konec[1])}" '
                         f'fill="none" stroke="{barva}" stroke-width="{sirka}"{carky}/>')
            sipky.append(hrot_smer(q, ux_, uy_, pul, barva))
        if vratna:
            o = ROZESTUP_VRATNE
            nx, ny = uy * o, -ux * o
            usecka((pa[0] + nx, pa[1] + ny), (pb[0] + nx, pb[1] + ny), ux, uy, "vlevo")
            usecka((pb[0] - nx, pb[1] - ny), (pa[0] - nx, pa[1] - ny), -ux, -uy, "vlevo")
        elif (a, b) in smery_r:
            usecka(pa, pb, ux, uy, None)
        else:
            usecka(pb, pa, -ux, -uy, None)
        pts = [(pa[0] + (pb[0] - pa[0]) * k / 30, pa[1] + (pb[1] - pa[1]) * k / 30) for k in range(31)]
        cary[j] = ("usecka", pts, None)

    # --- popisky enzymů (kandidáti kolem středu čáry, kontrola kolizí) ---
    popisky_reakci = [popisek_reakce(r[2], r[3], r[6], r[7], r[9]) for r in REAKCE]
    mrizka = defaultdict(list)
    for _, pts, _ in cary.values():
        for px, py in pts:
            mrizka[(int(px // 16), int(py // 16))].append((px, py))

    def body_v_boxu(box):
        n_ = 0
        for gx_ in range(int(box[0] // 16), int(box[2] // 16) + 1):
            for gy_ in range(int(box[1] // 16), int(box[3] // 16) + 1):
                for px, py in mrizka.get((gx_, gy_), ()):
                    if box[0] <= px <= box[2] and box[1] <= py <= box[3]:
                        n_ += 1
        return n_
    boxy_uzlu = [(u.x - u.w / 2 - 4, u.y - u.h / 2 - 4, u.x + u.w / 2 + 4, u.y + u.h / 2 + 4)
                 for u in uzly.values()]

    def prekryv(b1, b2):
        return not (b1[2] < b2[0] or b2[2] < b1[0] or b1[3] < b2[1] or b2[3] < b1[1])

    umistene, popisky_svg = [], []
    for j, (druh, pts, ustred) in cary.items():
        radky = popisky_reakci[j]
        if not radky:
            continue
        lw = max(t.w for t, _ in radky)
        lh = sum(t.h + 1 for t, _ in radky)
        kandidati = []
        if druh == "oblouk_in":
            for du in (0, -0.15, 0.15):
                t = ustred + du * (2 * math.pi / n)
                ox, oy = math.cos(t), math.sin(t)
                vzd = abs(ox) * lw / 2 + abs(oy) * lh / 2
                cxl = cx0 + (R_in - 10 - vzd) * ox
                cyl = cy0 + (R_in - 10 - vzd) * oy
                kandidati.append((cxl, cyl - lh / 2, abs(du) * 20))
        elif druh == "oblouk":
            for du in (0, -0.12, 0.12):
                t = ustred + du * (2 * math.pi / n)
                for rr, pref in ((R + 12, 0), (R - 12, 4)):
                    # střed popisku posunutý ven/dovnitř o polovinu jeho rozměru
                    ox, oy = math.cos(t), math.sin(t)
                    vzd = abs(ox) * lw / 2 + abs(oy) * lh / 2
                    znam = 1 if rr > R else -1
                    cxl = cx0 + (R + znam * (12 + vzd)) * ox
                    cyl = cy0 + (R + znam * (12 + vzd)) * oy
                    kandidati.append((cxl, cyl - lh / 2, pref + abs(du) * 20))
        else:
            (x0, y0), (x1, y1) = pts[0], pts[-1]
            L = math.hypot(x1 - x0, y1 - y0) or 1
            nx, ny = (y1 - y0) / L, -(x1 - x0) / L
            for f in (0.5, 0.35, 0.65):
                mx, my = x0 + (x1 - x0) * f, y0 + (y1 - y0) * f
                for znam, pref in ((1, 0), (-1, 2)):
                    vzd = abs(nx) * lw / 2 + abs(ny) * lh / 2 + 6
                    kandidati.append((mx + znam * nx * vzd, my + znam * ny * vzd - lh / 2,
                                      pref + abs(f - 0.5) * 10))
        nej = None
        for cxl, ytop, pref in kandidati:
            box = (cxl - lw / 2 - 3, ytop - 2, cxl + lw / 2 + 3, ytop + lh + 2)
            skore = pref + 15 * body_v_boxu(box)
            skore += 500 * sum(prekryv(box, b2) for b2 in umistene)
            skore += 500 * sum(prekryv(box, b2) for b2 in boxy_uzlu)
            if nej is None or skore < nej[0]:
                nej = (skore, cxl, ytop, box)
        _, cxl, ytop, box = nej
        umistene.append(box)
        for t, barva in radky:
            popisky_svg.append(t.path(cxl - t.w / 2, ytop, barva, halo=BARVA_POZADI))
            ytop += t.h + 1

    # skutečná křížení rovných čar (oblouky cyklu se nekříží)
    usecky = [(j, pts[0], pts[-1]) for j, (druh, pts, _) in cary.items() if druh == "usecka"]

    def protina(p1, p2, p3, p4):
        def o(a, b, c):
            return (b[0] - a[0]) * (c[1] - a[1]) - (b[1] - a[1]) * (c[0] - a[0])
        return (o(p1, p2, p3) * o(p1, p2, p4) < 0) and (o(p3, p4, p1) * o(p3, p4, p2) < 0)
    krizeni = sum(1 for x in range(len(usecky)) for y in range(x + 1, len(usecky))
                  if protina(usecky[x][1], usecky[x][2], usecky[y][1], usecky[y][2]))

    ymax = max([u.y + u.h / 2 for u in uzly.values()] + [b[3] for b in umistene])
    sirka = max([u.x + u.w / 2 for u in uzly.values()] + [b[2] for b in umistene]) + OKRAJ
    return dokonci_svg(uzly_drahy, uzly, data, popisky, sipky, popisky_svg, krizeni, ymax, sirka)


# =============================================================================
# 7c) NÁHLED: vykreslení vlastního SVG do obrázku (matplotlib Agg, bez dalších knihoven)
# =============================================================================
def _barva(hodnota):
    if not hodnota or hodnota == "none":
        return None
    from matplotlib.colors import to_rgba
    try:
        return to_rgba(hodnota)
    except ValueError:
        return None


def _cesta_z_d(d):
    """SVG atribut d -> matplotlib Path (příkazy M L H V C Q Z A, abs. i rel.)."""
    from matplotlib.path import Path
    tok = re.findall(r"[MLHVCQZAmlhvcqza]|[-+]?(?:\d*\.\d+|\d+\.?)(?:[eE][-+]?\d+)?", d)
    v, c = [], []
    i, cmd = 0, None
    x = y = sx = sy = 0.0

    def cislo():
        nonlocal i
        i += 1
        return float(tok[i - 1])
    while i < len(tok):
        if re.match(r"[A-Za-z]", tok[i]):
            cmd = tok[i]
            i += 1
            if cmd in "Zz":
                if v:
                    v.append((sx, sy))
                    c.append(Path.CLOSEPOLY)
                x, y = sx, sy
                continue
        rel = cmd.islower()
        C = cmd.upper()
        if C == "M":
            nx_, ny_ = cislo(), cislo()
            x, y = (x + nx_, y + ny_) if rel else (nx_, ny_)
            sx, sy = x, y
            v.append((x, y))
            c.append(Path.MOVETO)
            cmd = "l" if rel else "L"
        elif C == "L":
            nx_, ny_ = cislo(), cislo()
            x, y = (x + nx_, y + ny_) if rel else (nx_, ny_)
            v.append((x, y))
            c.append(Path.LINETO)
        elif C == "H":
            nx_ = cislo()
            x = x + nx_ if rel else nx_
            v.append((x, y))
            c.append(Path.LINETO)
        elif C == "V":
            ny_ = cislo()
            y = y + ny_ if rel else ny_
            v.append((x, y))
            c.append(Path.LINETO)
        elif C in "CQ":
            k = 3 if C == "C" else 2
            body_ = []
            for _ in range(k):
                px, py = cislo(), cislo()
                body_.append((x + px, y + py) if rel else (px, py))
            v += body_
            c += [Path.CURVE4 if C == "C" else Path.CURVE3] * k
            x, y = body_[-1]
        elif C == "A":                      # kruhový/eliptický oblouk -> úsečky
            rx, ry, rot, velky, smer_, ex, ey = [cislo() for _ in range(7)]
            if rel:
                ex, ey = x + ex, y + ey
            if rx == 0 or ry == 0:
                v.append((ex, ey))
                c.append(Path.LINETO)
            else:
                fi = math.radians(rot)
                cf, sf = math.cos(fi), math.sin(fi)
                dx2, dy2 = (x - ex) / 2, (y - ey) / 2
                x1p, y1p = cf * dx2 + sf * dy2, -sf * dx2 + cf * dy2
                lam = (x1p / rx) ** 2 + (y1p / ry) ** 2
                if lam > 1:
                    rx, ry = rx * math.sqrt(lam), ry * math.sqrt(lam)
                num = rx * rx * ry * ry - rx * rx * y1p * y1p - ry * ry * x1p * x1p
                den = rx * rx * y1p * y1p + ry * ry * x1p * x1p
                koef = math.sqrt(max(0, num / den)) if den else 0
                if velky == smer_:
                    koef = -koef
                cxp, cyp = koef * rx * y1p / ry, -koef * ry * x1p / rx
                ccx = cf * cxp - sf * cyp + (x + ex) / 2
                ccy = sf * cxp + cf * cyp + (y + ey) / 2
                t1 = math.atan2((y1p - cyp) / ry, (x1p - cxp) / rx)
                t2 = math.atan2((-y1p - cyp) / ry, (-x1p - cxp) / rx)
                dt = t2 - t1
                if smer_ and dt < 0:
                    dt += 2 * math.pi
                elif not smer_ and dt > 0:
                    dt -= 2 * math.pi
                for k in range(1, 25):
                    t = t1 + dt * k / 24
                    v.append((ccx + rx * math.cos(t) * cf - ry * math.sin(t) * sf,
                              ccy + rx * math.cos(t) * sf + ry * math.sin(t) * cf))
                    c.append(Path.LINETO)
            x, y = ex, ey
        else:
            i += 1
    if not v:
        return None
    return Path(v, c)


def svg_na_obrazek(svg, meritko=1.0, max_pixelu=14_000_000):
    """Vykreslí SVG vytvořené tímto programem do obrázku PIL (pro náhled v okně).
    Vrací (obrázek, skutečné měřítko)."""
    import numpy as np
    from PIL import Image
    from matplotlib.backends.backend_agg import RendererAgg
    from matplotlib.transforms import Affine2D
    koren = ET.fromstring(svg)
    W = float(koren.get("width"))
    H = float(koren.get("height"))
    if W * H * meritko * meritko > max_pixelu:
        meritko = math.sqrt(max_pixelu / (W * H))
    sw, sh = max(1, int(W * meritko)), max(1, int(H * meritko))
    rend = RendererAgg(sw, sh, 72)
    rend.clear()

    def styly(el):
        st = {}
        for k in ("fill", "stroke", "stroke-width", "stroke-dasharray", "stroke-linecap"):
            if el.get(k) is not None:
                st[k] = el.get(k)
        for cast in (el.get("style") or "").split(";"):
            if ":" in cast:
                k, hodn = cast.split(":", 1)
                st[k.strip()] = hodn.strip()
        return st

    def kresli(path, st, tx, ty):
        tr = Affine2D().translate(tx, ty).scale(meritko, -meritko).translate(0, sh)
        vypln = _barva(st.get("fill", "#000000"))
        obrys = _barva(st.get("stroke"))
        gc = rend.new_gc()
        gc.set_antialiased(True)
        gc.set_snap(False)                  # bez přichytávání k pixelům (tenké tahy písmen)
        path.should_simplify = False
        if obrys is not None:
            sirka_ = float(re.sub(r"[^\d.]", "", st.get("stroke-width", "1")) or 1) * meritko
            gc.set_linewidth(max(sirka_, 0.3))
            gc.set_foreground(obrys, isRGBA=True)
            gc.set_joinstyle("round")
            gc.set_capstyle("round" if st.get("stroke-linecap") == "round" else "butt")
            dash = st.get("stroke-dasharray")
            if dash and dash != "none":
                vals = [float(x) * meritko for x in re.split(r"[ ,]+", dash.strip()) if x]
                gc.set_dashes(0, [max(v_, 0.1) for v_ in vals])
        else:
            gc.set_linewidth(0)
        rend.draw_path(gc, path, tr, vypln)
        gc.restore()

    def projdi(el, tx, ty):
        tag = el.tag.split("}")[-1]
        if tag == "g":
            m = re.search(r"translate\(\s*([-\d.e]+)[ ,]+([-\d.e]+)", el.get("transform") or "")
            if m:
                tx, ty = tx + float(m.group(1)), ty + float(m.group(2))
            for ch in el:
                projdi(ch, tx, ty)
        elif tag == "path":
            p = _cesta_z_d(el.get("d", ""))
            if p is not None:
                kresli(p, styly(el), tx, ty)
        elif tag == "rect":
            x, y = float(el.get("x", 0)), float(el.get("y", 0))
            w, h = float(el.get("width", 0)), float(el.get("height", 0))
            r = min(float(el.get("rx", 0) or 0), w / 2, h / 2)
            if r > 0:
                d = (f"M{x + r},{y} H{x + w - r} Q{x + w},{y} {x + w},{y + r} V{y + h - r} "
                     f"Q{x + w},{y + h} {x + w - r},{y + h} H{x + r} Q{x},{y + h} {x},{y + h - r} "
                     f"V{y + r} Q{x},{y} {x + r},{y} Z")
            else:
                d = f"M{x},{y} H{x + w} V{y + h} H{x} Z"
            kresli(_cesta_z_d(d), styly(el), tx, ty)
        elif tag == "svg":
            for ch in el:
                projdi(ch, tx, ty)
    projdi(koren, 0.0, 0.0)
    pole = np.asarray(rend.buffer_rgba())
    obr = Image.fromarray(pole.copy(), "RGBA")
    pozadi_ = Image.new("RGBA", obr.size, (255, 255, 255, 255))
    pozadi_.alpha_composite(obr)
    return pozadi_.convert("RGB"), meritko


# =============================================================================
# 8) OKNO APLIKACE (Tkinter – součást Pythonu)
# =============================================================================
def otevri_soubor(cesta):
    try:
        if sys.platform.startswith("win"):
            os.startfile(cesta)                     # noqa: výchozí prohlížeč SVG
        elif sys.platform == "darwin":
            subprocess.Popen(["open", cesta])
        else:
            subprocess.Popen(["xdg-open", cesta])
    except Exception:
        import webbrowser
        webbrowser.open("file://" + os.path.abspath(cesta))


def spust_gui(pocatecni=None):
    import queue
    import tkinter as tk
    from tkinter import ttk, filedialog, messagebox

    ANO, NE = "☑", "☐"
    okno = tk.Tk()
    okno.title("KEGG dráha → SVG")
    ikona = os.path.join(ZDROJE, "assets", "ikona.ico")
    if os.path.exists(ikona):
        try:
            okno.iconbitmap(ikona)
        except tk.TclError:
            pass
    okno.geometry("1220x820")
    try:
        okno.tk.call("tk", "scaling", 1.25)
    except tk.TclError:
        pass
    fronta = queue.Queue()
    stav = dict(draha=None, otocene=set())

    # --- horní lišta: soubor ---
    horni = ttk.Frame(okno, padding=8)
    horni.pack(fill="x")
    ttk.Label(horni, text="KGML soubor:").pack(side="left")
    var_cesta = tk.StringVar(value=pocatecni or "")
    ttk.Entry(horni, textvariable=var_cesta, width=70).pack(side="left", padx=6, fill="x", expand=True)

    def vyber_soubor():
        start_dir = os.path.expanduser("~/Downloads")
        p = filedialog.askopenfilename(title="Vyber KGML (XML) soubor z KEGG",
                                       initialdir=start_dir if os.path.isdir(start_dir) else SLOZKA,
                                       filetypes=[("KGML / XML", "*.xml *.kgml"), ("Vše", "*.*")])
        if p:
            var_cesta.set(p)
            nacti()
    ttk.Button(horni, text="Otevřít…", command=vyber_soubor).pack(side="left")
    var_online = tk.BooleanVar(value=True)
    ttk.Checkbutton(horni, text="stahovat z KEGG", variable=var_online).pack(side="left", padx=8)
    tl_nacist = ttk.Button(horni, text="Načíst", command=lambda: nacti())
    tl_nacist.pack(side="left")

    # --- stavový řádek ---
    var_stav = tk.StringVar(value="Otevři KGML soubor (KEGG → stránka dráhy → Download KGML).")
    dolni = ttk.Frame(okno, padding=(8, 4))
    dolni.pack(side="bottom", fill="x")
    ttk.Label(dolni, textvariable=var_stav).pack(side="left")
    prubeh = ttk.Progressbar(dolni, mode="indeterminate", length=160)
    prubeh.pack(side="right")

    # --- střed: tabulky + volby ---
    stred = ttk.Panedwindow(okno, orient="horizontal")
    stred.pack(fill="both", expand=True, padx=8, pady=4)
    karty = ttk.Notebook(stred)
    stred.add(karty, weight=4)
    volby = ttk.Frame(stred, padding=10)
    stred.add(volby, weight=1)

    def tabulka(rodic, sloupce, sirky):
        ram = ttk.Frame(rodic)
        tv = ttk.Treeview(ram, columns=sloupce, show="headings", selectmode="extended")
        for s, w in zip(sloupce, sirky):
            tv.heading(s, text=s)
            tv.column(s, width=w, anchor="w", stretch=(w > 120))
        sb = ttk.Scrollbar(ram, orient="vertical", command=tv.yview)
        tv.configure(yscrollcommand=sb.set)
        tv.pack(side="left", fill="both", expand=True)
        sb.pack(side="right", fill="y")
        return ram, tv

    # karta Metabolity
    k_met = ttk.Frame(karty, padding=4)
    karty.add(k_met, text="Metabolity")
    ttk.Label(k_met, text="Klikni na řádek = zapnout/vypnout. Vypnutý metabolit uprostřed dráhy "
                          "se nahradí rámečkem „⋯ vynechány N metabolity“.\n"
                          "Přenašeče (ThPP, lipoyllysin …) a meziprodukty vázané na enzym jsou vypnuté: "
                          "přenašeč se nekreslí, vázané meziprodukty se sloučí do jedné šipky komplexu.",
              foreground="#555").pack(anchor="w", pady=(0, 4))
    ram, tv_met = tabulka(k_met, ("✓", "Název", "KEGG", "Struktura", "Napojení"), (36, 330, 80, 80, 90))
    ram.pack(fill="both", expand=True)

    # karta Reakce
    k_rx = ttk.Frame(karty, padding=4)
    karty.add(k_rx, text="Reakce")

    # karta Vzhled
    from tkinter import colorchooser
    k_vz = ttk.Frame(karty, padding=8)
    karty.add(k_vz, text="Vzhled")
    SOUBOR_VZHLEDU = os.path.join(SLOZKA, "nastaveni_vzhledu.json")
    vz_vars = {}
    vz_vzorky = {}
    levy = ttk.LabelFrame(k_vz, text="Barvy (klikni na čtvereček)", padding=8)
    levy.grid(row=0, column=0, sticky="nsew", padx=(0, 8))
    pravy = ttk.Frame(k_vz)
    pravy.grid(row=0, column=1, sticky="nsew")
    k_vz.columnconfigure(0, weight=1)
    k_vz.columnconfigure(1, weight=1)
    BARVY_GUI = [
        ("barva_pozadi", "pozadí obrázku"), ("barva_sipek", "šipky"), ("barva_nazvu", "názvy látek"),
        ("barva_genu", "geny u enzymů"), ("barva_ec", "EC čísla"), ("barva_vazeb", "vazby ve strukturách"),
        ("barva_nadpisu", "nadpis"), ("barva_popisu", "popisky, legenda"),
        ("barva_vychozi", "výchozí látka – výplň"), ("barva_vychozi_obrys", "výchozí látka – obrys"),
        ("barva_mapy", "navazující dráha – výplň"), ("barva_mapy_obrys", "navazující dráha – obrys"),
        ("barva_mapy_text", "navazující dráha – text"), ("barva_zkratky", "vynechané kroky – čára"),
        ("barva_mezery", "vynechané kroky – výplň"), ("barva_org", "enzymy organismu (zvýraznění)"),
    ]

    def vyber_barvu(klic):
        puv = vz_vars[klic].get()
        nova = colorchooser.askcolor(color=puv, title="Vyber barvu", parent=okno)[1]
        if nova:
            vz_vars[klic].set(nova)
            vz_vzorky[klic].configure(background=nova)

    for i, (klic, popis) in enumerate(BARVY_GUI):
        vz_vars[klic] = tk.StringVar(value=VYCHOZI_VZHLED[klic])
        vz = tk.Label(levy, width=3, relief="solid", borderwidth=1, background=VYCHOZI_VZHLED[klic],
                      cursor="hand2")
        vz.grid(row=i, column=0, padx=(0, 8), pady=2, sticky="w")
        vz.bind("<Button-1>", lambda e, k_=klic: vyber_barvu(k_))
        vz_vzorky[klic] = vz
        ttk.Label(levy, text=popis).grid(row=i, column=1, sticky="w")
    vz_vars["barevne_atomy"] = tk.BooleanVar(value=True)
    ttk.Checkbutton(levy, text="barevné atomy (O červeně, N modře …)",
                    variable=vz_vars["barevne_atomy"]).grid(row=len(BARVY_GUI), column=0, columnspan=2,
                                                            sticky="w", pady=(6, 0))

    ram_p = ttk.LabelFrame(pravy, text="Předvolba", padding=8)
    ram_p.pack(fill="x")
    var_predvolba = tk.StringVar(value="Klasický")
    ttk.Combobox(ram_p, textvariable=var_predvolba, values=list(PREDVOLBY), state="readonly",
                 width=20).pack(side="left")
    ttk.Button(ram_p, text="Použít", command=lambda: nastav_vzhled_gui(
        {**VYCHOZI_VZHLED, **PREDVOLBY[var_predvolba.get()]})).pack(side="left", padx=6)

    ram_t = ttk.LabelFrame(pravy, text="Tvary", padding=8)
    ram_t.pack(fill="x", pady=6)
    vz_vars["hrot"] = tk.StringVar(value="plny")
    radek_h = ttk.Frame(ram_t)
    radek_h.pack(anchor="w")
    ttk.Label(radek_h, text="Hrot šipky:").pack(side="left")
    for hodn, text in (("plny", "plný ➤"), ("trojuhelnik", "trojúhelník ▶"), ("otevreny", "otevřený >")):
        ttk.Radiobutton(radek_h, text=text, value=hodn, variable=vz_vars["hrot"]).pack(side="left", padx=4)

    def spin(rodic, klic, popis, od, do, krok):
        r_ = ttk.Frame(rodic)
        r_.pack(fill="x", pady=1)
        ttk.Label(r_, text=popis).pack(side="left")
        vz_vars[klic] = tk.DoubleVar(value=VYCHOZI_VZHLED[klic])
        ttk.Spinbox(r_, from_=od, to=do, increment=krok, textvariable=vz_vars[klic], width=6).pack(side="right")
    spin(ram_t, "tloustka", "Tloušťka šipek (px):", 0.5, 4, 0.1)
    spin(ram_t, "velikost_hrotu", "Velikost hrotu (px):", 5, 20, 1)
    spin(ram_t, "polomer_rohu", "Zaoblení rohů šipek (0 = ostré):", 0, 30, 1)
    spin(ram_t, "zaobleni_ramecku", "Zaoblení rámečků (0 = hranaté):", 0, 25, 1)

    ram_v = ttk.LabelFrame(pravy, text="Velikosti", padding=8)
    ram_v.pack(fill="x")
    spin(ram_v, "delka_vazby", "Velikost struktur (délka vazby, px):", 12, 40, 1)
    spin(ram_v, "velikost_atomu", "Písmo atomů (C, O, N …):", 0.3, 0.9, 0.05)
    spin(ram_v, "pismo_nazvu", "Písmo názvů látek (pt):", 8, 24, 0.5)
    spin(ram_v, "pismo_enzymu", "Písmo enzymů (pt):", 6, 18, 0.5)
    spin(ram_v, "rozestup", "Rozestupy (× výchozí):", 0.6, 2.0, 0.1)

    ram_f = ttk.LabelFrame(pravy, text="Písmo", padding=8)
    ram_f.pack(fill="x", pady=6)
    try:
        from matplotlib import font_manager
        fonty = sorted({f.name for f in font_manager.fontManager.ttflist})
    except Exception:
        fonty = ["DejaVu Sans"]
    vz_vars["font"] = tk.StringVar(value=VYCHOZI_VZHLED["font"])
    ttk.Combobox(ram_f, textvariable=vz_vars["font"], values=fonty, width=30).pack(anchor="w")

    ram_b = ttk.Frame(pravy)
    ram_b.pack(fill="x", pady=4)

    def ziskej_vzhled():
        v = {}
        for k_, var in vz_vars.items():
            try:
                v[k_] = var.get()
            except (tk.TclError, ValueError):
                v[k_] = VYCHOZI_VZHLED[k_]
        return v

    def nastav_vzhled_gui(v):
        for k_, var in vz_vars.items():
            if k_ in v:
                var.set(v[k_])
                if k_ in vz_vzorky:
                    vz_vzorky[k_].configure(background=v[k_])

    def uloz_vzhled():
        import json
        with open(SOUBOR_VZHLEDU, "w", encoding="utf-8") as f_:
            json.dump(ziskej_vzhled(), f_, ensure_ascii=False, indent=2)
        var_stav.set(f"Vzhled uložen do {SOUBOR_VZHLEDU} – načte se při příštím spuštění.")

    ttk.Button(ram_b, text="Uložit jako výchozí", command=uloz_vzhled).pack(side="left")
    ttk.Button(ram_b, text="Obnovit původní", command=lambda: nastav_vzhled_gui(VYCHOZI_VZHLED)).pack(side="left", padx=4)
    ttk.Button(ram_b, text="Náhled", command=lambda: nahled()).pack(side="left")
    if os.path.exists(SOUBOR_VZHLEDU):
        try:
            nastav_vzhled_gui(nacti_vzhled_ze_souboru(SOUBOR_VZHLEDU))
        except Exception:
            pass

    # karta Náhled
    k_nah = ttk.Frame(karty, padding=4)
    karty.add(k_nah, text="Náhled")
    lista = ttk.Frame(k_nah)
    lista.pack(fill="x", pady=(0, 4))
    ttk.Button(lista, text="Obnovit náhled", command=lambda: nahled()).pack(side="left")
    ttk.Label(lista, text="   Zvětšení:").pack(side="left")
    var_zoom = tk.StringVar(value="přizpůsobit šířce")
    cb_zoom = ttk.Combobox(lista, textvariable=var_zoom, state="readonly", width=17,
                           values=["přizpůsobit šířce", "25 %", "50 %", "75 %", "100 %", "150 %", "200 %"])
    cb_zoom.pack(side="left", padx=4)
    var_nah_info = tk.StringVar(value="Klikni na „Náhled“ (vpravo) – obrázek se vykreslí sem, nic se neukládá.")
    ttk.Label(lista, textvariable=var_nah_info, foreground="#555").pack(side="left", padx=8)
    ram_c = ttk.Frame(k_nah)
    ram_c.pack(fill="both", expand=True)
    platno = tk.Canvas(ram_c, background="#ffffff", highlightthickness=0)
    sb_y = ttk.Scrollbar(ram_c, orient="vertical", command=platno.yview)
    sb_x = ttk.Scrollbar(ram_c, orient="horizontal", command=platno.xview)
    platno.configure(yscrollcommand=sb_y.set, xscrollcommand=sb_x.set)
    sb_y.pack(side="right", fill="y")
    sb_x.pack(side="bottom", fill="x")
    platno.pack(side="left", fill="both", expand=True)
    platno.bind("<ButtonPress-1>", lambda e: platno.scan_mark(e.x, e.y))
    platno.bind("<B1-Motion>", lambda e: platno.scan_dragto(e.x, e.y, gain=1))
    platno.bind("<MouseWheel>", lambda e: platno.yview_scroll(int(-e.delta / 120) or (-1 if e.delta > 0 else 1), "units"))
    platno.bind("<Shift-MouseWheel>", lambda e: platno.xview_scroll(int(-e.delta / 120) or (-1 if e.delta > 0 else 1), "units"))
    platno.bind("<Button-4>", lambda e: platno.yview_scroll(-3, "units"))
    platno.bind("<Button-5>", lambda e: platno.yview_scroll(3, "units"))
    ttk.Label(k_rx, text="Vypnutá reakce se nekreslí. „Otočit směr“ obrátí šipku nevratné reakce "
                         "(KGML ji někdy uvádí obráceně).", foreground="#555").pack(anchor="w", pady=(0, 4))
    ram, tv_rx = tabulka(k_rx, ("✓", "Z", "Typ", "Do", "EC", "Geny organismu", "KEGG"),
                         (36, 220, 50, 220, 170, 120, 110))
    ram.pack(fill="both", expand=True)

    def prepni(tv, ids=None, hodnota=None):
        for iid in (ids or tv.selection()):
            v = list(tv.item(iid, "values"))
            v[0] = hodnota if hodnota is not None else (NE if v[0] == ANO else ANO)
            tv.item(iid, values=v)
        aktualizuj_pocty()

    def klik(event, tv):
        iid = tv.identify_row(event.y)
        if iid and tv.identify_region(event.x, event.y) == "cell":
            if not (event.state & 0x0005):          # bez Shift/Ctrl = přepnout
                prepni(tv, [iid])
    tv_met.bind("<ButtonRelease-1>", lambda e: klik(e, tv_met))
    tv_rx.bind("<ButtonRelease-1>", lambda e: klik(e, tv_rx))
    tv_met.bind("<space>", lambda e: prepni(tv_met))
    tv_rx.bind("<space>", lambda e: prepni(tv_rx))

    def tlacitka(rodic, tv, extra=()):
        r = ttk.Frame(rodic)
        r.pack(fill="x", pady=4)
        ttk.Button(r, text="Vybrat vše", command=lambda: prepni(tv, tv.get_children(), ANO)).pack(side="left")
        ttk.Button(r, text="Zrušit vše", command=lambda: prepni(tv, tv.get_children(), NE)).pack(side="left", padx=4)
        ttk.Button(r, text="Obrátit výběr", command=lambda: prepni(tv, tv.get_children())).pack(side="left")
        ttk.Label(r, text="  (Ctrl/Shift + klik = označit více řádků, mezerník = přepnout označené)",
                  foreground="#777").pack(side="left")
        for text, fn in extra:
            ttk.Button(r, text=text, command=fn).pack(side="right", padx=2)
    tlacitka(k_met, tv_met)

    def otoc():
        d = stav["draha"]
        for iid in tv_rx.selection():
            i = int(iid)
            stav["otocene"] ^= {i}
            r = d.reakce[i]
            a, b = (r["b"], r["a"]) if i in stav["otocene"] else (r["a"], r["b"])
            v = list(tv_rx.item(iid, "values"))
            v[1], v[3] = d.metabolity[a]["nazev"], d.metabolity[b]["nazev"]
            tv_rx.item(iid, values=v)
    tlacitka(k_rx, tv_rx, extra=[("Otočit směr označených", otoc)])

    # --- volby ---
    ttk.Label(volby, text="Nastavení obrázku", font=("", 11, "bold")).pack(anchor="w", pady=(0, 8))
    ttk.Label(volby, text="Výchozí metabolit (vlevo / nahoře):").pack(anchor="w")
    var_start = tk.StringVar()
    cb_start = ttk.Combobox(volby, textvariable=var_start, state="readonly", width=34)
    cb_start.pack(anchor="w", pady=(0, 8), fill="x")
    ttk.Label(volby, text="Název metabolitu:").pack(anchor="w")
    var_pop = tk.StringVar(value="pod")
    ttk.Radiobutton(volby, text="pod strukturou", value="pod", variable=var_pop).pack(anchor="w")
    ttk.Radiobutton(volby, text="vedle struktury", value="vpravo", variable=var_pop).pack(anchor="w", pady=(0, 6))
    ttk.Label(volby, text="Rozložení:").pack(anchor="w")
    var_rozl = tk.StringVar(value="radky")
    radek_r = ttk.Frame(volby)
    radek_r.pack(anchor="w", pady=(0, 2))
    ttk.Radiobutton(radek_r, text="řádkové", value="radky", variable=var_rozl).pack(side="left")
    ttk.Radiobutton(radek_r, text="kruhové (cyklus)", value="kruh", variable=var_rozl).pack(side="left", padx=8)
    var_cyklus = tk.StringVar(value="")
    ttk.Label(volby, textvariable=var_cyklus, foreground="#7b3fa0", wraplength=250).pack(anchor="w", pady=(0, 6))
    var_odkazy = tk.BooleanVar(value=True)
    var_geny = tk.BooleanVar(value=True)
    var_zvyr = tk.BooleanVar(value=False)
    var_zkr = tk.BooleanVar(value=True)
    var_otevrit = tk.BooleanVar(value=True)
    var_glyk = tk.BooleanVar(value=False)
    for text, var in (("navazující dráhy (rámečky)", var_odkazy),
                      ("názvy genů organismu u enzymů", var_geny),
                      ("zvýraznit enzymy organismu zeleně", var_zvyr),
                      ("šipky přes vynechané kroky", var_zkr),
                      ("glykany jako atomové vzorce", var_glyk)):
        ttk.Checkbutton(volby, text=text, variable=var).pack(anchor="w", pady=1)
    var_coa = tk.BooleanVar(value=True)
    ttk.Checkbutton(volby, text="zkracovat koenzym A na „CoA“", variable=var_coa).pack(anchor="w", pady=1)
    ttk.Checkbutton(volby, text="otevřít SVG po vytvoření", variable=var_otevrit).pack(anchor="w", pady=1)
    ttk.Label(volby, text="Vlastní nadpis (prázdné = podle dráhy):").pack(anchor="w", pady=(8, 0))
    var_nadpis = tk.StringVar()
    ttk.Entry(volby, textvariable=var_nadpis).pack(fill="x")
    ttk.Separator(volby).pack(fill="x", pady=10)
    ttk.Button(volby, text="Které dráhy to umí?",
               command=lambda: messagebox.showinfo("Podporované dráhy", PODPOROVANE_DRAHY)).pack(fill="x", pady=2)
    ttk.Button(volby, text="Jen reakce organismu",
               command=lambda: jen_organismus()).pack(fill="x", pady=2)
    var_pocty = tk.StringVar(value="")
    ttk.Label(volby, textvariable=var_pocty, foreground="#555", wraplength=230).pack(anchor="w", pady=6)
    radek_t = ttk.Frame(volby)
    radek_t.pack(fill="x", pady=(10, 2))
    tl_nahled = ttk.Button(radek_t, text="Náhled", command=lambda: nahled())
    tl_nahled.pack(side="left", fill="x", expand=True, ipady=8, padx=(0, 4))
    tl_svg = ttk.Button(radek_t, text="Vytvořit SVG…", command=lambda: vytvor())
    tl_svg.pack(side="left", fill="x", expand=True, ipady=8)
    var_pdf = tk.BooleanVar(value=False)
    ttk.Checkbutton(volby, text="navíc PDF", variable=var_pdf).pack(anchor="w")

    def vybrane():
        met = {tv_met.item(i, "values")[2] for i in tv_met.get_children()
               if tv_met.item(i, "values")[0] == ANO}
        rx = {int(i) for i in tv_rx.get_children() if tv_rx.item(i, "values")[0] == ANO}
        return met, rx

    def aktualizuj_pocty():
        if not stav["draha"]:
            return
        met, rx = vybrane()
        var_pocty.set(f"Vybráno {len(met)} z {len(tv_met.get_children())} metabolitů "
                      f"a {len(rx)} z {len(tv_rx.get_children())} reakčních kroků.")

    def jen_organismus():
        d = stav["draha"]
        if not d:
            return
        rx = {i for i, r in enumerate(d.reakce) if r["geny"]}
        if not rx:
            messagebox.showinfo("Jen organismus", "Tento KGML neobsahuje geny organismu.")
            return
        met = {c for i in rx for c in (d.reakce[i]["a"], d.reakce[i]["b"])
               if d.metabolity[c].get("role", "normal") == "normal"}
        for iid in tv_rx.get_children():
            prepni(tv_rx, [iid], ANO if int(iid) in rx else NE)
        for iid in tv_met.get_children():
            prepni(tv_met, [iid], ANO if tv_met.item(iid, "values")[2] in met else NE)

    # --- práce na pozadí ---
    def pracuj(fn, hotovo):
        prubeh.start(12)
        tl_nacist.state(["disabled"])
        tl_svg.state(["disabled"])
        tl_nahled.state(["disabled"])

        def beh():
            try:
                fronta.put(("ok", hotovo, fn()))
            except BaseException as e:  # jakákoli chyba se zobrazí v okně (nikdy tiché zaseknutí)
                import traceback
                traceback.print_exc()
                fronta.put(("chyba", None, e))
        threading.Thread(target=beh, daemon=True).start()

    def kontrola():
        try:
            while True:
                typ, hotovo, data = fronta.get_nowait()
                if typ == "log":
                    var_stav.set(data)
                    continue
                prubeh.stop()
                tl_nacist.state(["!disabled"])
                tl_svg.state(["!disabled"])
                tl_nahled.state(["!disabled"])
                if typ == "chyba":
                    var_stav.set(f"Chyba: {data}")
                    messagebox.showerror("Chyba", str(data))
                else:
                    hotovo(data)
        except queue.Empty:
            pass
        okno.after(100, kontrola)
    okno.after(100, kontrola)

    def nacti():
        p = var_cesta.get().strip()
        if p and not os.path.exists(p):
            messagebox.showerror("Soubor nenalezen", p)
            return
        online = var_online.get()
        pracuj(lambda: nacti_drahu(p or None, online, log=lambda t: fronta.put(("log", None, t))),
               naplnit)

    def naplnit(d):
        stav["draha"], stav["otocene"] = d, set()
        tv_met.delete(*tv_met.get_children())
        tv_rx.delete(*tv_rx.get_children())
        stupen = defaultdict(int)
        for r in d.reakce:
            stupen[r["a"]] += 1
            stupen[r["b"]] += 1
        for c, m in sorted(d.metabolity.items(), key=lambda kv: kv[1]["poradi"]):
            role = m.get("role", "normal")
            pozn = {"prenasec": " · přenašeč", "vazany": " · vázaný na enzym"}.get(role, "")
            tv_met.insert("", "end", iid="m" + c,
                          values=(ANO if role == "normal" else NE, m["nazev"], c,
                                  "ano" if m["zdroj"] else "—", f"{stupen[c]} reakcí{pozn}"))
        for i, r in enumerate(d.reakce):
            typ = {"<=>": "⇌", "->": "→", "->+": "→ (+)"}[r["typ"]]
            tv_rx.insert("", "end", iid=str(i),
                         values=(ANO, d.metabolity[r["a"]]["nazev"], typ, d.metabolity[r["b"]]["nazev"],
                                 " · ".join(r["ec"]) or "—", " · ".join(r["geny"]), " ".join(r["rn"])))
        moznosti = sorted(d.metabolity, key=lambda c: d.metabolity[c]["nazev"].lower())
        cb_start["values"] = [f'{d.metabolity[c]["nazev"]}  [{c}]' for c in moznosti]
        if d.start:
            var_start.set(f'{d.metabolity[d.start]["nazev"]}  [{d.start}]')
        bez = sum(1 for m in d.metabolity.values() if not m["zdroj"])
        var_stav.set(f"{d.org}{d.cislo} – {d.titul}: {len(d.metabolity)} metabolitů, {len(d.reakce)} "
                     f"reakčních kroků ({d.zdroj_dat})" + (f"; {bez} bez struktury" if bez else ""))
        aktualizuj_pocty()
        stav["nahled_svg"] = None
        cyk = cyklus_v_draze(d)
        if len(cyk) >= 4:
            var_cyklus.set(f"↻ Dráha obsahuje cyklus ({len(cyk)} látek: "
                           f"{d.metabolity[cyk[0]]['nazev']} → … ) – zkus kruhové rozložení.")
        else:
            var_cyklus.set("")
        if len(d.metabolity) > VELKA_MAPA:
            messagebox.showinfo(
                "Velká mapa",
                f"Dráha má {len(d.metabolity)} metabolitů. Pro čitelný obrázek doporučuji "
                f"na kartách Metabolity/Reakce vypnout části, které nepotřebuješ.")

    def parametry():
        m = re.search(r"\[([A-Z]\d{5})\]", var_start.get())
        return dict(start=m.group(1) if m else None, otocene=set(stav["otocene"]),
                    popisky=var_pop.get(), odkazy=var_odkazy.get(), geny=var_geny.get(),
                    zvyraznit=var_zvyr.get(), zkratky=var_zkr.get(),
                    glykany_atomy=var_glyk.get(), zkratit_coa=var_coa.get(),
                    vzhled=ziskej_vzhled(),
                    nadpis=var_nadpis.get(), rozlozeni=var_rozl.get())

    def meritko_nahledu(svg_text):
        z = var_zoom.get()
        if z.endswith("%"):
            return float(z.split()[0]) / 100
        sirka_svg = float(re.search(r'width="([\d.]+)"', svg_text).group(1))
        return max(0.05, (platno.winfo_width() - 4) / sirka_svg)

    def zobraz_nahled(v):
        from PIL import ImageTk
        obr, mer, n, kriz = v
        stav["nahled_obr"] = ImageTk.PhotoImage(obr)      # reference musí zůstat
        platno.delete("all")
        platno.create_image(0, 0, anchor="nw", image=stav["nahled_obr"])
        platno.configure(scrollregion=(0, 0, obr.width, obr.height))
        var_nah_info.set(f"{n} uzlů · křížení šipek: {kriz} · zobrazeno {mer * 100:.0f} % "
                         f"(SVG je vektorové, při uložení je ostré v jakékoli velikosti)")
        var_stav.set("Náhled hotový – pokud se líbí, ulož ho tlačítkem „Vytvořit SVG…“.")

    def nahled():
        d = stav["draha"]
        if not d:
            messagebox.showinfo("Nejdřív načti dráhu", "Otevři a načti KGML soubor.")
            return
        met, rx = vybrane()
        if not met:
            messagebox.showwarning("Prázdný výběr", "Není vybrán žádný metabolit.")
            return
        karty.select(k_nah)
        okno.update_idletasks()
        param = parametry()

        def prace():
            fronta.put(("log", None, "Počítám náhled…"))
            n, kriz = vytvor_svg(d, met, rx, None, **param)
            svg_text = POSLEDNI_SVG
            stav["nahled_svg"] = (svg_text, n, kriz)
            obr, mer = svg_na_obrazek(svg_text, meritko_nahledu(svg_text))
            return obr, mer, n, kriz
        pracuj(prace, zobraz_nahled)

    def zmena_zoomu(_e=None):
        if not stav.get("nahled_svg"):
            return
        svg_text, n, kriz = stav["nahled_svg"]
        mer = meritko_nahledu(svg_text)
        pracuj(lambda: (*svg_na_obrazek(svg_text, mer), n, kriz), zobraz_nahled)
    cb_zoom.bind("<<ComboboxSelected>>", zmena_zoomu)

    def vytvor():
        d = stav["draha"]
        if not d:
            messagebox.showinfo("Nejdřív načti dráhu", "Otevři a načti KGML soubor.")
            return
        met, rx = vybrane()
        if not met:
            messagebox.showwarning("Prázdný výběr", "Není vybrán žádný metabolit.")
            return
        if len(met) > VELKA_MAPA and not messagebox.askyesno(
                "Velká mapa",
                f"Vybráno {len(met)} metabolitů. Výpočet může trvat i několik minut a obrázek "
                f"bude hodně velký a hůř čitelný.\n\nDoporučuji vypnout části, které nepotřebuješ "
                f"(vynechané úseky se nahradí rámečkem „⋯ vynechány N metabolity“).\n\n"
                f"Pokračovat i tak?"):
            return
        cil = filedialog.asksaveasfilename(
            title="Uložit SVG", defaultextension=".svg", initialdir=SLOZKA,
            initialfile=f"{d.org}{d.cislo}_draha.svg", filetypes=[("SVG", "*.svg")])
        if not cil:
            return
        param = parametry()
        pdf = var_pdf.get()

        def prace():
            fronta.put(("log", None, "Kreslím struktury a počítám rozložení…"))
            n, kriz = vytvor_svg(d, met, rx, cil, **param)
            if pdf:
                preved(cil, "pdf")
            return cil, n, kriz

        def hotovo(v):
            cil_, n, kriz = v
            var_stav.set(f"Uloženo: {cil_}  ({n} uzlů, křížení šipek: {kriz})")
            if var_otevrit.get():
                otevri_soubor(cil_)
        pracuj(prace, hotovo)

    if pocatecni:
        okno.after(200, nacti)
    okno.mainloop()


# =============================================================================
# 9) PŘÍKAZOVÁ ŘÁDKA
# =============================================================================
def najdi_kgml(nazev):
    import glob
    if nazev and os.path.exists(nazev):
        return nazev
    kandidati = [os.path.join(SLOZKA, nazev)] if nazev else []
    cislo = re.sub(r"\D", "", os.path.basename(nazev or "")) or "00360"
    for d in (SLOZKA, os.getcwd(), os.path.expanduser("~/Downloads")):
        kandidati += sorted(glob.glob(os.path.join(d, f"*{cislo}*.xml")))
    for k in kandidati:
        if os.path.exists(k):
            return k
    return None


def preved(svg, fmt):
    cil = svg.rsplit(".", 1)[0] + "." + fmt
    if fmt == "emf":
        if shutil.which("inkscape"):
            subprocess.run(["inkscape", svg, f"--export-filename={cil}"], check=True)
            return cil
        print("EMF: Inkscape nenalezen – otevřete SVG v Inkscapu/Corelu a uložte jako EMF.")
        return None
    if shutil.which("rsvg-convert"):
        opt = ["-z", "1.5"] if fmt == "png" else []
        subprocess.run(["rsvg-convert", "-f", fmt, *opt, "-o", cil, svg], check=True)
        return cil
    try:
        import cairosvg
        getattr(cairosvg, f"svg2{fmt}")(url=svg, write_to=cil, **({"scale": 1.5} if fmt == "png" else {}))
        return cil
    except ImportError:
        print(f"{fmt.upper()}: chybí rsvg-convert i cairosvg (pip install cairosvg).")
        return None


def main():
    if len(sys.argv) == 1:          # bez argumentů (např. ▷ v PyCharmu) -> okno
        spust_gui(najdi_kgml("gga00360.xml"))
        return
    ap = argparse.ArgumentParser(description="KEGG dráha (KGML) -> vektorové SVG se strukturami.")
    ap.add_argument("kgml", nargs="?", default="gga00360.xml", help="KGML soubor z KEGG")
    ap.add_argument("-o", "--vystup", default=None, help="výstupní SVG")
    ap.add_argument("--gui", action="store_true", help="otevřít okno s výběrem")
    ap.add_argument("--offline", action="store_true", help="nestahovat nic z KEGG")
    ap.add_argument("--vynechat", default="", help="C-čísla metabolitů oddělená čárkou")
    ap.add_argument("--jen-organismus", action="store_true", help="jen reakce s geny organismu")
    ap.add_argument("--start", default=None, help="C-číslo výchozího metabolitu")
    ap.add_argument("--popisky", choices=["pod", "vpravo"], default="pod")
    ap.add_argument("--bez-odkazu", action="store_true", help="bez rámečků navazujících drah")
    ap.add_argument("--bez-genu", action="store_true", help="bez názvů genů u enzymů")
    ap.add_argument("--zvyraznit", action="store_true", help="enzymy organismu zeleně")
    ap.add_argument("--kruh", action="store_true", help="kruhové rozložení (cyklické dráhy)")
    ap.add_argument("--predvolba", choices=list(PREDVOLBY), default=None, help="předvolba vzhledu")
    ap.add_argument("--vzhled", default=None,
                    help="soubor s nastavením vzhledu (JSON, např. nastaveni_vzhledu.json z okna)")
    ap.add_argument("--s-kofaktory", action="store_true",
                    help="kreslit i přenašeče (ThPP, lipoyllysin …) a meziprodukty vázané na enzym")
    ap.add_argument("--plne-coa", action="store_true", help="kreslit celý koenzym A (bez zkratky CoA)")
    ap.add_argument("--delka-vazby", type=float, default=None, help="velikost struktur (px, výchozí 20)")
    ap.add_argument("--velikost-atomu", type=float, default=None, help="písmo atomů (výchozí 0.5)")
    ap.add_argument("--nadpis", default=None, help="vlastní nadpis obrázku")
    ap.add_argument("--podporovane", action="store_true", help="vypsat, které dráhy program umí")
    ap.add_argument("--glykany-atomy", action="store_true",
                    help="glykany jako plné atomové vzorce místo symbolů SNFG")
    ap.add_argument("--png", action="store_true")
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--emf", action="store_true", help="EMF metafile (vyžaduje Inkscape)")
    args = ap.parse_args()
    if args.podporovane:
        print(PODPOROVANE_DRAHY)
        return

    kgml = najdi_kgml(args.kgml)
    if args.gui:
        spust_gui(kgml)
        return
    if kgml is None:
        if re.sub(r"\D", "", args.kgml) in ("00360", ""):
            print("KGML nenalezen – používám vestavěnou mapu 00360 (bez genů organismu).")
        else:
            sys.exit(f"Soubor {args.kgml} nenalezen.")
    d = nacti_drahu(kgml, online=not args.offline)
    rx = set(range(len(d.reakce)))
    if args.jen_organismus:
        rx = {i for i in rx if d.reakce[i]["geny"]}
    met = {c for i in rx for c in (d.reakce[i]["a"], d.reakce[i]["b"])
           if args.s_kofaktory or d.metabolity[c].get("role", "normal") == "normal"}
    met -= {x.strip() for x in args.vynechat.split(",") if x.strip()}
    vystup = args.vystup or os.path.join(SLOZKA, f"{d.org}{d.cislo}_draha.svg")
    vzhled = dict(PREDVOLBY.get(args.predvolba, {})) if args.predvolba else {}
    soubor_vz = args.vzhled or os.path.join(SLOZKA, "nastaveni_vzhledu.json")
    if os.path.exists(soubor_vz) and not args.predvolba:
        vzhled.update(nacti_vzhled_ze_souboru(soubor_vz))
    if args.delka_vazby:
        vzhled["delka_vazby"] = args.delka_vazby
    if args.velikost_atomu:
        vzhled["velikost_atomu"] = args.velikost_atomu
    n, kriz = vytvor_svg(d, met, rx, vystup, start=args.start, popisky=args.popisky,
                         odkazy=not args.bez_odkazu, geny=not args.bez_genu,
                         zvyraznit=args.zvyraznit, glykany_atomy=args.glykany_atomy,
                         zkratit_coa=not args.plne_coa, delka_vazby=args.delka_vazby,
                         velikost_atomu=args.velikost_atomu, nadpis=args.nadpis,
                         rozlozeni="kruh" if args.kruh else "radky", vzhled=vzhled)
    print(f"Uloženo: {vystup}   ({n} uzlů, křížení šipek: {kriz})")
    for fmt in ("png", "pdf", "emf"):
        if getattr(args, fmt):
            c = preved(vystup, fmt)
            if c:
                print("Uloženo:", c)


if __name__ == "__main__":
    main()
