# KEGG_Pathway – metabolické dráhy z KEGG jako vektorové SVG se strukturami

Program vezme dráhu z databáze [KEGG](https://www.kegg.jp/kegg/pathway.html) (soubor KGML)
a nakreslí ji znovu, tentokrát se **strukturními vzorci** místo koleček. Výsledkem je čisté
**vektorové SVG**: struktury, šipky i veškerý text jsou křivky, takže obrázek zůstane ostrý
při libovolném zvětšení a dá se dál upravovat v Inkscapu, Illustratoru nebo Corelu.

![Ukázka: metabolismus fenylalaninu (detail)](docs/ukazka_detail.png)

Celá mapa: [priklady/gga00360_fenylalanin.svg](priklady/gga00360_fenylalanin.svg) ·
s vynechanými kroky: [priklady/gga00360_vynechane_kroky.svg](priklady/gga00360_vynechane_kroky.svg) ·
[náhled celé mapy (PNG)](docs/ukazka_celkem.png)

*English summary: draws any KEGG metabolic pathway (KGML file) as a clean, fully vector SVG
with chemical structures (RDKit), glycan SNFG symbols, orthogonal non-crossing arrows,
EC numbers / organism genes and links to neighbouring pathways. Includes a Tkinter GUI for
choosing which metabolites and reactions to draw; hidden intermediates are replaced by a
"steps omitted" box. UI and output are in Czech.*

---

## Co program umí

- **Strukturní vzorce** všech metabolitů (RDKit), všechny ve stejném měřítku, se zmenšenými
  atomovými symboly. Koenzym A se zkracuje na „CoA“ (dá se vypnout).
- **Glykany** (např. gangliosidy v dráze 00604) se kreslí standardními symboly
  [SNFG](https://www.ncbi.nlm.nih.gov/glycans/snfg.html) s popisky vazeb (β1-4, α2-3 …),
  nebo volitelně jako plné atomové vzorce.
- **Rozložení bez křížení:** vrstvový (Sugiyamův) algoritmus s minimalizací křížení,
  pravoúhlé šipky, každá s vlastní svislou „kolejí“, a popisky enzymů s kontrolou kolizí.
- **Správný směr šipek:** nevratné reakce mají plný hrot (→), vratné dvojici polovičních (⇌).
  Reakce bez známého enzymu jsou čárkovaně.
- **Enzymy:** u šipek jsou EC čísla a názvy genů organismu z KGML (např. kuře `gga`, člověk `hsa`).
  Volitelně se enzymy organismu zvýrazní zeleně.
- **Navazující dráhy** jako zaoblené rámečky (citrátový cyklus, metabolismus tyrosinu …).
- **Okno s výběrem:** zaškrtneš, které metabolity a reakce chceš. Když vypneš metabolit
  uprostřed dráhy, nahradí ho fialový rámeček **„⋯ vynechány 2 metabolity: …“**
  s tečkovanými šipkami, takže je vidět, že kroky chybí.
- Export do **SVG**, volitelně **PDF/PNG** (přes `rsvg-convert` nebo `cairosvg`)
  a **EMF** (Windows metafile, přes Inkscape).

## Instalace

Potřebuješ Python 3.9 nebo novější.

```bash
git clone https://github.com/Zacker452/KEGG_Pathway.git
cd KEGG_Pathway
pip install -r requirements.txt
```

Tkinter (okno aplikace) je součástí Pythonu na Windows a macOS. Na Linuxu ho případně
doinstaluj (`sudo apt install python3-tk`).

## Kde vzít KGML soubor

1. Otevři dráhu v KEGG, např. <https://www.kegg.jp/pathway/gga00360> (fenylalanin u kuřete).
   Kód organismu (`gga`, `hsa`, `mmu`, `eco` …) i číslo dráhy můžeš změnit přímo v adrese.
2. Klikni na **Download KGML** a soubor ulož (např. do složky Stažené soubory).

Pro dráhu **00360 (fenylalanin)** má program data vestavěná, takže funguje i bez KGML
a bez internetu. Jen pak nezná geny organismu.

## Použití

### Okno aplikace

```bash
python kegg_draha_svg.py
```

(V PyCharmu stačí soubor spustit zeleným ▷.)

![Okno aplikace](docs/okno.png)

1. **Otevřít…** – vyber KGML soubor. Program sám hledá i ve složce Stažené soubory.
2. Karta **Metabolity** – klikem na řádek látku zapneš nebo vypneš (☑/☐).
   Ctrl/Shift + klik označí více řádků, mezerník je pak přepne najednou.
3. Karta **Reakce** – vypnutá reakce se nekreslí. **Otočit směr označených**
   obrátí šipku (KGML občas uvádí nevratnou reakci obráceně).
4. Vpravo nastav vzhled a klikni na **Vytvořit SVG…**.

### Příkazová řádka

```bash
python kegg_draha_svg.py gga00360.xml                       # celá dráha -> gga00360_draha.svg
python kegg_draha_svg.py gga00360.xml --vynechat C00601,C07086 --pdf
python kegg_draha_svg.py hsa00010.xml --jen-organismus --zvyraznit
python kegg_draha_svg.py --podporovane                      # které dráhy program umí
python kegg_draha_svg.py --help                             # všechny volby
```

## Nastavení

| V okně | Na příkazové řádce | Co dělá |
|---|---|---|
| Výchozí metabolit | `--start C00079` | látka úplně vlevo, od ní se dráha rozkládá |
| Název pod / vedle struktury | `--popisky pod\|vpravo` | umístění názvu metabolitu |
| navazující dráhy | `--bez-odkazu` | rámečky s odkazy na další dráhy KEGG |
| názvy genů organismu | `--bez-genu` | geny z KGML u šipek |
| zvýraznit enzymy organismu zeleně | `--zvyraznit` | jako zelené boxy v KEGG |
| šipky přes vynechané kroky | – | rámeček „⋯ vynechány N metabolity“ místo skrytých látek |
| glykany jako atomové vzorce | `--glykany-atomy` | místo symbolů SNFG plná struktura |
| zkracovat koenzym A na „CoA“ | `--plne-coa` (vypne) | celý CoA má 48 atomů, se zkratkou je struktura přehlednější |
| Velikost struktur | `--delka-vazby 20` | délka vazby v px (12–40) |
| Velikost písma atomů | `--velikost-atomu 0.5` | relativní velikost symbolů C, O, N … |
| Vlastní nadpis | `--nadpis "…"` | jinak se použije název dráhy |
| Jen reakce organismu | `--jen-organismus` | jen reakce, pro které má organismus gen |
| navíc PDF | `--pdf`, `--png`, `--emf` | další formáty |
| stahovat z KEGG | `--offline` (vypne) | bez internetu se použije jen mezipaměť |

Další konstanty vzhledu (barvy, mezery, písmo, poloměr rohů šipek) jsou na začátku sekce
`VZHLED` ve skriptu.

## Které dráhy jsou podporované

| Typ dráhy | Podpora |
|---|---|
| Metabolické dráhy **00xxx** (glykolýza 00010, citrátový cyklus 00020, puriny 00230, aminokyseliny …) | ✔ |
| Libovolný organismus (`gga`, `hsa`, `mmu`, `eco` …) i referenční mapy (`map`, `ko`) | ✔ |
| Glykanové dráhy (00510, 00604 …) – symboly SNFG | ✔ |
| Signální a regulační dráhy **04xxx**, dráhy nemocí **05xxx** | ✘ – v KGML nejsou reakce mezi látkami, jen vztahy proteinů (program to oznámí) |
| Přehledové mapy **01100, 01200** … | ✘ – příliš velké, bez popsaných reakcí |

**Na co si dát pozor**

- Cyklické dráhy (citrátový cyklus) se kreslí zleva doprava jako řada s šipkou zpět, ne jako kruh.
- Látky bez struktury v KEGG (protein, akceptor, polymery) mají jen rámeček s názvem.
- Uzlové metabolity (pyruvát, acetyl-CoA), které KEGG kreslí na více místech, jsou v obrázku
  jednou. Vede z nich proto víc šipek.
- U map nad 120 látek program upozorní. Doporučuji vypnout části, které nepotřebuješ.
- České názvy látek jsou vestavěné pro dráhu 00360. U ostatních drah se použijí anglické
  názvy z KEGG (vlastní překlady lze doplnit do tabulky `CESKE_NAZVY`).

## Jak to funguje

1. **Načtení dat.** Z KGML se vezmou reakce a geny organismu. U map organismu (např. `gga`)
   se z KEGG REST API stáhne referenční mapa `ko<číslo>` se všemi reakcemi a k ní EC čísla,
   názvy a struktury (MOL, u glykanů strom cukrů KCF). Vše se ukládá do složky `kegg_cache/`,
   takže podruhé to funguje i offline.
2. **Výběr.** Vybrané látky a reakce se sloučí podle dvojic látek (paralelní enzymy
   = jedna šipka). Skryté látky mezi vybranými se nahradí rámečkem „vynechané kroky“.
3. **Rozložení.** Vrstvy podle vzdálenosti od výchozí látky, pomocné uzly pro dlouhé hrany,
   barycentrická metoda s výměnami sousedů (minimalizace křížení) a svislé polohy
   izotonickou regresí (PAV).
4. **Šipky.** Pravoúhlé vedení, pořadí svislých kolejí v každé mezeře optimalizované
   proti křížení, popisky s kontrolou kolizí.
5. **Výstup.** Struktury z RDKit (text atomů jako křivky), ostatní text převedený na obrysy
   přes matplotlib `TextPath`. V SVG není žádný `<text>`.

## Data KEGG a licence

Program pracuje s daty KEGG, která si každý uživatel stahuje sám (KGML soubory a
[KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html)). KEGG je zdarma pro akademické
použití. Komerční použití a hromadné šíření dat KEGG podléhá
[licenčním podmínkám KEGG](https://www.kegg.jp/kegg/legal.html). Repozitář proto
neobsahuje soubory KGML ani mezipaměť `kegg_cache/`.

Symboly glykanů podle SNFG: Varki A. et al., *Glycobiology* 25, 1323–1324 (2015).
