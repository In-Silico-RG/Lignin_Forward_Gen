#!/usr/bin/env python
"""
2026-09-14: append the candidate rows read from docs/13_candidate_13C_inventories.md to
data/reference_inventories.csv with use_fit = 0 and vetted = pending. RUN ONCE (commit
d92f1cc); kept for provenance. The numbers are exactly those transcribed in docs/13,
sections 2 and 4; mapping conventions are in each row's `notes`.
"""
import csv
import os

P = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, "reference_inventories.csv")
hdr = next(csv.reader(open(P)))
rows = []


def add(sid, source, doi, table, species, tissue, typ, iso, method, sg, sgm, h, denom,
        bo4=None, b5=None, bb=None, b1=None, c55=None, c4o5=None, dbdo=None, dc=None, ca=None,
        use_val=1, notes=""):
    r = dict.fromkeys(hdr, "")
    r.update(dict(sample_id=sid, source=source, doi=doi, table=table, species=species, tissue=tissue,
                  type=typ, isolation=iso, method=method, sg_ratio=sg, sg_method=sgm, h_frac=h,
                  denominator=denom, use_fit=0, use_val=use_val, vetted="pending"))
    for k, v in dict(bo4=bo4, b5=b5, bb=bb, b1=b1, c55=c55, c4o5=c4o5, dbdo=dbdo, dc=dc, ca_end=ca).items():
        if v is not None:
            r[k] = v
    r["notes"] = notes
    rows.append(r)


# --- Sette 2011 QQ-HSQC per 100 C9 (docs/13 sec 4.3, Table 2)
S = "Sette Wechselberger Crestini 2011 Chem Eur J 17:9529"; D = "10.1002/chem.201003045"
add("sette2011_norwayspruce_mwl", S, D, "Table 2", "Picea abies", "wood", "softwood", "MWL", "HSQC-quant", 0.0, "exact(no S)", 0, "per100Ar", 44.7, 10.6, 3.2, 1.3, dbdo=3.1, notes="QQ-HSQC acetylated MWL; stated max rel. errors bo4 2%, b5 19%, bb 1%, dbdo 27%, b1 2%; 5-5/4-O-5 not HSQC-detectable")
add("sette2011_blackspruce_mwl", S, D, "Table 2", "Picea mariana", "wood", "softwood", "MWL", "HSQC-quant", 0.0, "exact(no S)", 0, "per100Ar", 49.4, 13.1, 5.8, 1.1, dbdo=2.0, notes="QQ-HSQC")
add("sette2011_norwayspruce_emal", S, D, "Table 2", "Picea abies", "wood", "softwood", "EMAL", "HSQC-quant", 0.0, "exact(no S)", 0, "per100Ar", 47.0, 10.6, 1.6, 2.6, dbdo=3.8, notes="QQ-HSQC; COST E41 round-robin EMAL")
add("sette2011_beech_mwl", S, D, "Table 2", "Fagus sylvatica (beech)", "wood", "hardwood", "MWL", "HSQC-quant", 1.67, "HSQC", 0, "per100Ar", 60.3, 1.0, 8.0, 3.5, dbdo=0.5, notes="QQ-HSQC; S/G 1.67 by QQ-HSQC (text); H not reported (blank -> 0 assumed)")
# --- Sette 2013 spruce % C9 (sec 4.1)
add("sette2013_norwayspruce_qqhsqc", "Sette Lange Crestini 2013 CSBJ 6:e201303016", "10.5936/csbj.201303016", "Table 1", "Picea abies", "wood", "softwood", "MWL", "HSQC-quant", 0.0, "exact(no S)", 0, "per100Ar", 40.1, 11.8, 6.7, 1.6, dbdo=2.5, ca=2.2, notes="QQ-HSQC column (HSQC0 column: 38.2/9.1/7.1/1.4/2.3); G2 = 100 internal standard; acetylated")
# --- Rencoret 2019 olive per 100 Ar (sec 4.5a)
R19 = "Rencoret Gutierrez Castro del Rio 2019 Holzforschung 73:25"; D19 = "10.1515/hf-2018-0077"
add("rencoret2019_olive_mwl", R19, D19, "Table 2", "Olea europaea pruning wood", "wood", "hardwood", "MWL", "HSQC-quant", 0.9, "HSQC", 0.01, "per100Ar", 52, 9, 6, 1, dbdo=1, ca=2, notes="HSQC per 100 aromatic units; H 1 G 51 S 48; spirodienone 1 -> b1")
add("rencoret2019_olive_wcw", R19, D19, "Table 2", "Olea europaea pruning wood", "wood", "hardwood", "WCW", "HSQC-quant", 1.0, "HSQC", 0.02, "per100Ar", 53, 8, 6, ca=4, notes="whole cell wall gel; dbdo and spirodienone not detected (dash)")
# --- Rencoret 2009 Paulownia % of side chains (sec 4.4)
add("rencoret2009_paulownia_mwl", "Rencoret et al. 2009 Ind Crops Prod 30:137", "10.1016/j.indcrop.2009.03.004", "Table 3", "Paulownia fortunei", "wood", "hardwood", "MWL", "HSQC-rel", 0.66, "HSQC", 0.0, "pct_of_sidechains", 62, 11, 12, 4, dbdo=2, ca=4, notes="b1 = spirodienone 3 + beta-1 1; H 0 G 60 S 40; Py-GC/MS S/G 0.66-0.67")
# --- Rencoret 2011 E. globulus per 100 Ar (sec 2.2), MWL columns
R11 = "Rencoret et al. 2011 Plant Physiol 155:667"; D11 = "10.1104/pp.110.167254"
add("rencoret2011_eglobulus_1mo_mwl", R11, D11, "Table V", "Eucalyptus globulus 1-month", "wood (young)", "hardwood", "MWL", "HSQC-quant", 1.33, "HSQC", 0.02, "per100Ar", 51, 4, 13, 3, ca=2, notes="b1 = spirodienone 2 + beta-1 1; H 2 G 42 S 56; no dbdo row; alpha-oxidised bo4 1 not added")
add("rencoret2011_eglobulus_18mo_mwl", R11, D11, "Table V", "Eucalyptus globulus 18-month", "wood", "hardwood", "MWL", "HSQC-quant", 2.70, "HSQC", 0.0, "per100Ar", 51, 2, 14, 4, ca=1, notes="b1 = spirodienone 3 + beta-1 1")
add("rencoret2011_eglobulus_9y_mwl", R11, D11, "Table V", "Eucalyptus globulus 9-year", "wood", "hardwood", "MWL", "HSQC-quant", 3.17, "HSQC", 0.0, "per100Ar", 50, 1, 14, 4, ca=1, notes="b1 = spirodienone 3 + beta-1 1")
# --- Rosado 2023 citrus (sec 2.3)
add("rosado2023_citrus_mwl", "Rosado et al. 2023 Polymers 15:1840", "10.3390/polym15081840", "Table 4", "Citrus sinensis pruning wood", "wood", "hardwood", "MWL", "HSQC-rel", 0.60, "HSQC", 0.01, "pct_of_linkages", 70, 15, 9, 3, dbdo=3, ca=6, notes="H 1 G 62 S 37; Py-GC/MS S/G 0.73")
# --- Lourenco 2016 Quercus suber (sec 2.4)
L16 = "Lourenco et al. 2016 Front Plant Sci 7:1612"; DL = "10.3389/fpls.2016.01612"
add("lourenco2016_qsuber_xylem_ml", L16, DL, "Table 4", "Quercus suber xylem", "wood", "hardwood", "ML", "HSQC-rel", 1.2, "HSQC", 0.01, "pct_of_linkages", 77, 9, 8, 3, dbdo=1, ca=1, notes="b1 = spirodienone 2 + open beta-1 1; alpha-ox bo4 2 not added; Py-GC/MS S/G 1.66")
add("lourenco2016_qsuber_phloem_ml", L16, DL, "Table 4", "Quercus suber phloem", "phloem", "other", "ML", "HSQC-rel", 0.7, "HSQC", 0.01, "pct_of_linkages", 71, 13, 7, 5, dbdo=2, ca=2, notes="phloem, not xylem; b1 = spirodienone 3 + open beta-1 2")
add("lourenco2016_qsuber_cork_ml", L16, DL, "Table 4", "Quercus suber cork (phellem)", "cork", "other", "ML", "HSQC-rel", 0.1, "HSQC", 0.02, "pct_of_linkages", 68, 20, 4, 3, dbdo=5, ca=8, notes="cork; 48% gamma-acylation; ferulates; validation-only")
# --- Rosado 2021 rice (grass)
R21 = "Rosado et al. 2021 Front Plant Sci 12:640475"; D21 = "10.3389/fpls.2021.640475"
add("rosado2021_rice_husk_mwl", R21, D21, "Table 4", "Oryza sativa husk", "husk", "grass", "MWL", "HSQC-rel", 0.15, "HSQC", 0.07, "pct_of_linkages", 65, 23, 4, 3, dbdo=5, ca=6, notes="grass: tricin 7, pCA 12, gamma-acyl 10%; bb = resinol 4 + THF 0")
add("rosado2021_rice_straw_mwl", R21, D21, "Table 4", "Oryza sativa straw", "straw", "grass", "MWL", "HSQC-rel", 0.34, "HSQC", 0.05, "pct_of_linkages", 78, 12, 5, 1, dbdo=4, ca=6, notes="grass: tricin 27, pCA 16, gamma-acyl 12%; bb = resinol 4 + THF 1")
# --- Wang & Deuss 2021 MWL rows (mol % of C9; bo4,b5,bb only)
W21 = "Wang & Deuss 2021 ChemSusChem 14:5186"; DW = "10.1002/cssc.202101527"
add("wangdeuss2021_birch_mwl", W21, DW, "Table 1", "Betula (birch)", "wood", "hardwood", "MWL", "HSQC-rel", 2.0, "HSQC", 0.03, "per100Ar", 53, 4, 11, notes="semi-quantitative HSQC, mol% of total C9 units (S+S'+G+G'+H); treated as HSQC-rel sigma")
add("wangdeuss2021_pine_mwl", W21, DW, "Table 1", "Pinus (pine)", "wood", "softwood", "MWL", "HSQC-rel", 0.0, "exact(no S)", 0.02, "per100Ar", 34, 17, 6, notes="as above")
add("wangdeuss2021_walnut_mwl", W21, DW, "Table 1", "Juglans (walnut) shell", "shell", "other", "MWL", "HSQC-rel", 1.1, "HSQC", 0.15, "per100Ar", 48, 13, 10, notes="nut shell, H 15%; validation-only")
add("wangdeuss2021_reed_mwl", W21, DW, "Table 1", "Phragmites (reed)", "stem", "grass", "MWL", "HSQC-rel", 0.8, "HSQC", 0.03, "per100Ar", 45, 9, 4, notes="grass: FA 3, pCA 20; validation-only")
# --- Stewart 2009 poplar WT and F5H (AcML; % of quantified units)
S09 = "Stewart et al. 2009 Plant Physiol 150:621"; DS = "10.1104/pp.109.137059"
add("stewart2009_poplar_wt_acml", S09, DS, "Table II", "Populus tremula x alba WT", "wood", "hardwood", "ML(acetylated)", "HSQC-rel", 2.24, "HSQC", 0.0, "pct_of_linkages", 87.4, 2.0, 7.9, 2.1, ca=0.7, notes="semi-quantitative HSQC; bo4 = A + A3; b1 = spirodienone; dbdo ~0.3% of bo4; PB 11.8% not a linkage")
add("stewart2009_poplar_f5h_acml", S09, DS, "Table II", "Populus tremula x alba C4H::F5H", "wood", "hardwood (transgenic)", "ML(acetylated)", "HSQC-rel", 38, "HSQC", 0.0, "pct_of_linkages", 86.1, 0.1, 9.7, 3.3, ca=0.9, notes="transgenic high-S (S/G 38 NMR, 14-20 by thio/DFRC); no dbdo detectable")
# --- Ralph 2012 C3H poplar (% of A+B+C)
R12 = "Ralph et al. 2012 BioEnergy Res 5:1009"; DR = "10.1007/s12155-012-9218-y"
add("ralph2012_poplar_wt_cel", R12, DR, "Table 1", "Populus WT", "wood", "hardwood", "CEL(acetylated)", "HSQC-rel", 1.70, "HSQC", 0.003, "pct_of_linkages", 88.8, 4.9, 6.3, ca=1.1, notes="% of A+B+C; H:G:S 0.3:37:63")
add("ralph2012_poplar_c3h14_cel", R12, DR, "Table 1", "Populus C3H-14 (C3H-downregulated)", "wood", "hardwood (transgenic)", "CEL(acetylated)", "HSQC-rel", 3.6, "HSQC", 0.31, "pct_of_linkages", 88.9, 4.3, 6.8, ca=1.1, notes="transgenic: H 31%, G 15, S 54 -> S/G 3.6 with pct_h 0.31")
# --- Peng 2021 (per 100 Ar, quantitative 2D-HSQC); bb later halved by the defaults script
P21 = "Peng et al. 2021 Front Bioeng Biotechnol 9:790539"; DP = "10.3389/fbioe.2021.790539"
add("peng2021_poplar_ck_del", P21, DP, "Table 4", "Populus alba x glandulosa control", "wood", "hardwood", "DEL", "HSQC-quant", 2.82, "HSQC", 0.006, "per100Ar", 57.79, 1.18, 13.02, notes="double enzymatic lignin; quantitative 2D-HSQC per 100 Ar; bb likely units-in-resinol convention -- CHECK (13.02 per 100 Ar is high for bonds)")
add("peng2021_poplar_c3h_del", P21, DP, "Table 4", "Populus alba x glandulosa C3H-down", "wood", "hardwood (transgenic)", "DEL", "HSQC-quant", 2.48, "HSQC", 0.011, "per100Ar", 56.78, 1.14, 10.76, notes="as above")
# --- Dou 2018 willow wood (% of A+B+C)
add("dou2018_willow_wood_wcw", "Dou et al. 2018 J Agric Food Chem 66:7294", "10.1021/acs.jafc.8b02014", "Table 4", "Salix (willow) wood", "wood", "hardwood", "WCW", "HSQC-rel", 2.0, "HSQC", 0.0, "pct_of_linkages", 89, 3, 7, notes="% of A+B+C; EL S/G 1.9; thio 1.4-2.1")
# --- del Rio 2012 elephant grass MWL (grass)
E12 = "del Rio et al. 2012 J Agric Food Chem 60:3619"; DE = "10.1021/jf300099g"
add("delrio2012_elephantgrass_cortex_mwl", E12, DE, "Table 6", "Pennisetum purpureum cortex", "stem", "grass", "MWL", "HSQC-rel", 1.4, "HSQC", 0.03, "pct_of_linkages", 82, 8, 5, 2, ca=6, notes="grass: bb = resinol 2 + THF 3; pCA 29%, gamma-acyl 39%; alpha-ox bo4 2 not added")
add("delrio2012_elephantgrass_pith_mwl", E12, DE, "Table 6", "Pennisetum purpureum pith", "stem", "grass", "MWL", "HSQC-rel", 1.5, "HSQC", 0.03, "pct_of_linkages", 82, 7, 6, 2, ca=7, notes="grass: bb = resinol 1 + THF 5; pCA 40%, gamma-acyl 55%")
# --- Happs 2021: poplar genotypes (% of bo4+bb+b5), whole biomass
H21 = "Happs et al. 2021 Biotechnol Biofuels 14:1897"; DH = "10.1186/s13068-021-01897-y"
happs = [("004", 1.6, 70, 24, 6.1), ("021", 0.89, 66, 25, 9.4), ("036", 1.9, 74, 21, 4.8), ("075", 2.1, 72, 23, 4.6),
         ("095", 1.5, 75, 20, 5.2), ("096", 2.1, 76, 20, 3.4), ("140", 1.7, 74, 22, 4.3), ("169", 1.2, 73, 21, 6.5),
         ("173", 1.1, 69, 23, 8.2), ("182", 1.8, 73, 21, 5.9), ("217", 1.7, 72, 23, 5.3), ("219", 2.3, 72, 25, 2.8),
         ("255", 1.8, 71, 22, 6.2), ("282", 1.1, 62, 25, 14), ("322", 2.2, 75, 22, 3.2), ("334", 2.3, 69, 27, 4.5),
         ("388", 0.94, 67, 22, 11), ("841", 2.4, 73, 23, 3.8), ("853", 2.4, 73, 23, 3.3), ("863", 1.4, 73, 23, 4.1),
         ("883", 2.3, 76, 21, 2.5)]
for g, sg, bo4, bb, b5 in happs:
    add("happs2021_besc%s" % g, H21, DH, "Table 3", "Populus trichocarpa BESC-%s" % g, "wood", "hardwood", "whole biomass", "HSQC-rel", sg, "HSQC", "", "pct_of_linkages", bo4, b5, bb, notes="% of (bo4+bb+b5), whole ball-milled biomass; S/G by HSQC differs from thio/py-MBMS for the same sample")
add("happs2021_nist_pine", H21, DH, "Table 3", "Pinus (NIST 8493)", "wood", "softwood", "whole biomass", "HSQC-rel", 0.0, "exact(no S)", "", "pct_of_linkages", 57, 33, 10, notes="% of (bo4+bb+b5)")

if __name__ == "__main__":
    raise SystemExit("Already applied on 2026-09-14 (commit d92f1cc). Kept for provenance; do not re-run.")
    with open(P, "a", newline="") as fh:
        csv.DictWriter(fh, fieldnames=hdr).writerows(rows)
