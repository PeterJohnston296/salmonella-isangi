#!/usr/bin/env python3
from pathlib import Path
import argparse, csv, gzip, hashlib, json, re, sys
import numpy as np
import pandas as pd
from Bio import Phylo

ROOT = Path(__file__).resolve().parents[1]
AGENTS = "AMP10 AZM15 CIP5 FOX30 MEM10 C30 PEF5 FEP30 SXT25 TCG15 IPM10 AK30 CPD10 ATM30 CN10 TZP36 CTX5".split()
EXTRA = {"ERR9939684", "ERR11201713", "ERR9939690", "SRR3049609"}

def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()

def table(name):
    return pd.read_csv(ROOT / "data/data_clean" / name, dtype=str, keep_default_na=False)

def matrix(name):
    d = pd.read_csv(ROOT / "data/matrices" / name, index_col=0)
    if len(d.index) != len(d.columns) or d.index.has_duplicates or d.columns.has_duplicates:
        raise ValueError(f"Invalid matrix identifiers: {name}")
    if set(d.index) != set(d.columns):
        raise ValueError(f"Matrix row/column membership: {name}")
    d = d.loc[d.index, d.index]
    x = d.to_numpy(dtype=float)
    if not (np.isfinite(x).all() and (x >= 0).all() and np.array_equal(x, x.T)
            and (np.diag(x) == 0).all() and np.equal(x, np.round(x)).all()):
        raise ValueError(f"Invalid pairwise distances: {name}")
    return d

def tree_ids(name):
    return [x.name for x in Phylo.read(ROOT / "data/trees" / name, "newick").get_terminals()]

def checked(condition, name, records):
    records.append({"check": name, "result": "pass" if bool(condition) else "fail"})
    if not condition:
        raise ValueError(name)

def validate():
    results = []
    m = table("isolate_accession_line_list_345.csv").set_index("isolate_id", drop=False)
    checked(len(m) == 345 and m.index.is_unique, "345 unique publication isolates", results)
    checked(m.run_accession.str.fullmatch(r"[EDS]RR\d+").all() and m.run_accession.is_unique,
            "345 unique run accessions", results)
    e = table("ENA_deposited_isolates_80.csv")
    run_col = next(c for c in e if c.lower().replace('_',' ') in ("run accession", "run level accession"))
    checked(len(e) == 80 and e[run_col].is_unique and set(e[run_col]) <= set(m.run_accession),
            "80 deposited isolates link to publication cohort", results)
    cc_ids = tree_ids("CC25_tree.nwk")
    st_ids = tree_ids("ST335_tree.nwk")
    pub_ids = tree_ids("CC25_publication_tree.nwk")
    checked(len(cc_ids) == 327 and len(set(cc_ids)) == 327, "327 unique CC25 tree tips", results)
    checked(len(st_ids) == 224 and len(set(st_ids)) == 224 and set(st_ids) < set(cc_ids),
            "224 ST335 tips are an exact CC25 subset", results)
    checked(set(pub_ids) == set(cc_ids), "Approved circular tree has the CC25 tip set", results)
    checked(m.loc[cc_ids].sequence_type.value_counts().to_dict() == {"335":224,"216":102,"5028":1},
            "CC25 sequence-type counts", results)
    checked(len(set(m.index)-set(cc_ids)) == 18, "18 publication isolates outside the CC25 tree", results)
    full = matrix("CC25_complete_331.csv")
    checked(len(full) == 331 and set(full.index)-set(cc_ids) == EXTRA, "Four explicit matrix-only isolates", results)
    checked(not (EXTRA & set(m.index)), "Matrix-only isolates not inserted in the publication cohort", results)
    sizes = {"CC25_tree_aligned_327.csv":327,"ST335_tree_aligned_224.csv":224,
             "Malawi_ST335_74.csv":74,"South_Africa_outbreak_ST335_37.csv":37,
             "Malawi_74_source_matrix.csv":74}
    for n, size in sizes.items():
        d = matrix(n)
        checked(len(d) == size and np.array_equal(d, full.loc[d.index, d.columns]), f"{n}: exact symmetric subset", results)
    mw = matrix("Malawi_ST335_74.csv")
    sa = matrix("South_Africa_outbreak_ST335_37.csv")
    checked(int(mw.to_numpy().max()) == 35 and int(sa.to_numpy().max()) == 108, "Malawi and South Africa SNP ranges", results)
    checked((m.loc[mw.index].country == "Malawi").all() and (m.loc[sa.index].country == "South Africa").all(),
            "Country-specific matrix membership", results)
    for p in sorted((ROOT/"data/matrices").glob("*-*.csv")):
        d = matrix(p.name)
        checked(np.array_equal(d, full.loc[d.index, d.columns]), f"{p.stem}: exact detail submatrix", results)
    panels = table("panel_isolate_index.csv")
    for pop, ids in (("CC25", cc_ids),("ST335", st_ids)):
        q = panels[panels.population == pop]
        checked(q.isolate_id.is_unique and set(q.isolate_id) == set(ids), f"{pop} detail panels cover every tip once", results)
    for pid, q in panels.groupby("panel_id", sort=False):
        lineages = m.loc[q.isolate_id].lineage.unique()
        pos = m.loc[q.isolate_id, "CC25_order" if pid.startswith("CC25") else "ST335_order"].astype(float).sort_values()
        checked(len(lineages) == 1 and np.all(np.diff(pos) == 1), f"{pid}: consecutive within-lineage panel", results)
    calls = table("phenotype_categories_as_recorded_42.csv").set_index("isolate_id", drop=False)
    zones = table("phenotype_zone_diameters_42.csv").set_index("isolate_id", drop=False)
    checked(len(calls) == 42 and calls.index.is_unique and set(calls.index) == set(zones.index), "42 paired AST records", results)
    checked(all(a in calls and a in zones for a in AGENTS), "All 17 recorded antibiotics retained", results)
    checked(set(calls[AGENTS].to_numpy().flat) <= {"S","I","R"}, "Recorded AST calls preserved", results)
    checked(np.isfinite(zones[AGENTS].astype(float)).all().all() and (zones[AGENTS].astype(float)>=0).all().all(),
            "714 numeric inhibition-zone measurements", results)
    i = calls[calls.group == "Isangi"]
    xdr = (i[["AMP10","C30","SXT25","PEF5"]] == "R").all(axis=1) & (i[["CPD10","CTX5"]] == "R").any(axis=1)
    checked(len(i)==19 and len(calls)-len(i)==23 and int(xdr.sum())==12, "19 Isangi, 23 comparators, 12 phenotypic XDR", results)
    checked(set(i.index) <= set(st_ids), "Every phenotyped Isangi links to ST335", results)
    checked((i.MEM10 == "S").all() and (i.IPM10 == "S").all() and (i.AK30 == "S").all(), "Recorded retained-susceptibility counts", results)
    sc = json.loads((ROOT/"data/data_clean/scientific_config.json").read_text())
    checked(sc["antibiotics"] == AGENTS, "Antibiotic dictionary order", results)
    g = table("genomic_determinants_345.csv").set_index("isolate_id")[sc["gene_order"]].astype(int)
    checked(len(g)==345 and set(g.index)==set(m.index) and len(g.columns)==55 and g.isin([0,1]).all().all(),
            "345 by 55 explicit curated genomic calls", results)
    flags = {k:g[[v for v in genes if v in g]].any(axis=1) for k,genes in sc["genomic_category_markers"].items()}
    mdr = flags["Penicillin"] & flags["Chloramphenicol"] & flags["Trimethoprim"] & flags["Sulfonamide"]
    gxdr = mdr & flags["Fluoroquinolone"] & flags["Third-generation cephalosporin"]
    checked(int(gxdr.loc[st_ids].sum()) == 199, "199 ST335 genomic XDR-associated profiles", results)
    checked(np.array_equal(gxdr.loc[m.index].astype(int), m.XDR_genomic_profile.astype(int)), "Genomic profiles match all 345 stored flags", results)
    checked(int(flags["Carbapenem"].loc[st_ids].sum()) == 6, "Six ST335 carbapenemase-marker profiles", results)
    checked(int(gxdr.loc[i.index].sum()) == 19, "Matched 19 genomic profiles retained separately from 12 phenotypes", results)
    annual=table("BSI_annual_resistance_summary_2015_2023.csv")
    for yr,q in annual.groupby("year"):
        total=q.annual_total.astype(int).unique()
        checked(len(total)==1 and q.n.astype(int).sum()==total[0] and np.allclose(q.percent.astype(float),q.n.astype(int)/total[0]*100),
                f"Annual surveillance counts and percentages: {yr}", results)
    raw=pd.read_csv(ROOT/"data/data_raw/ebg_25_snpmatrix.tsv",sep="\t",index_col=0)
    clean=lambda s:re.sub(r"_S\d+_L\d+$","",s.replace("#","_"))
    raw.index=[clean(str(s)) for s in raw.index];raw.columns=[clean(str(s)) for s in raw.columns]
    checked(set(raw.index)==set(full.index) and np.array_equal(raw.loc[full.index,full.columns],full), "Original 331-source matrix matches canonical numeric archive", results)
    return pd.DataFrame(results)

def burden():
    raw=ROOT/"data/data_raw"
    a=pd.read_csv(raw/"all_isangi_assemblies_amrfinder2_20062024.csv",dtype=str,keep_default_na=False)
    lookup=pd.read_csv(raw/"isangi_st_lookup.csv",dtype=str,keep_default_na=False)
    a["isolate_id"]=a.Isolate_ID.str.replace(r"_assembly_contigs$","",regex=True)
    known=set(a.isolate_id)
    if not known <= set(lookup.isolate_id):raise ValueError("AMRFinder isolate membership differs from lookup")
    a=a[(a["Element type"]=="AMR")&(a["Element subtype"]=="AMR")&(a["Gene symbol"]!="")]
    counts=a.drop_duplicates(["isolate_id","Gene symbol"]).groupby("isolate_id").size()
    study=lookup.copy()
    study["present_in_feature_file"]=study.isolate_id.isin(known)
    study["amr_gene_count"]=study.isolate_id.map(counts).fillna(0).astype(int)
    study["group"]=np.where(study.ST=="335","Isangi ST335","Other Isangi")
    frequencies={};n_all=0;n_excluded=0;n_unresolved=0;n_blank=0;n_empty_rows=0
    for c in pd.read_csv(raw/"2026.01.13.ncbi_pathogen_salmonella_amr_serovar.tsv.tsv.gz",sep="\t",dtype=str,keep_default_na=False,chunksize=100000):
        empty=c.eq("").all(axis=1)
        n_empty_rows+=int(empty.sum())
        c=c.loc[~empty]
        n_all+=len(c)
        is_isangi=c["Computed types"].str.contains("serotype=Isangi",case=False,regex=False)
        n_excluded+=int(is_isangi.sum());c=c[~is_isangi]
        n_unresolved+=int((~c["Computed types"].str.contains(r"serotype=[^,\s]+",regex=True)).sum())
        s=c["#AMR genotypes"]
        n_blank+=int((s=="").sum())
        v=s.str.count("=COMPLETE")-s.str.contains(r"(?:^|,)mdsA=COMPLETE(?:,|$)",regex=True).astype(int)-s.str.contains(r"(?:^|,)mdsB=COMPLETE(?:,|$)",regex=True).astype(int)
        for val,n in v.clip(lower=0).value_counts().items():frequencies[int(val)]=frequencies.get(int(val),0)+int(n)
    summary=[]
    values=np.repeat(sorted(frequencies),[frequencies[k] for k in sorted(frequencies)])
    for label, vals in [("NCBI non-Isangi",values)]+[(g,t.amr_gene_count.to_numpy()) for g,t in study.groupby("group",sort=False)]:
        q=np.quantile(vals,[.25,.5,.75])
        summary.append(dict(group=label,n=len(vals),minimum=int(min(vals)),q1=float(q[0]),median=float(q[1]),q3=float(q[2]),maximum=int(max(vals))))
    freq=[dict(group="NCBI non-Isangi",amr_gene_count=k,n=frequencies[k]) for k in sorted(frequencies)]
    for g,t in study.groupby("group",sort=False):
        freq += [dict(group=g,amr_gene_count=int(k),n=int(n)) for k,n in t.amr_gene_count.value_counts().sort_index().items()]
    note=dict(NCBI_total_records=n_all,wholly_empty_TSV_records_ignored=n_empty_rows,NCBI_Isangi_excluded=n_excluded,NCBI_comparator_records=len(values),
        NCBI_unresolved_serovar_in_comparator=n_unresolved,NCBI_blank_genotype_cells=n_blank,
        study_records=len(study),study_records_absent_from_hit_file=int((~study.present_in_feature_file).sum()),
        other_Isangi_unassigned_ST=int(((study.group=="Other Isangi")&~study.ST.str.fullmatch(r"\d+")).sum()))
    return study,pd.DataFrame(summary),pd.DataFrame(freq),note

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument("--burden",action="store_true")
    args=parser.parse_args()
    out=ROOT/"results/validation";out.mkdir(parents=True,exist_ok=True)
    results=validate();results.to_csv(out/"data_checks.csv",index=False)
    if args.burden:
        study,summary,frequency,note=burden()
        study.to_csv(out/"Figure_3B_study_counts.csv",index=False)
        summary.to_csv(out/"Figure_3B_summary.csv",index=False)
        frequency.to_csv(out/"Figure_3B_frequencies.csv",index=False)
        (out/"Figure_3B_input_audit.json").write_text(json.dumps(note,indent=2)+"\n")
        print(summary.to_string(index=False))
    print(f"{len(results)} data checks passed.")
if __name__ == "__main__":main()
