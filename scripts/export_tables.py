#!/usr/bin/env python3
from pathlib import Path
import json
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
DATA=ROOT/'data/data_clean'
OUT=ROOT/'results/tables'
OUT.mkdir(parents=True,exist_ok=True)
def get(n):return pd.read_csv(DATA/n,dtype=str,keep_default_na=False)
def save(n,x):x.to_csv(OUT/f'Table_S{n:02d}.csv',index=False)
m=get('isolate_accession_line_list_345.csv').set_index('isolate_id',drop=False)
p=get('phenotype_categories_as_recorded_42.csv').set_index('isolate_id',drop=False)
z=get('phenotype_zone_diameters_42.csv').set_index('isolate_id',drop=False)
agents=json.loads((DATA/'scientific_config.json').read_text())['antibiotics']
ids=p.index[p.group=='Isangi']
def ast(x, cols):
 out=x.loc[ids,['isolate_id']+cols].copy();out.insert(1,'run_accession',m.loc[ids,'run_accession']);return out
save(1,get('neonatal_unit_isolates_37.csv'))
save(2,get('disinfectant_MIC_MBC_original.csv'))
b=get('Kenneth_recorded_biocide_results_42.csv')
save(3,b[b['group']=='Isangi'])
save(4,b[b['group']=='Comparator'])
save(5,m.reset_index(drop=True))
save(6,get('ENA_deposited_isolates_80.csv'))
save(7,ast(z,agents[:9]));save(8,ast(z,agents[9:]))
save(9,ast(p,agents[:9]));save(10,ast(p,agents[9:]))
x=ast(p,['AMP10','C30','SXT25','PEF5','CPD10','CTX5'])
x['MDR_phenotype']=(x[['AMP10','C30','SXT25']]=='R').all(axis=1)
x['XDR_phenotype']=x.MDR_phenotype&(x.PEF5=='R')&(x[['CPD10','CTX5']]=='R').any(axis=1)
x['genomic_XDR_profile']=m.loc[ids,'XDR_genomic_profile'].to_numpy()
save(11,x)
save(12,get('antibiotic_dictionary_17.csv'))
save(13,get('genomic_category_marker_dictionary.csv'))
save(14,m[['isolate_id','run_accession','XDR_genomic_profile','genomic_determinants']])
save(15,pd.DataFrame({'population':['Study collection','CC25 tree','CC25 complete matrix','ST335 tree','Malawi ST335','South African outbreak-associated ST335'], 'n':[345,327,331,224,74,37]}))
save(16,get('isolates_outside_CC25_tree_18.csv'))
q=get('panel_summary.csv');save(17,q[q.population=='CC25'])
save(18,m.loc[ids,['isolate_id','run_accession','sample_accession','CC25_panel','ST335_panel','XDR_phenotype']])
save(19,m.loc[['CAAXQV','CAAYM5','CIV143E9'],['isolate_id','run_accession','source_group','sequence_type','CC25_panel','ST335_panel']])
print('19 numbered source tables written to results/tables/.')
