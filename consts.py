import pandas as pd
import os
from importlib import reload

PATH = os.getcwd()
MT_GENES = ("ND1",  "ND2",  "ND3",  "ND4", "ND4L", "ND5", "ND6", "CO1",  "CO2",  "CO3",  "CYB", 'CYTB',  "ATP6", "ATP8", 'MT-RNR1', 'MT-RNR2')
MT_GENES = ['MT-' + gene if 'MT-' not in gene else gene for gene in MT_GENES]

# Make MT_GENES without MT- prefix
MT_GENES_NO_PREFIX = [gene.replace('MT-', '') for gene in MT_GENES]
MT_GENES = MT_GENES + MT_GENES_NO_PREFIX

def genelist(file):
    """
    Generate a dataframe based on a list of genes from a file    
    """
    prefix = os.path.split(file)[-1].replace('.txt', '').replace('.tsv', '')
    if file.endswith('.txt'):
        with open(file) as f:
            genes = f.read().splitlines()
        df = pd.DataFrame({'gene':genes, 'Description' : [prefix for _ in genes]})
    elif file.endswith('.tsv'):
        df = pd.read_csv(file, delimiter='\t')
        df['Description'] = prefix
    return df

def genelist_on_folder(folder, mode = '.txt'):
    """
    Generate a dataframe based on a list of genes from a folder
    """
    files = os.listdir(rf'{folder}')
    if mode == '.txt':
        files = [f for f in files if f.endswith('.txt')]
    elif mode == '.tsv':
        files = [f for f in files if f.endswith('.tsv')]
    df = pd.DataFrame()
    for i, file in enumerate(files):
        if i == 0:
            df = genelist(os.path.join(folder, file))
        else:
            df = pd.concat([df, genelist(os.path.join(folder, file))])
    return df

PEPPRO_RANGES = {
    'Degradation_ratio' : [0,1],
    'Pause_index' : [10, 999],
    'mRNA_contamination' : [0, 3.5],
    'Pct_uninformative_adapter_reads' : [0, 0.25],
    'TSS_coding_score' : [5, 999],
    'TSS_non-coding_score' : [2, 999]
    }

COMPLEX_ORDER = ['I', 'II', 'III', 'IV', 'V', 'Ribosome']

DROS_FB_MITO = pd.read_excel(os.path.join(PATH, 'data', 'genes', 'mito.genes.xlsx'))

GENELIST_DF = genelist_on_folder(os.path.join(PATH, 'data', 'genes'), mode ='.tsv')[['Name', 'Complex', 'Description', 'Protein_type']].dropna()
GENELIST_DF['Genome'] = GENELIST_DF['Name'].apply(lambda x: 'MT' if x in MT_GENES else 'NUC')
UNIFY_DICT = {'Control': 'Normoxia', '24h' : 'Hypoxia', ' normoxia' : 'Normoxia', ' hypoxia' : 'Hypoxia', '1%hypoxia':'Hypoxia', 'Hpx' : 'Hypoxia', 'DMSO' : 'Normoxia'}

AP1 = 'CCN1,MAFG,BCL2L11,MAF,NFATC1,PLAU,AGT,TIMP1,FOS,MYC,TGFB1,NPPA,PENK,IFNG,MT2A,ESR1,MMP1,CSF2,NR3C1,HLAA,TP53,IL4,IL5,IL6,EDN1,JUN,CDK1,TH,SP1,COL1A2,CXCL8,MYB,CCL2,MMP9,ETS1,FABP4,ATF2,FOSL1,FOSL2,CREB1,JUNB,GJA1,JUND,EGR1,ATF3,IL10,GATA2,CCND1,DUSP1,NTS,ELF1,CTNNB1,CDKN2A,CDKN1B,FOSB,PTEN,IL2,ACTA1,EP300,NFATC3,DMP1,NFATC2,CBFB,TRIP6,HIF1A,CRTC1,COPS5,BAG1,TCF7L2,DMTF1'.split(',')

N_NUC_OXPHOS_GENES = 88
N_MT_OXPHOS_GENES = 13
N_NUC_RIBOSOME_GENES = 78
N_MT_RIBOSOME_GENES = 2

MTDNA_GENE_ORDER = ['RNR1', 'RNR2', 'ND1', 'ND2', 'CO1', 'CO2', 'ATP8', 'ATP6', 'CO3', 'ND3', 'ND4L', 'ND4', 'ND5', 'CYB', 'ND6']
# Add MT to genes without MT- prefix
MTDNA_GENE_ORDER = ['MT-' + gene if 'MT-' not in gene else gene for gene in MTDNA_GENE_ORDER]
