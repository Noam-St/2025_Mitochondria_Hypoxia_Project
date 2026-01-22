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
PH_COLORS = {
    'Arthropoda':'tab:orange',
    'Chordata': 'tab:red',
    'Cnidaria': 'tab:blue',
    'Echinodermata': 'tab:green',
    'Mollusca': 'tab:purple',
    'Nematoda': 'tab:brown'
    }

HS = 'Homo sapiens'
SP = 'Strongylocentrotus purpuratus'
DM = 'Drosophila melanogaster'
MM = 'Macaca mulatta'
PT = 'Pan trogolodytes'
DP = 'Daphnia pulex'
NS = 'Nematostella sp.'
AG = 'Anopheles gambiae'
BL = 'Blattella germanica'
PM = 'Petromyzon marinus'
CE = 'Caenorhabditis elegans'
OB = 'Octopus bimaculoides'
BT = 'Bemisia tabaci'
TC = 'Tribolium castaneum'
SO = 'Sepia officinalis'
BM = 'Bombyx mori'
PC = 'Pomacea canaliculata'
SC = 'Saccharomyces cerevisiae'
SP = 'Schizosaccharomyces pombe'
ZN = 'Zootermopsis nevadensis'
PR = 'Pieris rapae'
PA = 'Periplaneta americana'
DK = 'Dendrolimus kikuchii'
CR = 'Carcinoscorpius rotundicauda'
VD = 'Varroa destructor'
PC2 = 'Paralithodes camtschaticus'

DM_TRANSCRIPTS = [
  [19405, 6405, True],
  [9530, 12134, True],
  [3716, 93, False],
  [10891, 6349, False],
  [15077, 11894, False]
]
MT_GENES = ("ND1",  "ND2",  "ND3",  "ND4", "ND4L", "ND5", "ND6", "CO1",  "CO2",  "CO3",  "CYB", 'CYTB',  "ATP6", "ATP8", 'MT-RNR1', 'MT-RNR2')
MT_GENES = ['MT-' + gene if 'MT-' not in gene else gene for gene in MT_GENES]

TRNA_DICT = {
    'tRNA-Ala': 'trnA', 'tRNA-Cys':'trnC','tRNA-Asp':'trnD','tRNA-Glu':'trnE',
    'tRNA-Phe':'trnF','tRNA-Gly':'trnG','tRNA-His':'trnH','tRNA-Ile':'trnI',
    'tRNA-Lys':'trnK','tRNA-Leu':'trnL','tRNA-Met':'trnM','tRNA-Asn':'trnN',
    'tRNA-Pro':'trnP','tRNA-Gln':'trnQ','tRNA-Arg':'trnR','tRNA-Ser':'trnS',
    'tRNA-Thr':'trnT','tRNA-Val':'trnV','tRNA-Trp':'trnW','tRNA-Tyr':'trnY'}

REPLACEMENT_DICT = {'rrnS':['small subunit ribsomal RNA','small ribosomal RNA subunit RNA','12S small subunit ribosomal RNA','small ribosomal subunit RNA','12S large subunit ribosomal RNA''small ribosomal RNA subunit RNA','12S small subunit ribosomal RNA','small ribosomal RNA','small subunit rRNA','rRNA-12S','Small subunit ribosomal RNA','12S ribosomal RNA','12s ribosomal RNA','s-rRNA','rrnS','small subunit ribosomal RNA','12S rRNA','srRNA',
                                       '12S-rRNA','12S','ssu rRNA','rns','12srRNA','12s ribosomal RNA',
                                       'rnS','SSU','small ribosomal RNA subunit RNA', 'rrs', 'RNS', 'rns', 'RNS', 'rnsA', 'rrnaS', '12srRNA', '12SrRNA', 'rrn12S', 'RNR1'],
                  'rrnL':['large subunit ribsomal RNA','large ribosomal RNA subunit RNA','16S large subunit ribosomal RNA','large ribosomal subunit RNA','16S large subunit ribosomal RNA''large ribosomal RNA subunit RNA','16S small subunit ribosomal RNA','large ribosomal RNA','large subunit rRNA','rRNA-16S','Large subunit ribosomal RNA','16S ribosomal RNA','16s ribosomal RNA','l-rRNA','rrnL','large subunit ribosomal RNA','16S rRNA',
                                       'lrRNA','16S-rRNA','16S','lsu rRNA','rnl','l6S ribosomal RNA',
                                       '16SrRNA','rrs','LSU','16S ribosomal RNA','rnL','large ribosonal RNA subunit RNA', 'rrl', 'RNL', 'rnl', 'rrnl', 'RNL', 'rnlA', 'rrnl', 'rnl_b','16srRNA', '16SrRNA', 'rrn16S', 'RNR2'],
                  'cox1':['COX1','coI','cytochrome oxidase subunit-1','cytochrome c subunit I','cytochrome oxidase I','cytochrome C oxidase subunit 1','cytochrome oxidase c subunit 1','Cytochrome c oxidase subunit 1','cytochrome oxidase subunit I','cytochrome oxidase sununit 1','Cytochrome oxidase subunit I','cytochrome c oxidase subunit 1','COI','cytochrome c oxidase subunit I','CO1',
                          'cox1','COXI','coi','cytochrome oxidase subunit 1','cytochrome oxidase subunit I',
                          'Cox1','CO-I','coxI','COI protein'],
                  'cox2':['COX2','cytochrome oxidase subunit-2','cytochrome c subunit II','cytochrome oxidase II','cytochrome C oxidase subunit 2','cytochrome oxidase c subunit 2','Cytochrome c oxidase subunit 2','cytochrome oxidase subunit II','cytochrome oxidase sununit 2','Cytochrome oxidase subunit II','cytochrome c oxidase subunit 2','COII','cytochrome c oxidase subunit II','CO2',
                          'cox2','COXII','coii','cytochrome oxidase subunit 2','cytochrome oxidase subunit II',
                          'Cox2','CO-II','coxII','COII protein'],
                  'cox3':['COX3','CoxIII','cytochrome oxidase subunit-3','cytochrome c subunit III','cytochrome oxidase III','cytochrome C oxidase subunit 3','cytochrome oxidase c subunit 3','Cytochrome c oxidase subunit 3','cytochrome oxidase subunit III','cytochrome oxidase sununit 3','Cytochrome oxidase subunit III','cytochrome c oxidase subunit 3','COIII','cytochrome c oxidase subunit III','CO3',
                          'cox3','COXIII','coiii','cytochrome oxidase subunit 3','cytochrome c oxidase subunit 3',
                          'Cox3','CO-III','coxIII','COIII protein'],
                  'cob':['cytochrome b oxidase','cytochrome B','CYTB','apocytochrome b','cytochrome c oxidase subunit b','CytB','Cytochroome b','cob','cytb','cyt b','Cyt b',
                          'Cytb','Cyt B','cytochrome b','Cb','cyt-B','cytB','CYB','Cyt-b','Cytb protein'],
                  'nad1':['NADH1 dehydrogenase subunit 1','NADH dehydogenase subunit 1','ND1','NADH dehydrogenase subunit-1','cytochrome c oxidase subunits I','NADH dehydrogenase subunits 1','NADH dehydrogenase subunit1','NADH dehydrogenase subunit 1','NADH-ubiquinone oxidoreductase subunit 1','NADH dehyrogenase subunit 1','NADH dehydrogenase subunit I','NADH  dehydrogenase 1','NADH dehydrogenase 1','nad1-1','nad1-0','NADH dehydrogenase subunit 1','nad1','NADH1','nd1','nadh1','NAD1',
                         'NADH dehydrognase subunit I','nadh1','ndh1','Nad1','ND-1','ND1 protein','NADH hydrognase subunit 1','NDI'],
                  'nad2':['nad2_a','NADH2 dehydrogenase subunit 2','NADH dehydogenase subunit 2','ND2','NADH dehydrogenase subunit-2','cytochrome c oxidase subunits II','NADH dehydrogenase subunits 2','NADH dehydrogenase subunit2','NADH dehydrogenase subunit 2','NADH-ubiquinone oxidoreductase subunit 2','NADH dehyrogenase subunit 2','NADH dehydrogenase subunit II','NADH  dehydrogenase 2','NADH dehydrogenase 2','nad2-1','nad2-0','NADH dehydrogenase subunit 2','nad2','NADH2','nd2','nadh2','NAD2',
                         'NADH dehydrogenase subunit II','nadh2','ndh2','Nad2','ND-2','ND2 protein','NADH hydrogenase subunit 2','NDII'],
                  'nad3':['NADH3 dehydrogenase subunit 3','NADH dehydogenase subunit 3','ND3','NADH dehydrogenase subunit-3','cytochrome c oxidase subunits III','NADH dehydrogenase subunits 3','NADH dehydrogenase subunit3','NADH dehydrogenase subunit 3','NADH-ubiquinone oxidoreductase subunit 3','NADH dehyrogenase subunit 3','NADH dehydrogenase subunit III','NADH  dehydrogenase 3','NADH dehydrogenase 3','nad3-1','nad3-0','NADH dehydrogenase subunit 3','nad3','NADH3','nd3','nadh3','NAD3',
                         'NADH dehydrogenase subunit III','nadh3','ndh3','Nad3','ND-3','ND3 protein','NADH hydrogenase subunit 3','NDIII'],
                  'nad4':['NADH4 dehydrogenase subunit 4','NADH dehydogenase subunit 4','ND4','NADH dehydrogenase subunit-4','cytochrome c oxidase subunits IV','NADH dehydrogenase subunits 4','NADH dehydrogenase subunit4','NADH dehydrogenase subunit 4','NADH-ubiquinone oxidoreductase subunit 4','NADH dehyrogenase subunit 4','NADH dehydrogenase subunit IV','NADH  dehydrogenase 4','NADH dehydrogenase 4','nad4-1','nad4-0','NADH dehydrogenase subunit 4','nad4','NADH4','nd4','nadh4','NAD4',
                         'NADH dehydrogenase subunit IV','nadh4','ndh4','Nad4','ND-4','ND4 protein','NADH hydrogenase subunit 4','NDIV'],
                  'nad5':['NAHD dehydrogenase subunit 5','NADH5 dehydrogenase subunit 5','NADH dehydogenase subunit 5','ND5','NADH dehydrogenase subunit-5','cytochrome c oxidase subunits V','NADH dehydrogenase subunits 5','NADH dehydrogenase subunit5','NADH dehydrogenase subunit 5','NADH-ubiquinone oxidoreductase subunit 5','NADH dehyrogenase subunit 5','NADH dehydrogenase subunit V','NADH  dehydrogenase 5','NADH dehydrogenase 5','nad5-1','nad5-0','NADH dehydrogenase subunit 5','nad5','NADH5','nd5','nadh5','NAD5',
                         'NADH dehydrogenase subunit V','nadh5','ndh5','Nad5','ND-5','ND5 protein','NADH hydrogenase subunit 5','NDV'],
                  'nad6':['NADH6 dehydrogenase subunit 6','NADH dehydogenase subunit 6','ND6','NADH dehydrogenase subunit-6','cytochrome c oxidase subunits VI','NADH dehydrogenase subunits 6','NADH dehydrogenase subunit6','NADH dehydrogenase subunit 6','NADH-ubiquinone oxidoreductase subunit 6','NADH dehyrogenase subunit 6','NADH dehydrogenase subunit VI','NADH  dehydrogenase 6','NADH dehydrogenase 6','nad6-1','nad6-0','NADH dehydrogenase subunit 6','nad6','NADH6','nd6','nadh6','NAD6',
                         'NADH dehydrogenase subunit VI','nadh6','ndh6','Nad6','ND-6','ND6 protein','NADH hydrogenase subunit 6','NDVI'],
                  'atp8':['ATP synthetase F0 subunit 8','ATP synthasee subunit 8','ATP8','ATPase subunit-8','ATPase subunits 8','ATPase subunit-8','ATPase subunits 8','ATP synthase 8','ATP synthase FO subunit 8','ATP synthase F0 subunit 8','atp8','ATPase8','ATPase 8',
                          'ATP synthase subunit 8','ATP 8','AT8','atpase8','ATP synthetase subunit 8','adenosine triphosphate subunit 8','Atp8','ATPase subunit 8','ATPase sunuint 8'],
                  'atp6':['ATP synthetase F0 subunit 6','ATP synthasee subunit 6','ATP6','ATPase subunit-6','ATPase subunits 6','ATPase subunit-6','ATPase subunits 6','ATP synthase 6','ATP synthase FO subunit 6','ATP synthase F0 subunit 6','atp6','ATPase6','ATPase 6',
                          'ATP synthase subunit 6','ATP 6','AT6','atpase6','ATP synthetase subunit 6','adenosine triphosphate subunit 6','Atp6','ATPase subunit 6','ATPase subuint 6'],
                  'nad4L':['NADH4L dehydrogenase subunit 4L','NADH dehydogenase subunit 4L','ND4L','NADH dehydrogenase subunit-4L','cytochrome c oxidase subunits IV L','NADH dehydrogenase subunits 4L','NADH dehydrogenase subunit4L','NADH dehydrogenase subunit 4 L','NADH-ubiquinone oxidoreductase subunit 4L','NADH dehyrogenase subunit 4L','NADH dehydrogenase subunit IV L','NADH  dehydrogenase 4L','NADH dehydrogenase 4L','NADH dehydrogenase subunit 4L','nad4l','nad4L','NADH4L','nd4L','NADH dehydrognase subunit 4L',
                          'nadh4L','NAD4L','NADH dehydrognase subunit IV L','nadh4l','Nad4L','ND-4L','ND4L protein','NADH hydrogenase subunit 4L','ND4l','NDIVL','NADH dehydrogenase subunit 4l'],
                  'atp9':['ATP9','ATPase subunits 9','ATPase subunit-9','ATPase subunits 9','ATP synthase 9','ATP synthase F0 subunit 9','atp9','ATPase9','ATPase 9','ATP synthase subunit 9','ATP 9',
                          'AT9','atpase9','ATP synthetase subunit 9','adenosine triphosphate subunit 9','Atp9',
                          'ATPase subunit 9','ATPase subuint 9'],
                  'mutS':['msh1','DNA mismatch repair protein MutS','MutS-like protein','putative mismatch repair protein','DNA mismatch repair protein mutS',
                          'DNA mismatch repair protein','MutS-like protein','mismatch repair protein'],
                  'heg':['homing endonuclease','truncated homing endonuclease','homing endonuclease'],
                  'secY':['SecY-independent transporter protein'],
                  'Reph':['replication helicase subunit'],
                  'ORFX':[],
                  'nad7':['ND7'],
                  'nad9':['ND9', 'nad9_2'],
                  'RNAX':[],
                  'nad8':['ND8'],
                  'nad10':['ND10'],
                  'nad11':['ND11'],
                  'cox11':['COX11'],
                  'atp1':['ATP1'],
                  'atp3':['ATP3']}