import csv
import GEOparse
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import pearsonr
from sklearn.decomposition import PCA
from sklearn.feature_selection import VarianceThreshold, SelectKBest
from sklearn.feature_selection import chi2
from sklearn.preprocessing import StandardScaler
import seaborn as sns


SMALL_SIZE = 2
MEDIUM_SIZE = 10
BIGGER_SIZE = 12

plt.rc('font', size=SMALL_SIZE)  # controls default text sizes
plt.rc('axes', titlesize=SMALL_SIZE)  # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)  # fontsize of the x and y labels
plt.rc('xtick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('ytick', labelsize=SMALL_SIZE)  # fontsize of the tick labels
plt.rc('legend', fontsize=SMALL_SIZE)  # legend fontsize
plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title
sns.set_theme()
# Khoda nane babaye sklearn ro biamorze vaghean !!!

control_samples = list()
test_samples = list()
gse = GEOparse.get_GEO(filepath="GSE48558_family.soft.gz", destdir="GSE48558")

exprs = []
gsmNames = []
gene_idx = []
metadata = {}


def read_expressions(gsm, exprs):
    if len(gsm.table) > 0:
        tmp = gsm.table['VALUE']
        tmp.index = gsm.table['ID_REF']
        gsmNames.append(name)
        if len(exprs) == 0:
            exprs = tmp.to_frame()
        else:
            exprs = pd.concat([exprs, tmp.to_frame()], axis=1)
        return exprs


file = open("Test/GSM1180750.txt", "r")
file = (file.read())
lines = file.split("\n")
for line in lines[1:-1]:
    line = line.split()
    gene_idx.append(line[0])

for name, gsm in gse.gsms.items():
    name = name.strip()
    sample = gse.gsms[name]
    if str(sample.metadata['source_name_ch1']) == str(['AML Patient']):
        # sample.table.to_csv('Test/' + name + '.txt', index=None, sep='\t', mode='w', )
        test_samples.append(sample)
        exprs = read_expressions(gsm, exprs)

    if str(sample.metadata['characteristics_ch1']) == str(['phenotype: Normal']):
        # sample.table.to_csv('Control/' + name + '.txt', index=None, sep='\t', mode='w', )
        control_samples.append(sample)
        exprs = read_expressions(gsm, exprs)

all_samples = control_samples + test_samples

exprs.columns = gsmNames

# Plot boxplot of expression data
with PdfPages('GSE_boxplot.pdf') as pdf:
    plt.boxplot(exprs, showfliers=False)
    plt.title(' Box Plot')
    # pdf.savefig()
    plt.savefig('boxplot.png')
    plt.close()
    '''try:
        plt.boxplot(numpy.log2(exprs).transpose(), showfliers=False)
        plt.title('log2(' + 'Box Plot' + ')')
        pdf.savefig()
        plt.close()
    except:
        pass'''

exprs_norm = (exprs - exprs.mean()) / exprs.std()

plt.boxplot(exprs_norm, showfliers=False)
plt.title(' Box Plot Normalized')
plt.savefig('boxplot_norm.png')
plt.close()

########################################################################################################################
genes = []
genes_dict = dict()


class Gene:
    def __init__(self, ID, adj_P_Val, P_Value, t, B, logFC, Gene_symbol, Gene_title):
        self.id = ID
        self.adj_P_val = adj_P_Val
        self.P_value = P_Value
        self.t = t
        self.B = B
        self.logFC = logFC
        self.Gene_symbol = Gene_symbol
        self.Gene_title = Gene_title


with open('GSE48558.top.table2.tsv') as tsv_file:
    csv_reader = csv.reader(tsv_file, delimiter="\t")
    line_count = 0
    for row in csv_reader:
        if line_count != 0:
            gene = Gene(int(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4]), float(row[5]),
                        str(row[6]), str(row[7]))
            genes.append(gene)
            genes_dict[gene.id] = gene
        line_count += 1

    print(f'Processed {line_count - 1} lines of data.')

trans = exprs_norm.transpose()

n_components = 50
# pca = PCA(n_components=0.95, svd_solver='full')
pca = PCA(n_components)
fitted = pca.fit(trans)
pca_exprs = pd.DataFrame(pca.transform(trans), columns=['PCA%i' % i for i in range(n_components)], index=trans.index)
print(pca_exprs)
n_pcs = pca.components_.shape[0]
###############################################################################################
# get the index of the most important feature on EACH component i.e. largest absolute value
# using LIST COMPREHENSION HERE
most_important = [np.abs(pca.components_[i]).argmax() for i in range(n_pcs)]

initial_feature_names = gene_idx

# get the names
most_important_names = [initial_feature_names[most_important[i]] for i in range(n_pcs)]
most_important_gene_names = [genes_dict[int(i)].Gene_symbol for i in most_important_names]
print(most_important_gene_names)
# using LIST COMPREHENSION HERE AGAIN
pca_dic = {'PC{}'.format(i + 1): most_important_names[i] for i in range(n_pcs)}
# build the dataframe
df = pd.DataFrame(sorted(pca_dic.items()))
###############################################################################################
heat_exprs = pd.DataFrame()
pca_trans = pca_exprs.transpose()

for name, gsm in gse.gsms.items():
    name = name.strip()
    sample = gse.gsms[name]
    if sample in control_samples:
        col_id = name + '_' + str(sample.metadata['source_name_ch1'])[2:-2]
        heat_exprs[col_id] = pca_trans[name]
    elif sample in test_samples:
        col_id = name + '_' + str('AML')
        heat_exprs[col_id] = pca_trans[name]

corr = heat_exprs.corr()
'''
pv_df = pd.DataFrame()
group = list()
i = 0
for name, gsm in gse.gsms.items():
    name = name.strip()
    sample = gse.gsms[name]
    if sample in control_samples and str(sample.metadata['source_name_ch1'])[2:-2] == 'CD34+HSPC':
        col_id = name + '_' + str(sample.metadata['source_name_ch1'])[2:-2]
        pv_df[col_id] = pca_trans[name]
        group.append('Normal')
    elif sample in test_samples:
        col_id = name + '_' + str('AML')
        pv_df[col_id] = pca_trans[name]
        group.append('AML')
pv_df = pv_df.transpose()
pv_df['group'] = group


grouped = pv_df.groupby('group')
'''


def get_pca_attribute(pca_id):
    return most_important_gene_names[i]


gene_set = set()

gene_list = open('gene_list.txt', 'w')
for gene in gene_set:
    gene_list.write(gene + '\n')

'''
ax = sns.heatmap(data=corr, xticklabels=corr.columns, yticklabels=corr.columns, vmax=1, vmin = -1,
                 cmap = 'Reds', linewidths=0.1, linecolor='Black')
#plt.show()
for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
             ax.get_xticklabels() + ax.get_yticklabels()):
    item.set_fontsize(4)
plt.savefig("foo.png", bbox_inches='tight')

ax = sns.clustermap(data=corr, xticklabels=corr.columns, yticklabels=corr.columns, vmax=1, vmin = -1,
                 cmap = 'PiYG', linewidths=0.1, linecolor='Black')
#plt.show()
plt.savefig("cluster.png", bbox_inches='tight')
file = open('gene_symbols.txt', 'w')
for col in pca_exprs.columns:
    print(col)
    gene_id = int(gene_idx[col])
    file.write(genes_dict[gene_id].Gene_symbol + '\t' + genes_dict[gene_id].Gene_title + '\n')
    exprs_red[gene_id] = trans[gene_id]
file.close()

# print(exprs_red)'''

#######################DUMP LAND
'''
selector = VarianceThreshold()
selector.fit(trans)
vars = selector.variances_
vars = list(vars)
vars.sort()
vars.reverse()
#print(vars)

threshold_n = 2
sel = VarianceThreshold(threshold=threshold_n)
selected = sel.fit(trans)
idx = np.where(selected.variances_ > threshold_n)[0]
file = open('gene_symbols.txt', 'w')
for i in idx:
    gene_id = int(gene_idx[i])
    file.write(genes_dict[gene_id].Gene_symbol + '\t' + genes_dict[gene_id].Gene_title + '\n')
    exprs_red[gene_id] = trans[gene_id]
file.close()
'''
gene_up = list()
gene_down = list()
for gene in genes:
    if gene.adj_P_val < 0.05 and gene.logFC < -1:
        gene_down.append(gene.Gene_symbol.split('///')[0])
    elif gene.adj_P_val < 0.05 and gene.logFC > 1:
        gene_up.append(gene.Gene_symbol.split('///')[0])
aml_up_gene = open('aml_up_gene.txt', 'w')
for i in gene_up:
    aml_up_gene.write(i + '\n')
aml_down_gene = open('aml_down_gene.txt', 'w')
for i in gene_down:
    aml_down_gene.write(i + '\n')
