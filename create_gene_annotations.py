import pandas as pd

biomart_genes = pd.read_csv('/sc/service/analysis/tmp/rat_neuron_clustering_analysis/reference_data/biomart_rat_genes.csv', index_col=0)

biomart_genes['mitochondrial'] = biomart_genes['Chromosome/scaffold name'].str.contains('mt', case=False, na=False)
biomart_genes['hemoglobin'] = biomart_genes['Gene description'].str.contains('hemoglobin', case=False, na=False)
biomart_genes['sex_linked'] = (biomart_genes['Chromosome/scaffold name'].str.match('y', case=False, na=False)) | (biomart_genes['Gene name'].str.contains('xist', case=False, na=False))
biomart_genes['canonical_neuronal'] = biomart_genes['Gene name'].str.match(pat='(Snap25)|(Syp)|(Tubb3)|(Elavl2)', case=False, na=False)
biomart_genes['glutamatergic'] = biomart_genes['Gene name'].str.match(pat='(Slc17a6)|(Slc17a7)|(Grin1)', case=False, na=False)

biomart_genes.to_csv('/sc/service/analysis/tmp/rat_neuron_clustering_analysis/reference_data/biomart_rat_genes.csv')