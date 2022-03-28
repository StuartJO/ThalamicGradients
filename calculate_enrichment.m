function [enrichment,p] = calculate_enrichment(hit_list, top_genes, full_gene_list)

x = sum(ismember(top_genes,hit_list));
n = sum(ismember(hit_list,full_gene_list));

N = length(top_genes);
M = length(full_gene_list);

enrichment = (x/N) / ((n-x) /(M-N));

if enrichment == Inf
    enrichment = 0;
end

p = 1-hygecdf(x-1,M,n,N);