library("STRINGdb")
library("biomaRt")

v <- read.delim("~/oskar/results/filtered-variants/StA.tsv")

db <- STRINGdb$new(version="9_1", species=9606, score_threshold=400)

result <- data.frame(gene=character(), stringsAsFactors=F)
for(disease_gene in c("BUB1B", "CEP57")) {
	neighbors <- db$get_neighbors(db$mp(disease_gene))
	neighbors <- sapply(strsplit(neighbors, "9606."), "[[", 2)
	
	# convert back to gene symbol using Ensembl biomart
	mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl") # GRCh37, v75
	neighbors.hugo <- getBM(attributes=c("ensembl_peptide_id", "hgnc_symbol"), filters="ensembl_peptide_id", values=neighbors, mart)
	neighbors.hugo[,disease_gene] <- disease_gene
	
	result <- merge(result, neighbors.hugo[,c("hgnc_symbol", disease_gene)], by.x="gene", by.y="hgnc_symbol", all=T)
}

result$string_1st_degree <- apply(result, 1, function(x) paste(na.omit(x[2:3]),collapse=","))

v <- merge(v, result[,c("gene", "string_1st_degree")], by.x="gene", by.y="gene", all.x=T)

write.table(v, file="~/oskar/results/filtered-variants/StA.string.tsv", col.names=T, row.names=F, sep="\t", quote=F)
