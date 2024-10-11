library(limma)
library(edgeR)

count <- read.table("/Users/daehwa/Library/CloudStorage/OneDrive-Personal/Junlab/Projects/Adipocyte/Results/Tables/CPM_RD_calculation/v20240213_RD-edgeR/P-R-counts_afterFilt.tsv", header = T, sep = '\t')
head(count)

y <- DGEList(counts=count[3:20], genes=count[1:2])
y$samples

y$samples$lib.size <- colSums(y$counts)

rownames(y$counts) <- rownames(y$genes) <- y$genes$gene_id
head(y$counts)

## TMM normalization
y <- normLibSizes(y)
y$samples
# plotMDS(y)

LIB      <- factor(substring(colnames(count)[3:20],1,1))
Treat    <- factor(paste0(substring(colnames(count)[3:20],2,2), 'd'))
Rep      <- factor(substring(colnames(count)[3:20],3,3))
targets  <- data.frame(Sample=colnames(y), LIB, Treat, Rep)
targets

Group <- factor(paste(targets$LIB, targets$Treat, sep="."))
cbind(targets,Group=Group)

design <- model.matrix(~0+Group)
colnames(design) <- levels(Group)
design

my.contrasts <- makeContrasts(
                        PvsR.4d_vs_0d = (P.4d - P.0d) - (R.4d - R.0d),
                        PvsR.8d_vs_0d = (P.8d - P.0d) - (R.8d - R.0d),
                        levels=design)

#estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
y$common.dispersion
# plotBCV(y)

fit <- glmQLFit(y, design)


# qlf <- glmQLFTest(fit, contrast=my.contrasts[,paste0("R.Tc.30m_vs_Ctl")])
# # topTags(qlf)

# plotMD(qlf)
# abline(h=c(0), col="red")

for (Cd in c('4d_vs_0d','8d_vs_0d')) {
        qlf <- glmQLFTest(fit, contrast=my.contrasts[,paste0("PvsR.",Cd)])
        # topTags(qlf)

        # plotMD(qlf)
        # abline(h=c(0), col="red")

        FDR <-p.adjust(qlf$table$PValue, method='fdr')
        # result <- cbind(qlf$genes, qlf$table, FDR, cpm(y))
        result <- cbind(qlf$genes, qlf$table, FDR)

        write.table(result, paste0("/Users/daehwa/Library/CloudStorage/OneDrive-Personal/Junlab/Projects/Adipocyte/Results/Tables/CPM_RD_calculation/v20240213_RD-edgeR/adi_PvsR_",Cd,".tsv"), sep = '\t', col.names = T, row.names = F, quote = FALSE)
}
