##################

### EdgeR differential expression anlysis 
### use first
###     source("https://bioconductor.org/biocLite.R")
###     biocLite("edgeR")
###

##################

### This program is an example of edgeR when a single data for each time point is available
### control allows to use te different time points o estimate dispersion
### See edgeR user guide for other estimate methods

######
library(edgeR)

rawData = read.table('Expression_Genes.txt', row.names = 1, skip=1)


time = factor(c('t0', 't1', 't2', 't3','t4'))

design = model.matrix(~0+time)
colnames(design) = levels(time)
print(time)
DGEData = DGEList(counts=rawData,group = time)
#
DGEData_control = DGEList(counts=rawData,group=rep(1,each=5))  

DGEData  = calcNormFactors(DGEData)
#
DGEData_control  = calcNormFactors(DGEData_control)
normData3 = cpm(DGEData ,normalized.lib.sizes=FALSE)
normData5 = cpm(DGEData_control,normalized.lib.sizes=FALSE)
#
DGEData_control = estimateGLMCommonDisp(DGEData_control, method="deviance",robust=TRUE, subset=NULL)
DGEData$common.dispersion <- DGEData_control$common.dispersion



fit = glmFit(DGEData, design)
#fit = glmFit(DGEData, design, dispersion=0.1) # dispersion can be assume a priori

comps = makeContrasts(Time0_1 = t1-t0,
                      Time0_2 = t2-t0,
                      Time0_3 = t3-t0,
                      Time0_4 = t4-t0,
                      levels=design)

Time0_1_fit = glmLRT(fit, contrast=comps[,'Time0_1'])
write.table(topTags(Time0_1,n=30000,adjust.method = 'BH', sort.by='p.value'),'Time0_1.de',sep='\t',quote=FALSE)

Time0_2_fit = glmLRT(fit, contrast=comps[,'Time0_2'])
write.table(topTags(Time0_2,n=30000,adjust.method = 'BH', sort.by='p.value'),'Time0_2.de',sep='\t',quote=FALSE)

Time0_3_fit = glmLRT(fit, contrast=comps[,'Time0_3'])
write.table(topTags(Time0_3,n=30000,adjust.method = 'BH', sort.by='p.value'),'Time0_3.de',sep='\t',quote=FALSE)

Time0_4_fit = glmLRT(fit, contrast=comps[,'Time0_4'])
write.table(topTags(Time0_4,n=30000,adjust.method = 'BH', sort.by='p.value'),'Time0_4.de',sep='\t',quote=FALSE)





