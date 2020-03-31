library(readr)
library(tximport)
library(rhdf5)
library(gdata)
library(readxl)
library(openxlsx)
library(CoreGx)
library(affy)
library(Biobase)
library(devtools)
install_url("http://mbni.org/customcdf/16.0.0/ensg.download/hgu133plus2hsensgcdf_16.0.0.tar.gz")
library(hgu133plus2hsensgcdf)

#read in sample info
message("Read sample information")
sampleinfo <- read.csv("/pfs/downAnnotations/CCLE_sample_info_file_2012-10-18.txt", sep="\t")
sampleinfo[sampleinfo == "" | sampleinfo == " "] <- NA
sampleinfo <- data.frame("cellid"=as.character(sampleinfo[ , "Cell.line.primary.name"]), sampleinfo, check.names=FALSE)
## remove duplicated cell line hybridization
## only the most recent experiment (as determine by hyridization date or CEL file timestamp) will be kept for each cell line
# sampleinfo <- sampleinfo[!duplicated(sampleinfo[ , "cellid"]), ,drop=FALSE]
#     sampleinfo[ , "cellid"] <- as.character(sampleinfo[ , "cellid"])
rownames(sampleinfo) <- as.character(sampleinfo[ , "cellid"])

`celfileChip` <-
  function (filename) {
    h <- affyio::read.celfile.header(filename, info="full")
    return(as.character(h$cdfName))
  }

`celfileDateHour` <-
function (filename) {
	h <- affyio::read.celfile.header(filename, info="full")
	#ddate <- grep("/", strsplit(h$DatHeader, " ")[[1]], value=TRUE)
	#ddate <- strsplit(ddate, split="/")[[1]]
	#CC <- ifelse(substr(ddate[3],1,1)=="9", "19", "20")
	if(length(h$ScanDate) > 0) {
	    h$ScanDate <- gsub(pattern="T", replacement=" ", x=h$ScanDate)
	    ddate <- strsplit(h$ScanDate, " ")[[1]]
    } else { ddate <- rep(NA, 2)}
    names(ddate) <- c("day", "hour")
	return(ddate)
}

load("/pfs/downloadCCLE_CELArray/celfile.timestamp.RData")
message("Keeping all replicates")
celfn <- list.celfiles(file.path("/pfs/downloadCCLE_CELArray"), full.names=TRUE)
celfns <- list.celfiles(file.path("/pfs/downloadCCLE_CELArray"), full.names=FALSE)
print(head(celfn))
print(head(celfns))
## experiments' names
names(celfn) <- names(celfns) <- gsub(".CEL.gz", "", celfns)
## chip type and date
chipt <- sapply(celfn, celfileChip)
chipd <- t(sapply(celfn, celfileDateHour))

## keep only the CEL files present in sampleinfo
fn <- sampleinfo[ , "Expression.arrays"]
myx <- intersect(fn[!is.na(fn)], names(celfns))
celfn <- celfn[myx]
celfns <- celfns[myx]
chipt <- chipt[myx]
chipd <- chipd[myx, , drop=FALSE]
celfile.timestamp <- celfile.timestamp[paste(myx, ".CEL", sep=""), , drop=FALSE]
sampleinfo <- sampleinfo[match(myx, sampleinfo[ , "Expression.arrays"]), ,drop=FALSE]
sampleinfo <- data.frame("samplename"=names(celfns), "filename"=celfns, "chiptype"=chipt, "hybridization.date"=chipd[ , "day"], "hybridization.hour"=chipd[ , "hour"], "file.day"=celfile.timestamp[ , "file.day"], "file.hour"=celfile.timestamp[ , "file.hour"], "batch"=NA, sampleinfo, check.names=FALSE)
rownames(sampleinfo) <- toupper(names(celfns))

##annotate genes/probes using Ensembl v99 (Jan 2020)
load("/pfs/downAnnotations/Ensembl.v99.annotation.RData")
#annot <- read.csv("/pfs/downAnnotations/annot_ensembl_all_genes.csv", stringsAsFactors=FALSE, check.names=FALSE, header=TRUE, row.names=1)

#create brainarray eset
eset <- just.rma(filenames=celfn, verbose=TRUE, cdfname="hgu133plus2hsensgcdf")
pData(eset) <- as.data.frame(sampleinfo[match(toupper(gsub("[.]CEL[.]gz$", "", rownames(pData(eset)))), rownames(sampleinfo)), , drop=FALSE])
colnames(exprs(eset)) <- rownames(pData(eset)) <- toupper(gsub("[.]CEL[.]gz$", "", colnames(exprs(eset))))
controls <- rownames(exprs(eset))[grep("AFFX", rownames(exprs(eset)))]
fData(eset) <- fData(eset)[which(!rownames(fData(eset)) %in% controls), , drop=FALSE]
exprs_sub <- exprs(eset)[which(!rownames(exprs(eset)) %in% controls), , drop=FALSE]
new_eset <- ExpressionSet(assayData = exprs_sub, phenoData = eset@phenoData, featureData = eset@featureData)
eset <- new_eset
ensemblIds <- sapply(strsplit(rownames(exprs(eset)), "_"), function (x) { return (x[[1]]) }) 
fData(eset) <- data.frame("Probe"=rownames(exprs(eset)), 
                                    "EnsemblGeneId"=ensemblIds,
                                    "Symbol"=features_gene[ensemblIds, "gene_name"],
                                    "GeneBioType"=features_gene[ensemblIds, "gene_biotype"])

rownames(fData(eset)) <- rownames(exprs(eset))
pData(eset)[ , "batchid"] <- NA
annotation(eset) <- "rna"
experimentData(eset)@preprocessing <- list(Normalisation=list(name="rma", package="affy", version=as.character(packageVersion("affy")))) 
#save microarray eset
save(eset, file="/pfs/out/ccle_ge_brainarray_rma.RData")



