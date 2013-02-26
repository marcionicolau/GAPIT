\name{GAPIT}
\alias{GAPIT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT(Y = NULL, G = NULL, GD = NULL, GM = NULL, KI = NULL, Z = NULL, CV = NULL, CV.Inheritance = NULL, GP = NULL, GK = NULL, group.from = 30, group.to = 1e+06, group.by = 10, DPP = 1e+05, kinship.cluster = "average", kinship.group = "Mean", kinship.algorithm = "VanRaden", bin.from = 10000, bin.to = 10000, bin.by = 10000, inclosure.from = 10, inclosure.to = 10, inclosure.by = 10, SNP.P3D = TRUE, SNP.effect = "Add", SNP.impute = "Middle", PCA.total = 0, SNP.fraction = 1, seed = 123, BINS = 20, SNP.test = TRUE, SNP.MAF = 0, FDR.Rate = 1, SNP.FDR = 1, SNP.permutation = FALSE, SNP.CV = NULL, SNP.robust = "GLM", file.from = 1, file.to = 1, file.total = NULL, file.fragment = 99999, file.path = NULL, file.G = NULL, file.Ext.G = NULL, file.GD = NULL, file.GM = NULL, file.Ext.GD = NULL, file.Ext.GM = NULL, ngrid = 100, llim = -10, ulim = 10, esp = 1e-10, LD.chromosome = NULL, LD.location = NULL, LD.range = NULL, sangwich.top = NULL, sangwich.bottom = NULL, QC = TRUE, GTindex = NULL, LD = 0.01, file.output = TRUE, cutOff = 0.01, Model.selection = FALSE, output.numerical = FALSE, output.hapmap = FALSE, Create.indicator = FALSE, QTN = NULL, QTN.round = 1, QTN.limit = 0, QTN.update = TRUE, QTN.method = "Penalty", Major.allele.zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{G}{
%%     ~~Describe \code{G} here~~
}
  \item{GD}{
%%     ~~Describe \code{GD} here~~
}
  \item{GM}{
%%     ~~Describe \code{GM} here~~
}
  \item{KI}{
%%     ~~Describe \code{KI} here~~
}
  \item{Z}{
%%     ~~Describe \code{Z} here~~
}
  \item{CV}{
%%     ~~Describe \code{CV} here~~
}
  \item{CV.Inheritance}{
%%     ~~Describe \code{CV.Inheritance} here~~
}
  \item{GP}{
%%     ~~Describe \code{GP} here~~
}
  \item{GK}{
%%     ~~Describe \code{GK} here~~
}
  \item{group.from}{
%%     ~~Describe \code{group.from} here~~
}
  \item{group.to}{
%%     ~~Describe \code{group.to} here~~
}
  \item{group.by}{
%%     ~~Describe \code{group.by} here~~
}
  \item{DPP}{
%%     ~~Describe \code{DPP} here~~
}
  \item{kinship.cluster}{
%%     ~~Describe \code{kinship.cluster} here~~
}
  \item{kinship.group}{
%%     ~~Describe \code{kinship.group} here~~
}
  \item{kinship.algorithm}{
%%     ~~Describe \code{kinship.algorithm} here~~
}
  \item{bin.from}{
%%     ~~Describe \code{bin.from} here~~
}
  \item{bin.to}{
%%     ~~Describe \code{bin.to} here~~
}
  \item{bin.by}{
%%     ~~Describe \code{bin.by} here~~
}
  \item{inclosure.from}{
%%     ~~Describe \code{inclosure.from} here~~
}
  \item{inclosure.to}{
%%     ~~Describe \code{inclosure.to} here~~
}
  \item{inclosure.by}{
%%     ~~Describe \code{inclosure.by} here~~
}
  \item{SNP.P3D}{
%%     ~~Describe \code{SNP.P3D} here~~
}
  \item{SNP.effect}{
%%     ~~Describe \code{SNP.effect} here~~
}
  \item{SNP.impute}{
%%     ~~Describe \code{SNP.impute} here~~
}
  \item{PCA.total}{
%%     ~~Describe \code{PCA.total} here~~
}
  \item{SNP.fraction}{
%%     ~~Describe \code{SNP.fraction} here~~
}
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{BINS}{
%%     ~~Describe \code{BINS} here~~
}
  \item{SNP.test}{
%%     ~~Describe \code{SNP.test} here~~
}
  \item{SNP.MAF}{
%%     ~~Describe \code{SNP.MAF} here~~
}
  \item{FDR.Rate}{
%%     ~~Describe \code{FDR.Rate} here~~
}
  \item{SNP.FDR}{
%%     ~~Describe \code{SNP.FDR} here~~
}
  \item{SNP.permutation}{
%%     ~~Describe \code{SNP.permutation} here~~
}
  \item{SNP.CV}{
%%     ~~Describe \code{SNP.CV} here~~
}
  \item{SNP.robust}{
%%     ~~Describe \code{SNP.robust} here~~
}
  \item{file.from}{
%%     ~~Describe \code{file.from} here~~
}
  \item{file.to}{
%%     ~~Describe \code{file.to} here~~
}
  \item{file.total}{
%%     ~~Describe \code{file.total} here~~
}
  \item{file.fragment}{
%%     ~~Describe \code{file.fragment} here~~
}
  \item{file.path}{
%%     ~~Describe \code{file.path} here~~
}
  \item{file.G}{
%%     ~~Describe \code{file.G} here~~
}
  \item{file.Ext.G}{
%%     ~~Describe \code{file.Ext.G} here~~
}
  \item{file.GD}{
%%     ~~Describe \code{file.GD} here~~
}
  \item{file.GM}{
%%     ~~Describe \code{file.GM} here~~
}
  \item{file.Ext.GD}{
%%     ~~Describe \code{file.Ext.GD} here~~
}
  \item{file.Ext.GM}{
%%     ~~Describe \code{file.Ext.GM} here~~
}
  \item{ngrid}{
%%     ~~Describe \code{ngrid} here~~
}
  \item{llim}{
%%     ~~Describe \code{llim} here~~
}
  \item{ulim}{
%%     ~~Describe \code{ulim} here~~
}
  \item{esp}{
%%     ~~Describe \code{esp} here~~
}
  \item{LD.chromosome}{
%%     ~~Describe \code{LD.chromosome} here~~
}
  \item{LD.location}{
%%     ~~Describe \code{LD.location} here~~
}
  \item{LD.range}{
%%     ~~Describe \code{LD.range} here~~
}
  \item{sangwich.top}{
%%     ~~Describe \code{sangwich.top} here~~
}
  \item{sangwich.bottom}{
%%     ~~Describe \code{sangwich.bottom} here~~
}
  \item{QC}{
%%     ~~Describe \code{QC} here~~
}
  \item{GTindex}{
%%     ~~Describe \code{GTindex} here~~
}
  \item{LD}{
%%     ~~Describe \code{LD} here~~
}
  \item{file.output}{
%%     ~~Describe \code{file.output} here~~
}
  \item{cutOff}{
%%     ~~Describe \code{cutOff} here~~
}
  \item{Model.selection}{
%%     ~~Describe \code{Model.selection} here~~
}
  \item{output.numerical}{
%%     ~~Describe \code{output.numerical} here~~
}
  \item{output.hapmap}{
%%     ~~Describe \code{output.hapmap} here~~
}
  \item{Create.indicator}{
%%     ~~Describe \code{Create.indicator} here~~
}
  \item{QTN}{
%%     ~~Describe \code{QTN} here~~
}
  \item{QTN.round}{
%%     ~~Describe \code{QTN.round} here~~
}
  \item{QTN.limit}{
%%     ~~Describe \code{QTN.limit} here~~
}
  \item{QTN.update}{
%%     ~~Describe \code{QTN.update} here~~
}
  \item{QTN.method}{
%%     ~~Describe \code{QTN.method} here~~
}
  \item{Major.allele.zero}{
%%     ~~Describe \code{Major.allele.zero} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y = NULL, G = NULL, GD = NULL, GM = NULL, KI = NULL, 
    Z = NULL, CV = NULL, CV.Inheritance = NULL, GP = NULL, GK = NULL, 
    group.from = 30, group.to = 1e+06, group.by = 10, DPP = 1e+05, 
    kinship.cluster = "average", kinship.group = "Mean", kinship.algorithm = "VanRaden", 
    bin.from = 10000, bin.to = 10000, bin.by = 10000, inclosure.from = 10, 
    inclosure.to = 10, inclosure.by = 10, SNP.P3D = TRUE, SNP.effect = "Add", 
    SNP.impute = "Middle", PCA.total = 0, SNP.fraction = 1, seed = 123, 
    BINS = 20, SNP.test = TRUE, SNP.MAF = 0, FDR.Rate = 1, SNP.FDR = 1, 
    SNP.permutation = FALSE, SNP.CV = NULL, SNP.robust = "GLM", 
    file.from = 1, file.to = 1, file.total = NULL, file.fragment = 99999, 
    file.path = NULL, file.G = NULL, file.Ext.G = NULL, file.GD = NULL, 
    file.GM = NULL, file.Ext.GD = NULL, file.Ext.GM = NULL, ngrid = 100, 
    llim = -10, ulim = 10, esp = 1e-10, LD.chromosome = NULL, 
    LD.location = NULL, LD.range = NULL, sangwich.top = NULL, 
    sangwich.bottom = NULL, QC = TRUE, GTindex = NULL, LD = 0.01, 
    file.output = TRUE, cutOff = 0.01, Model.selection = FALSE, 
    output.numerical = FALSE, output.hapmap = FALSE, Create.indicator = FALSE, 
    QTN = NULL, QTN.round = 1, QTN.limit = 0, QTN.update = TRUE, 
    QTN.method = "Penalty", Major.allele.zero = FALSE) 
{
    print("--------------------- Welcome to GAPIT ----------------------------")
    echo = TRUE
    GAPIT.Version = GAPIT.0000()
    Timmer = GAPIT.Timmer(Infor = "GAPIT")
    Memory = GAPIT.Memory(Infor = "GAPIT")
    myGenotype <- GAPIT.Genotype(G = G, GD = GD, GM = GM, KI = KI, 
        kinship.algorithm = kinship.algorithm, PCA.total = PCA.total, 
        SNP.fraction = SNP.fraction, SNP.test = SNP.test, file.path = file.path, 
        file.from = file.from, file.to = file.to, file.total = file.total, 
        file.fragment = file.fragment, file.G = file.G, file.Ext.G = file.Ext.G, 
        file.GD = file.GD, file.GM = file.GM, file.Ext.GD = file.Ext.GD, 
        file.Ext.GM = file.Ext.GM, SNP.MAF = SNP.MAF, FDR.Rate = FDR.Rate, 
        SNP.FDR = SNP.FDR, SNP.effect = SNP.effect, SNP.impute = SNP.impute, 
        LD.chromosome = LD.chromosome, LD.location = LD.location, 
        LD.range = LD.range, GP = GP, GK = GK, bin.size = NULL, 
        inclosure.size = NULL, Timmer = Timmer, Memory = Memory, 
        sangwich.top = sangwich.top, sangwich.bottom = sangwich.bottom, 
        GTindex = NULL, file.output = file.output, Create.indicator = Create.indicator, 
        Major.allele.zero = Major.allele.zero)
    Timmer = myGenotype$Timmer
    Memory = myGenotype$Memory
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Genotype for all")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Genotype for all")
    KI = myGenotype$KI
    PC = myGenotype$PC
    genoFormat = myGenotype$genoFormat
    hasGenotype = myGenotype$hasGenotype
    byFile = myGenotype$byFile
    fullGD = myGenotype$fullGD
    GD = myGenotype$GD
    GI = myGenotype$GI
    GT = myGenotype$GT
    G = myGenotype$G
    rownames(GD) = GT
    colnames(GD) = GI[, 1]
    if (output.numerical) 
        write.table(GD, "GAPIT.Genotype.Numerical.txt", quote = FALSE, 
            sep = "\t", row.names = TRUE, col.names = NA)
    if (output.hapmap) 
        write.table(myGenotype$G, "GAPIT.Genotype.hmp.txt", quote = FALSE, 
            sep = "\t", row.names = FALSE, col.names = FALSE)
    if (is.null(Y) & is.null(GP)) 
        return(list(GWAS = NULL, GPS = NULL, Pred = NULL, compression = NULL, 
            kinship.optimum = NULL, kinship = myGenotype$KI, 
            PCA = myGenotype$PC, GD = data.frame(cbind(as.data.frame(GT), 
                as.data.frame(GD))), GI = GI, G = myGenotype$G))
    if (is.null(Y)) 
        return(list(GWAS = NULL, GPS = NULL, Pred = NULL, compression = NULL, 
            kinship.optimum = NULL, kinship = myGenotype$KI, 
            PCA = myGenotype$PC, GD = data.frame(cbind(as.date.frame(GT), 
                as.data.frame(GD))), Gi = GI, G = myGenotype$G))
    rm(myGenotype)
    gc()
    print("--------------------Processing traits----------------------------------")
    if (!is.null(Y)) {
        print("Phenotype provided!")
        if (ncol(Y) < 2) 
            stop("Phenotype should have taxa name and one trait at least. Please correct phenotype file!")
        for (trait in 2:ncol(Y)) {
            print(paste("Processing trait: ", colnames(Y)[trait], 
                sep = ""))
            gapitMain <- GAPIT.Main(Y = Y[, c(1, trait)], G = G, 
                GD = GD, GM = GM, KI = KI, Z = Z, CV = CV, CV.Inheritance = CV.Inheritance, 
                GP = GP, GK = GK, SNP.P3D = SNP.P3D, kinship.algorithm = kinship.algorithm, 
                bin.from = bin.from, bin.to = bin.to, bin.by = bin.by, 
                inclosure.from = inclosure.from, inclosure.to = inclosure.to, 
                inclosure.by = inclosure.by, group.from = group.from, 
                group.to = group.to, group.by = group.by, kinship.cluster = kinship.cluster, 
                kinship.group = kinship.group, name.of.trait = colnames(Y)[trait], 
                file.path = file.path, file.from = file.from, 
                file.to = file.to, file.total = file.total, file.fragment = file.fragment, 
                file.G = file.G, file.Ext.G = file.Ext.G, file.GD = file.GD, 
                file.GM = file.GM, file.Ext.GD = file.Ext.GD, 
                file.Ext.GM = file.Ext.GM, SNP.MAF = SNP.MAF, 
                FDR.Rate = FDR.Rate, SNP.FDR = SNP.FDR, SNP.effect = SNP.effect, 
                SNP.impute = SNP.impute, PCA.total = PCA.total, 
                GAPIT.Version = GAPIT.Version, GT = GT, SNP.fraction = SNP.fraction, 
                seed = seed, BINS = BINS, SNP.test = SNP.test, 
                DPP = DPP, SNP.permutation = SNP.permutation, 
                LD.chromosome = LD.chromosome, LD.location = LD.location, 
                LD.range = LD.range, SNP.CV = SNP.CV, SNP.robust = SNP.robust, 
                genoFormat = genoFormat, hasGenotype = hasGenotype, 
                byFile = byFile, fullGD = fullGD, PC = PC, GI = GI, 
                Timmer = Timmer, Memory = Memory, sangwich.top = sangwich.top, 
                sangwich.bottom = sangwich.bottom, QC = QC, GTindex = GTindex, 
                LD = LD, file.output = file.output, cutOff = cutOff, 
                Model.selection = Model.selection, Create.indicator = Create.indicator, 
                QTN = QTN, QTN.round = QTN.round, QTN.limit = QTN.limit, 
                QTN.update = QTN.update, QTN.method = QTN.method, 
                Major.allele.zero = Major.allele.zero)
        }
        if (ncol(Y > 2) & file.output) {
            Timmer = gapitMain$Timmer
            Memory = gapitMain$Memory
            file = paste("GAPIT.", "All", ".Timming.csv", sep = "")
            write.table(Timmer, file, quote = FALSE, sep = ",", 
                row.names = FALSE, col.names = TRUE)
            file = paste("GAPIT.", "All", ".Memory.Stage.csv", 
                sep = "")
            write.table(Memory, file, quote = FALSE, sep = ",", 
                row.names = FALSE, col.names = TRUE)
        }
        if (ncol(Y) == 2) {
            h2 = as.matrix(as.numeric(as.vector(gapitMain$Compression[, 
                5]))/(as.numeric(as.vector(gapitMain$Compression[, 
                5])) + as.numeric(as.vector(gapitMain$Compression[, 
                6]))), length(gapitMain$Compression[, 6]), 1)
            colnames(h2) = c("Heritability")
            print("GAPIT accomplished successfully for single trait. Results are saved. GWAS and GPS are returned!")
            return(list(QTN = gapitMain$QTN, GWAS = gapitMain$GWAS, 
                GPS = gapitMain$GPS, Pred = gapitMain$Pred, compression = as.data.frame(cbind(gapitMain$Compression, 
                  h2)), kinship.optimum = gapitMain$kinship.optimum, 
                kinship = gapitMain$kinship, PCA = gapitMain$PC))
        }
        else {
            print("GAPIT accomplished successfully for multiple traits. Results are saved")
            return(list(GWAS = NULL, GPS = NULL, Pred = NULL, 
                compression = NULL, kinship.optimum = NULL, kinship = gapitMain$KI, 
                PCA = gapitMain$PC))
        }
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
