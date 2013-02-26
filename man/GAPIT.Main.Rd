\name{GAPIT.Main}
\alias{GAPIT.Main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Main(Y, G = NULL, GD = NULL, GM = NULL, KI = NULL, Z = NULL, CV = NULL, CV.Inheritance = NULL, SNP.P3D = TRUE, GP = NULL, GK = NULL, group.from = 1e+06, group.to = 1, group.by = 10, kinship.cluster = "average", kinship.group = "Mean", kinship.algorithm = NULL, DPP = 50000, ngrid = 100, llin = -10, ulim = 10, esp = 1e-10, file.path = NULL, file.from = NULL, file.to = NULL, file.total = NULL, file.fragment = 512, file.G = NULL, file.Ext.G = NULL, file.GD = NULL, file.GM = NULL, file.Ext.GD = NULL, file.Ext.GM = NULL, SNP.MAF = 0, FDR.Rate = 1, SNP.FDR = 1, SNP.effect = "Add", SNP.impute = "Middle", PCA.total = 0, GAPIT.Version = GAPIT.Version, name.of.trait, GT = NULL, SNP.fraction = 1, seed = 123, BINS = 20, SNP.test = TRUE, SNP.robust = "FaST", LD.chromosome = NULL, LD.location = NULL, LD.range = NULL, bin.from = 10000, bin.to = 5e+06, bin.by = 1000, inclosure.from = 10, inclosure.to = 1000, inclosure.by = 10, SNP.permutation = FALSE, SNP.CV = NULL, genoFormat = NULL, hasGenotype = NULL, byFile = NULL, fullGD = NULL, PC = NULL, GI = NULL, Timmer = NULL, Memory = NULL, sangwich.top = NULL, sangwich.bottom = NULL, QC = TRUE, GTindex = NULL, LD = 0.05, file.output = TRUE, cutOff = 0.01, Model.selection = FALSE, Create.indicator = FALSE, QTN = NULL, QTN.round = 1, QTN.limit = 0, QTN.update = TRUE, QTN.method = "Penalty", Major.allele.zero = FALSE)
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
  \item{SNP.P3D}{
%%     ~~Describe \code{SNP.P3D} here~~
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
  \item{kinship.cluster}{
%%     ~~Describe \code{kinship.cluster} here~~
}
  \item{kinship.group}{
%%     ~~Describe \code{kinship.group} here~~
}
  \item{kinship.algorithm}{
%%     ~~Describe \code{kinship.algorithm} here~~
}
  \item{DPP}{
%%     ~~Describe \code{DPP} here~~
}
  \item{ngrid}{
%%     ~~Describe \code{ngrid} here~~
}
  \item{llin}{
%%     ~~Describe \code{llin} here~~
}
  \item{ulim}{
%%     ~~Describe \code{ulim} here~~
}
  \item{esp}{
%%     ~~Describe \code{esp} here~~
}
  \item{file.path}{
%%     ~~Describe \code{file.path} here~~
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
  \item{SNP.MAF}{
%%     ~~Describe \code{SNP.MAF} here~~
}
  \item{FDR.Rate}{
%%     ~~Describe \code{FDR.Rate} here~~
}
  \item{SNP.FDR}{
%%     ~~Describe \code{SNP.FDR} here~~
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
  \item{GAPIT.Version}{
%%     ~~Describe \code{GAPIT.Version} here~~
}
  \item{name.of.trait}{
%%     ~~Describe \code{name.of.trait} here~~
}
  \item{GT}{
%%     ~~Describe \code{GT} here~~
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
  \item{SNP.robust}{
%%     ~~Describe \code{SNP.robust} here~~
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
  \item{SNP.permutation}{
%%     ~~Describe \code{SNP.permutation} here~~
}
  \item{SNP.CV}{
%%     ~~Describe \code{SNP.CV} here~~
}
  \item{genoFormat}{
%%     ~~Describe \code{genoFormat} here~~
}
  \item{hasGenotype}{
%%     ~~Describe \code{hasGenotype} here~~
}
  \item{byFile}{
%%     ~~Describe \code{byFile} here~~
}
  \item{fullGD}{
%%     ~~Describe \code{fullGD} here~~
}
  \item{PC}{
%%     ~~Describe \code{PC} here~~
}
  \item{GI}{
%%     ~~Describe \code{GI} here~~
}
  \item{Timmer}{
%%     ~~Describe \code{Timmer} here~~
}
  \item{Memory}{
%%     ~~Describe \code{Memory} here~~
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
function (Y, G = NULL, GD = NULL, GM = NULL, KI = NULL, Z = NULL, 
    CV = NULL, CV.Inheritance = NULL, SNP.P3D = TRUE, GP = NULL, 
    GK = NULL, group.from = 1e+06, group.to = 1, group.by = 10, 
    kinship.cluster = "average", kinship.group = "Mean", kinship.algorithm = NULL, 
    DPP = 50000, ngrid = 100, llin = -10, ulim = 10, esp = 1e-10, 
    file.path = NULL, file.from = NULL, file.to = NULL, file.total = NULL, 
    file.fragment = 512, file.G = NULL, file.Ext.G = NULL, file.GD = NULL, 
    file.GM = NULL, file.Ext.GD = NULL, file.Ext.GM = NULL, SNP.MAF = 0, 
    FDR.Rate = 1, SNP.FDR = 1, SNP.effect = "Add", SNP.impute = "Middle", 
    PCA.total = 0, GAPIT.Version = GAPIT.Version, name.of.trait, 
    GT = NULL, SNP.fraction = 1, seed = 123, BINS = 20, SNP.test = TRUE, 
    SNP.robust = "FaST", LD.chromosome = NULL, LD.location = NULL, 
    LD.range = NULL, bin.from = 10000, bin.to = 5e+06, bin.by = 1000, 
    inclosure.from = 10, inclosure.to = 1000, inclosure.by = 10, 
    SNP.permutation = FALSE, SNP.CV = NULL, genoFormat = NULL, 
    hasGenotype = NULL, byFile = NULL, fullGD = NULL, PC = NULL, 
    GI = NULL, Timmer = NULL, Memory = NULL, sangwich.top = NULL, 
    sangwich.bottom = NULL, QC = TRUE, GTindex = NULL, LD = 0.05, 
    file.output = TRUE, cutOff = 0.01, Model.selection = FALSE, 
    Create.indicator = FALSE, QTN = NULL, QTN.round = 1, QTN.limit = 0, 
    QTN.update = TRUE, QTN.method = "Penalty", Major.allele.zero = FALSE) 
{
    shortcut = FALSE
    LL.save = 1e+10
    thisY = Y[, 2]
    thisY = thisY[!is.na(thisY)]
    if (length(thisY) < 3) {
        shortcut = TRUE
    }
    else {
        if (var(thisY) == 0) 
            shortcut = TRUE
    }
    if (shortcut) {
        print(paste("Y is empty. No GWAS/GS performed for ", 
            name.of.trait, sep = ""))
        return(list(compression = NULL, kinship.optimum = NULL, 
            kinship = KI, PC = PC, GWAS = NULL, GPS = NULL, Pred = NULL, 
            REMLs = NULL, Timmer = Timmer, Memory = Memory))
    }
    print("------------Examining data (QC)------------------------------------------")
    if (is.null(Y)) 
        stop("GAPIT says: Phenotypes must exist.")
    if (is.null(KI) & missing(GD) & kinship.algorithm != "SUPER") 
        stop("GAPIT says: Kinship is required. As genotype is not provided, kinship can not be created.")
    if (is.null(GD) & is.null(GT)) {
        GT = as.matrix(Y[, 1])
        GD = matrix(1, nrow(Y), 1)
        GI = as.data.frame(matrix(0, 1, 3))
        colnames(GI) = c("SNP", "Chromosome", "Position")
    }
    if (PCA.total > 0 & !is.null(CV)) 
        CV = GAPIT.CVMergePC(CV, PC)
    if (PCA.total > 0 & is.null(CV)) 
        CV = PC
    if (kinship.algorithm != "None" & kinship.algorithm != "SUPER" & 
        is.null(Z)) {
        taxa = as.character(Y[, 1])
        Z = as.data.frame(diag(1, nrow(Y)))
        Z = rbind(taxa, Z)
        taxa = c("Taxa", as.character(taxa))
        Z = cbind(taxa, Z)
    }
    if (kinship.algorithm != "None" & kinship.algorithm != "SUPER" & 
        !is.null(Z)) {
        if (nrow(Z) - 1 < nrow(Y)) 
            Z = GAPIT.ZmatrixFormation(Z = Z, Y = Y)
    }
    noCV = FALSE
    if (is.null(CV)) {
        noCV = TRUE
        CV = Y[, 1:2]
        CV[, 2] = 1
        colnames(CV) = c("taxa", "overall")
    }
    print("QC is in process...")
    CVI <- CV
    if (QC) {
        qc <- GAPIT.QC(Y = Y, KI = KI, GT = GT, CV = CV, Z = Z, 
            GK = GK)
        GTindex = qc$GTindex
        Y = qc$Y
        KI = qc$KI
        CV = qc$CV
        Z = qc$Z
        GK = qc$GK
        if (noCV) 
            CVI = qc$CV
    }
    if (kinship.algorithm == "None") {
        if (min(CV[, 2]) == max(CV[, 2])) 
            CV = NULL
        theTDP = GAPIT.TDP(Y = Y, CV = CV, SNP = as.data.frame(cbind(GT[GTindex], 
            as.matrix(as.data.frame(GD[GTindex, ])))), QTN = QTN, 
            Round = QTN.round, QTN.limit = QTN.limit, QTN.update = QTN.update, 
            Method = QTN.method)
        theGWAS = cbind(GM, theTDP$p, NA, NA, NA)
        return(list(Compression = NULL, kinship.optimum = NULL, 
            kinship = NULL, PC = NULL, GWAS = theGWAS, GPS = NULL, 
            Pred = NULL, REMLs = NULL, QTN = theTDP$QTN, Timmer = Timmer, 
            Memory = Memory))
    }
    print("The value of QC is")
    print(QC)
    rm(qc)
    gc()
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "QC")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "QC")
    byPass.top = FALSE
    byPass = FALSE
    if (!is.null(sangwich.top)) 
        byPass.top = ((sangwich.top == "FaST" | sangwich.top == 
            "SUPER" | sangwich.top == "DC"))
    if (!is.null(sangwich.bottom)) 
        byPass = ((sangwich.bottom == "FaST" | sangwich.bottom == 
            "SUPER" | sangwich.bottom == "DC"))
    print("Try to group from and to were set to 1")
    if (byPass) {
        print("group from and to were set to 1")
        group.from = 1
        group.to = 1
    }
    print("------------Examining data (QC) done-------------------------------------")
    if (!is.null(sangwich.top) & is.null(GP)) {
        print("-------------------Sandwich top bun-----------------------------------")
        if (is.null(GK)) {
            set.seed(1)
            nY = floor(nrow(Y) * 0.9)
            nG = ncol(GD)
            if (nG > nY) {
                snpsam = sample(1:nG, nY)
            }
            else {
                snpsam = 1:nG
            }
            GK = GD[GTindex, snpsam]
            SNPVar = apply(as.matrix(GK), 2, var)
            GK = GK[, SNPVar > 0]
            GK = cbind(as.data.frame(GT[GTindex]), as.data.frame(GK))
        }
        GP = GAPIT.Bread(Y = Y, CV = CV, Z = Z, KI = KI, GK = GK, 
            GD = cbind(as.data.frame(GT), as.data.frame(GD)), 
            GM = GI, method = sangwich.top, GTindex = GTindex, 
            LD = LD)$GWAS
        GK = NULL
        print("-------------------Sagnwich top bun: done-----------------------------")
    }
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "SagnwichTop")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "SagnwichTop")
    print("-------------------Sandwich burger and dressing------------------------")
    if (group.from > group.to) 
        stop("GAPIT says: group.to should  be larger than group.from. Please correct them!")
    if (is.null(CV) | (!is.null(CV) & group.to < ncol(CV))) {
        group.from = 1
        group.to = 1
        warning("The upper bound of groups (group.to) is not sufficient. both boundries were set to a and GLM is performed!")
    }
    if (!is.null(CV) & group.from < 1) {
        group.from = 1
        warning("The lower bound of groups should be 1 at least. It was set to 1!")
    }
    nk = 1e+09
    if (!is.null(KI)) 
        nk = min(nk, nrow(KI))
    if (!is.null(GK)) 
        nk = min(nk, nrow(GK))
    if (!is.null(KI)) {
        if (group.to > nk) {
            group.to = nk
            warning("The upper bound of groups is too high. It was set to the size of kinship!")
        }
        if (group.from > nk) {
            group.from = nk
            warning("The lower bound of groups is too high. It was set to the size of kinship!")
        }
    }
    if (!is.null(CV)) {
        if (group.to <= ncol(CV) + 1) {
            group.from = ncol(CV) + 2
            group.to = ncol(CV) + 2
            warning("The upper bound of groups (group.to) is not sufficient. both boundries were set to their minimum and GLM is performed!")
        }
    }
    if (inclosure.to > nrow(Y)) 
        inclosure.to = nrow(Y) - 1
    bin.level = seq(bin.from, bin.to, by = bin.by)
    inclosure = seq(inclosure.from, inclosure.to, by = inclosure.by)
    GROUP = seq(group.to, group.from, by = -group.by)
    if (missing("kinship.cluster")) 
        kinship.cluster = c("ward", "single", "complete", "average", 
            "mcquitty", "median", "centroid")
    if (missing("kinship.group")) 
        kinship.group = c("Mean", "Max", "Min", "Median")
    numSetting = length(GROUP) * length(kinship.cluster) * length(kinship.group) * 
        length(bin.level) * length(inclosure)
    ys = as.matrix(Y[2])
    X0 = as.matrix(CV[, -1])
    CV.taxa = CVI[, 1]
    count = 0
    Compression = matrix(, numSetting, 6)
    colnames(Compression) = c("Type", "Cluster", "Group", "REML", 
        "VA", "VE")
    if (min(X0[, 1]) != max(X0[, 1])) 
        X0 <- cbind(1, X0)
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "DataProcessing")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "DataProcessing")
    print("-------------------------Iteration in process--------------------------")
    print(paste("Total iterations: ", numSetting, sep = ""))
    for (bin in bin.level) {
        for (inc in inclosure) {
            if (!byPass & (!is.null(GK) | !is.null(GP))) {
                print("Grilling KI...")
                myGenotype <- GAPIT.Genotype(G = NULL, GD = cbind(as.data.frame(GT), 
                  as.data.frame(GD)), GM = GI, KI = NULL, kinship.algorithm = kinship.algorithm, 
                  PCA.total = 0, SNP.fraction = SNP.fraction, 
                  SNP.test = SNP.test, file.path = file.path, 
                  file.from = file.from, file.to = file.to, file.total = file.total, 
                  file.fragment = file.fragment, file.G = file.G, 
                  file.Ext.G = file.Ext.G, file.GD = file.GD, 
                  file.GM = file.GM, file.Ext.GD = file.Ext.GD, 
                  file.Ext.GM = file.Ext.GM, SNP.MAF = SNP.MAF, 
                  FDR.Rate = FDR.Rate, SNP.FDR = SNP.FDR, SNP.effect = SNP.effect, 
                  SNP.impute = SNP.impute, LD.chromosome = LD.chromosome, 
                  LD.location = LD.location, LD.range = LD.range, 
                  GP = GP, GK = GK, bin.size = bin, inclosure.size = inc, 
                  SNP.CV = SNP.CV, Timmer = Timmer, Memory = Memory, 
                  GTindex = GTindex, sangwich.top = NULL, sangwich.bottom = sangwich.bottom, 
                  file.output = file.output, Create.indicator = Create.indicator, 
                  Major.allele.zero = Major.allele.zero)
                Timmer = myGenotype$Timmer
                Memory = myGenotype$Memory
                KI = myGenotype$KI
                nk = nrow(KI)
                GROUP = GROUP[GROUP <= nk]
            }
            for (ca in kinship.cluster) {
                for (group in GROUP) {
                  for (kt in kinship.group) {
                    if (numSetting == 1 & hasGenotype) {
                      optOnly = FALSE
                    }
                    else {
                      optOnly = TRUE
                    }
                    if (!SNP.test) 
                      optOnly = TRUE
                    if (optOnly | Model.selection) {
                      colInclude = 1
                      optOnly = TRUE
                    }
                    else {
                      colInclude = c(1:ncol(GD))
                    }
                    if (!optOnly) 
                      print("Compressing and Genome screening...")
                    count = count + 1
                    if (!byPass) {
                      if (count == 1) 
                        print("-------Mixed model with Kinship-----------------------------")
                      if (group < ncol(X0) + 1) 
                        group = 1
                      cp <- GAPIT.Compress(KI = KI, kinship.cluster = ca, 
                        kinship.group = kt, GN = group, Timmer = Timmer, 
                        Memory = Memory)
                      Timmer = cp$Timmer
                      Memory = cp$Memory
                      Timmer = GAPIT.Timmer(Timmer = Timmer, 
                        Infor = "PreP3D 2_cp")
                      Memory = GAPIT.Memory(Memory = Memory, 
                        Infor = "PreP3D 2_cp")
                      bk <- GAPIT.Block(Z = Z, GA = cp$GA, KG = cp$KG)
                      Timmer = GAPIT.Timmer(Timmer = Timmer, 
                        Infor = "PreP3D 2_bk")
                      Memory = GAPIT.Memory(Memory = Memory, 
                        Infor = "PreP3D 2 bk")
                      zc <- GAPIT.ZmatrixCompress(Z = Z, GAU = bk$GA)
                      Timmer = GAPIT.Timmer(Timmer = Timmer, 
                        Infor = "PreP3D 2_zc")
                      Memory = GAPIT.Memory(Memory = Memory, 
                        Infor = "PreP3D 2 zc")
                      zrow = nrow(zc$Z)
                      zcol = ncol(zc$Z) - 1
                      Timmer = GAPIT.Timmer(Timmer = Timmer, 
                        Infor = "Prio PreP3D")
                      Memory = GAPIT.Memory(Memory = Memory, 
                        Infor = "Prio PreP3D")
                      p3d <- GAPIT.EMMAxP3D(ys = ys, xs = as.matrix(as.data.frame(GD[GTindex, 
                        colInclude])), K = as.matrix(bk$KW), 
                        Z = matrix(as.numeric(as.matrix(zc$Z[, 
                          -1])), nrow = zrow, ncol = zcol), X0 = X0, 
                        CVI = CVI, CV.Inheritance = CV.Inheritance, 
                        GI = GI, SNP.P3D = SNP.P3D, Timmer = Timmer, 
                        Memory = Memory, fullGD = fullGD, SNP.permutation = SNP.permutation, 
                        GP = GP, file.path = file.path, file.from = file.from, 
                        file.to = file.to, file.total = file.total, 
                        file.fragment = file.fragment, byFile = byFile, 
                        file.G = file.G, file.Ext.G = file.Ext.G, 
                        file.GD = file.GD, file.GM = file.GM, 
                        file.Ext.GD = file.Ext.GD, file.Ext.GM = file.Ext.GM, 
                        GTindex = GTindex, genoFormat = genoFormat, 
                        optOnly = optOnly, SNP.effect = SNP.effect, 
                        SNP.impute = SNP.impute, name.of.trait = name.of.trait, 
                        Create.indicator = Create.indicator, 
                        Major.allele.zero = Major.allele.zero)
                      Timmer = p3d$Timmer
                      Memory = p3d$Memory
                      Timmer = GAPIT.Timmer(Timmer = Timmer, 
                        Infor = "Post PreP3D")
                      Memory = GAPIT.Memory(Memory = Memory, 
                        Infor = "Post PreP3D")
                      print(paste(count, "of", numSetting, "--", 
                        "Vg=", round(p3d$vgs, 4), "VE=", round(p3d$ves, 
                          4), "-2LL=", round(p3d$REMLs, 2), "  Clustering=", 
                        ca, "  Group number=", group, "  Group kinship=", 
                        kt, sep = " "))
                      if (count == 1) {
                        KI.save = KI
                        LL.save = p3d$REMLs
                      }
                      else {
                        if (p3d$REMLs < LL.save) {
                          KI.save = KI
                          LL.save = p3d$REMLs
                        }
                      }
                      Compression[count, 1] = kt
                      Compression[count, 2] = ca
                      Compression[count, 3] = group
                      Compression[count, 4] = p3d$REMLs
                      Compression[count, 5] = p3d$vgs
                      Compression[count, 6] = p3d$ves
                    }
                    else {
                      if (count == 1) 
                        print("-------The burger is SNP-----------------------------------")
                      if (!is.null(GP)) {
                        myGenotype <- GAPIT.Genotype(G = NULL, 
                          GD = NULL, GM = GI, KI = NULL, kinship.algorithm = "SUPER", 
                          PCA.total = 0, SNP.fraction = SNP.fraction, 
                          SNP.test = SNP.test, file.path = file.path, 
                          file.from = file.from, file.to = file.to, 
                          file.total = file.total, file.fragment = file.fragment, 
                          file.G = file.G, file.Ext.G = file.Ext.G, 
                          file.GD = file.GD, file.GM = file.GM, 
                          file.Ext.GD = file.Ext.GD, file.Ext.GM = file.Ext.GM, 
                          SNP.MAF = SNP.MAF, FDR.Rate = FDR.Rate, 
                          SNP.FDR = SNP.FDR, SNP.effect = SNP.effect, 
                          SNP.impute = SNP.impute, LD.chromosome = LD.chromosome, 
                          LD.location = LD.location, LD.range = LD.range, 
                          GP = GP, GK = NULL, bin.size = bin, 
                          inclosure.size = inc, SNP.CV = SNP.CV, 
                          GTindex = GTindex, sangwich.top = NULL, 
                          sangwich.bottom = sangwich.bottom, 
                          Timmer = Timmer, Memory = Memory, file.output = file.output, 
                          Create.indicator = Create.indicator, 
                          Major.allele.zero = Major.allele.zero)
                        Timmer = myGenotype$Timmer
                        Memory = myGenotype$Memory
                        Timmer = GAPIT.Timmer(Timmer = Timmer, 
                          Infor = "Genotype for burger")
                        Memory = GAPIT.Memory(Memory = Memory, 
                          Infor = "Genotype for burger")
                        GK = GD[GTindex, myGenotype$SNP.QTN]
                        SNPVar = apply(as.matrix(GK), 2, var)
                        GK = GK[, SNPVar > 0]
                        GK = cbind(as.data.frame(GT[GTindex]), 
                          as.data.frame(GK))
                      }
                      if (!is.null(GK) & numSetting > 1) {
                        print("-------Calculating likelihood-----------------------------------")
                        myBurger = GAPIT.Burger(Y = Y, CV = CV, 
                          GK = GK)
                        myREML = myBurger$REMLs
                        myVG = myBurger$vg
                        myVE = myBurger$ve
                      }
                      else {
                        myREML = NA
                        myVG = NA
                        myVE = NA
                      }
                      if (count == 1) {
                        GK.save = GK
                        LL.save = myREML
                      }
                      else {
                        if (myREML < LL.save) {
                          GK.save = GK
                          LL.save = myREML
                        }
                      }
                      Compression[count, 1] = 1
                      Compression[count, 2] = bin
                      Compression[count, 3] = inc
                      Compression[count, 4] = myREML
                      Compression[count, 5] = myVG
                      Compression[count, 6] = myVG
                      print(Compression[count, ])
                    }
                  }
                  if (group == 1) 
                    break
                }
            }
        }
    }
    if (Model.selection == TRUE) {
        print("------------------------Model selection for optimal number of PCs and Covariates-------------------------------------------------")
        KI = KI.save
        if (numSetting > 1) {
            Compression = Compression[order(as.numeric(Compression[, 
                4]), decreasing = FALSE), ]
            kt = Compression[1, 1]
            ca = Compression[1, 2]
            group = Compression[1, 3]
        }
        cp <- GAPIT.Compress(KI = KI, kinship.cluster = ca, kinship.group = kt, 
            GN = group, Timmer = Timmer, Memory = Memory)
        Timmer = cp$Timmer
        Memory = cp$Memory
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_cp")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2_cp")
        bk <- GAPIT.Block(Z = Z, GA = cp$GA, KG = cp$KG)
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_bk")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 bk")
        zc <- GAPIT.ZmatrixCompress(Z = Z, GAU = bk$GA)
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_zc")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 zc")
        z0 = as.matrix(zc$Z[, -1])
        Z1 = matrix(as.numeric(z0), nrow = nrow(z0), ncol = ncol(z0))
        BIC <- rep(NA, ncol(X0))
        LogLike <- rep(NA, ncol(X0))
        for (i in 1:ncol(X0)) {
            X0.test <- as.matrix(X0[, 1:i])
            p3d <- GAPIT.EMMAxP3D(ys = ys, xs = as.matrix(as.data.frame(GD[, 
                1])), K = as.matrix(bk$KW), Z = Z1, X0 = X0.test, 
                CVI = CVI, CV.Inheritance = CV.Inheritance, GI = GI, 
                SNP.P3D = SNP.P3D, Timmer = Timmer, Memory = Memory, 
                fullGD = fullGD, SNP.permutation = SNP.permutation, 
                GP = GP, file.path = file.path, file.from = file.from, 
                file.to = file.to, file.total = file.total, file.fragment = file.fragment, 
                byFile = byFile, file.G = file.G, file.Ext.G = file.Ext.G, 
                file.GD = file.GD, file.GM = file.GM, file.Ext.GD = file.Ext.GD, 
                file.Ext.GM = file.Ext.GM, GTindex = GTindex, 
                genoFormat = genoFormat, optOnly = TRUE, SNP.effect = SNP.effect, 
                SNP.impute = SNP.impute, name.of.trait = name.of.trait, 
                Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)
            k.num.param <- 2 + i
            LogLike[i] <- p3d$logLM
            BIC[i] <- p3d$logLM - (0.5 * k.num.param * log(nrow(Z1)))
        }
        Optimum.from.BIC <- which(BIC == max(BIC))
        print(paste("-----------------------The optimal number of PCs/covariates is ", 
            (Optimum.from.BIC - 1), " -------------------------", 
            sep = ""))
        BIC.Vector <- cbind(as.matrix(rep(0:(ncol(X0) - 1))), 
            as.matrix(BIC), as.matrix(LogLike))
        colnames(BIC.Vector) <- c("Number of PCs/Covariates", 
            "BIC (larger is better) - Schwarz 1978", "log Likelihood Function Value")
        write.table(BIC.Vector, paste("GAPIT.", name.of.trait, 
            ".BIC.Model.Selection.Results.csv", sep = ""), quote = FALSE, 
            sep = ",", row.names = FALSE, col.names = TRUE)
        X0 <- X0[, 1:(Optimum.from.BIC)]
        if (Optimum.from.BIC == 1) {
            X0 <- as.matrix(X0)
        }
        print("The dimension of X0 after model selection is:")
        print(dim(X0))
        print("The head of X0 after model selection is")
        print(head(X0))
    }
    print("---------------------Sandwich bottom bun-------------------------------")
    print("Compression")
    print(Compression)
    if (numSetting == 1) {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GWAS")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GWAS")
    }
    if ((numSetting > 1) | (!is.null(sangwich.bottom) & !byPass) | 
        Model.selection) {
        print("Genomic screening...")
        optOnly = FALSE
        if (!hasGenotype) 
            optOnly = TRUE
        if (!SNP.test) 
            optOnly = TRUE
        if (optOnly) {
            colInclude = 1
        }
        else {
            colInclude = c(1:ncol(GD))
        }
        if (numSetting > 1) {
            Compression = Compression[order(as.numeric(Compression[, 
                4]), decreasing = FALSE), ]
            kt = Compression[1, 1]
            ca = Compression[1, 2]
            group = Compression[1, 3]
            print(paste("Optimum: ", Compression[1, 2], Compression[1, 
                1], Compression[1, 3], Compression[1, 5], Compression[1, 
                6], Compression[1, 4], sep = " "))
        }
        print("--------------  Sandwich bottom ------------------------")
        if (!byPass) {
            print("--------------  Sandwich bottom with raw burger------------------------")
            if (Model.selection == FALSE) {
                if (is.null(sangwich.bottom)) 
                  KI = KI.save
                cp <- GAPIT.Compress(KI = KI, kinship.cluster = ca, 
                  kinship.group = kt, GN = group, Timmer = Timmer, 
                  Memory = Memory)
                Timmer = cp$Timmer
                Memory = cp$Memory
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_cp")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2_cp")
                bk <- GAPIT.Block(Z = Z, GA = cp$GA, KG = cp$KG)
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_bk")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 bk")
                zc <- GAPIT.ZmatrixCompress(Z = Z, GAU = bk$GA)
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PreP3D 2_zc")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "PreP3D 2 zc")
                z0 = as.matrix(zc$Z[, -1])
                Z1 = matrix(as.numeric(z0), nrow = nrow(z0), 
                  ncol = ncol(z0))
            }
            print("--------------EMMAxP3D with the optimum setting-----------------------")
            p3d <- GAPIT.EMMAxP3D(ys = ys, xs = as.matrix(as.data.frame(GD[GTindex, 
                colInclude])), K = as.matrix(bk$KW), Z = Z1, 
                X0 = as.matrix(X0), CVI = CVI, CV.Inheritance = CV.Inheritance, 
                GI = GI, SNP.P3D = SNP.P3D, Timmer = Timmer, 
                Memory = Memory, fullGD = fullGD, SNP.permutation = SNP.permutation, 
                GP = GP, file.path = file.path, file.from = file.from, 
                file.to = file.to, file.total = file.total, file.fragment = file.fragment, 
                byFile = byFile, file.G = file.G, file.Ext.G = file.Ext.G, 
                file.GD = file.GD, file.GM = file.GM, file.Ext.GD = file.Ext.GD, 
                file.Ext.GM = file.Ext.GM, GTindex = GTindex, 
                genoFormat = genoFormat, optOnly = optOnly, SNP.effect = SNP.effect, 
                SNP.impute = SNP.impute, name.of.trait = name.of.trait, 
                Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GWAS")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "GWAS")
            print("--------------EMMAxP3D with the optimum setting done------------------")
        }
    }
    if (byPass) {
        print("---------------Sandwich bottom with grilled burger---------------------")
        print("---------------Sandwich bottom: reload bins ---------------------------")
        GK = GK.save
        myBread = GAPIT.Bread(Y = Y, CV = CV, Z = Z, GK = GK, 
            GD = cbind(as.data.frame(GT), as.data.frame(GD)), 
            GM = GI, method = sangwich.bottom, GTindex = GTindex, 
            LD = LD)
        print("SUPER saving results...")
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GWAS")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GWAS")
    }
    print("--------------------Final results presentations------------------------")
    if (!byPass) {
        if (length(bk$KW) > 1 & length(bk$KW) < length(KI) & 
            length(bk$KW) < 1000 & file.output) {
            pdf(paste("GAPIT.", name.of.trait, ".Kin.Optimum.pdf", 
                sep = ""), width = 12, height = 12)
            par(mar = c(25, 25, 25, 25))
            heatmap.2(as.matrix(bk$KW), cexRow = 0.2, cexCol = 0.2, 
                col = rev(heat.colors(256)), scale = "none", 
                symkey = FALSE, trace = "none")
            dev.off()
        }
    }
    if (byFile & !fullGD) {
        print("Loading GWAS results from file...")
        for (file in file.from:file.to) {
            frag = 1
            numSNP = file.fragment
            while (numSNP == file.fragment) {
                print(paste("Current file ", file, "Fragment: ", 
                  frag))
                theGI <- try(read.table(paste("GAPIT.TMP.GI.", 
                  name.of.trait, file, ".", frag, ".txt", sep = ""), 
                  head = TRUE), silent = TRUE)
                theP <- try(read.table(paste("GAPIT.TMP.ps.", 
                  name.of.trait, file, ".", frag, ".txt", sep = ""), 
                  head = FALSE), silent = TRUE)
                theMAF <- try(read.table(paste("GAPIT.TMP.maf.", 
                  name.of.trait, file, ".", frag, ".txt", sep = ""), 
                  head = FALSE), silent = TRUE)
                thenobs <- try(read.table(paste("GAPIT.TMP.nobs.", 
                  name.of.trait, file, ".", frag, ".txt", sep = ""), 
                  head = FALSE), silent = TRUE)
                thersquare_base <- try(read.table(paste("GAPIT.TMP.rsquare.base.", 
                  name.of.trait, file, ".", frag, ".txt", sep = ""), 
                  head = FALSE), silent = TRUE)
                thersquare <- try(read.table(paste("GAPIT.TMP.rsquare.", 
                  name.of.trait, file, ".", frag, ".txt", sep = ""), 
                  head = FALSE), silent = TRUE)
                theeffect.est <- try(read.table(paste("GAPIT.TMP.effect.est.", 
                  name.of.trait, file, ".", frag, ".txt", sep = ""), 
                  head = FALSE), silent = TRUE)
                if (inherits(theGI, "try-error")) {
                  numSNP = 0
                }
                else {
                  numSNP = nrow(theGI)
                  colnames(theP) = "P"
                  colnames(theMAF) = "MAF"
                  colnames(thenobs) = "nobs"
                  colnames(thersquare_base) = "Base.Model.R.square"
                  colnames(thersquare) = "Model.R.square"
                  colnames(theeffect.est) = "Effect.Est"
                  colnames(theGI) = colnames(GI)
                  if (file == file.from & frag == 1) {
                    GI = theGI
                    allP = theP
                    allMAF = theMAF
                    allnobs = thenobs
                    allrsquare_base = thersquare_base
                    allrsquare = thersquare
                    alleffect.est = theeffect.est
                  }
                  else {
                    allP = as.data.frame(rbind(as.matrix(allP), 
                      as.matrix(theP)))
                    allMAF = as.data.frame(rbind(as.matrix(allMAF), 
                      as.matrix(theMAF)))
                    allnobs = as.data.frame(rbind(as.matrix(allnobs), 
                      as.matrix(thenobs)))
                    allrsquare_base = as.data.frame(rbind(as.matrix(allrsquare_base), 
                      as.matrix(thersquare_base)))
                    allrsquare = as.data.frame(rbind(as.matrix(allrsquare), 
                      as.matrix(thersquare)))
                    alleffect.est = as.data.frame(rbind(as.matrix(alleffect.est), 
                      as.matrix(theeffect.est)))
                    GI = as.data.frame(rbind(as.matrix(GI), as.matrix(theGI)))
                  }
                }
                frag = frag + 1
            }
        }
        p3d$ps = allP
        p3d$maf = allMAF
        p3d$nobs = allnobs
        p3d$rsquare_base = allrsquare_base
        p3d$rsquare = allrsquare
        p3d$effect.est = alleffect.est
        theFile = paste("GAPIT.TMP.", name.of.trait, ".*")
        system("cmd /c del \"GAPIT.TMP*.*\"")
        system("cmd /c del \"GAPIT.TMP*.*\"")
        print("GWAS results loaded from all files succesfully!")
    }
    print("Generating summary")
    GWAS = NULL
    GPS = NULL
    rm(zc)
    gc()
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Final")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Final")
    print("Genomic Breeding Values (GBV) ...")
    if (!byPass) {
        if (length(bk$KW) > ncol(X0)) {
            gs <- GAPIT.GS(KW = bk$KW, KO = bk$KO, KWO = bk$KWO, 
                GAU = bk$GAU, UW = cbind(p3d$BLUP, p3d$PEV))
        }
        print("Writing GBV and Acc...")
        GPS = NULL
        if (length(bk$KW) > ncol(X0)) 
            GPS = gs$BLUP
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GPS")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GPS")
        print("GBV and accuracy distribution...")
        if (length(bk$KW) > ncol(X0) & file.output) {
            GAPIT.GS.Visualization(gsBLUP = gs$BLUP, BINS = BINS, 
                name.of.trait = name.of.trait)
        }
        print("Compression portfolios...")
        if (file.output) 
            GAPIT.Compression.Visualization(Compression = Compression, 
                name.of.trait = name.of.trait)
        print("Compression Visualization done")
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Compression.Visualization")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Compression.Visualization")
        ps = p3d$ps
        nobs = p3d$nobs
        maf = p3d$maf
        rsquare_base = p3d$rsquare_base
        rsquare = p3d$rsquare
        effect.est = p3d$effect.est
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Extract p3d results")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Extract p3d results")
    }
    else {
        print("The head of myBread$GWAS is")
        print(head(myBread$GWAS))
        GPS = myBread$BLUP
        ps = myBread$GWAS[, 4]
        nobs = myBread$GWAS[, 6]
        maf = myBread$GWAS[, 5] * 0 + 0.5
        rsquare_base = rep(NA, length(ps))
        rsquare = rep(NA, length(ps))
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Extract bread results")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Extract bread results")
    }
    if ((!byPass) & (!Model.selection)) {
        BLUE = data.frame(cbind(data.frame(CV.taxa), data.frame(p3d$BLUE)))
        colnames(BLUE) = c("Taxa", "BLUE")
        BB = merge(gs$BLUP, BLUE, by.x = "Taxa", by.y = "Taxa")
        Prediction = BB[, 5] + BB[, 7]
        Pred = data.frame(cbind(BB, data.frame(Prediction)))
    }
    if (!byPass & file.output) {
        print("Exporting BLUP and Pred")
        try(write.table(gs$BLUP, paste("GAPIT.", name.of.trait, 
            ".BLUP.csv", sep = ""), quote = FALSE, sep = ",", 
            row.names = FALSE, col.names = TRUE))
        try(write.table(Pred, paste("GAPIT.", name.of.trait, 
            ".PRED.csv", sep = ""), quote = FALSE, sep = ",", 
            row.names = FALSE, col.names = TRUE))
    }
    if (byPass) {
        theK.return = NULL
    }
    else {
        theK.return = cp$KG
    }
    if (byPass) 
        Compression[1, 4] = 0
    if (hasGenotype & SNP.test & !is.na(Compression[1, 4])) {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Extract GWAS start")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Extract GWAS start")
        print("Filtering SNPs with MAF...")
        index = maf >= SNP.MAF
        PWI.Filtered = cbind(GI, ps, maf, nobs, rsquare_base, 
            rsquare)[index, ]
        colnames(PWI.Filtered) = c("SNP", "Chromosome", "Position ", 
            "P.value", "maf", "nobs", "Rsquare.of.Model.without.SNP", 
            "Rsquare.of.Model.with.SNP")
        if (!byPass) {
            if (Create.indicator) {
                GI.counter <- cbind(GI, seq(1:nrow(GI)))
                GI.counter.data.frame <- data.frame(GI.counter)
                colnames(GI.counter.data.frame) <- c("X1", "X2", 
                  "X3", "X4")
                effect.est.data.frame <- data.frame(effect.est)
                colnames(effect.est.data.frame) <- c("X1", "X2", 
                  "X3")
                GWAS.2 <- merge(GI.counter.data.frame, effect.est.data.frame, 
                  by.x = "X4", by.y = "X1")
                GWAS.2 <- GWAS.2[, -1]
                colnames(GWAS.2) <- c("SNP", "Chromosome", "Position ", 
                  "Genotype", "Allelic Effect Estimate")
            }
            if (!Create.indicator) {
                GWAS.2 <- cbind(GI, effect.est)
                colnames(GWAS.2) <- c("SNP", "Chromosome", "Position ", 
                  "Allelic Effect Estimate")
            }
        }
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "MAF filtered")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "MAF filtered")
        print("SNPs filtered with MAF")
        if (!is.null(PWI.Filtered)) {
            print("Calculating FDR...")
            PWIP <- GAPIT.Perform.BH.FDR.Multiple.Correction.Procedure(PWI = PWI.Filtered, 
                FDR.Rate = FDR.Rate, FDR.Procedure = "BH")
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Multiple Correction")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "Multiple Correction")
            print("QQ plot...")
            if (file.output) 
                GAPIT.QQ(P.values = PWIP$PWIP[, 4], name.of.trait = name.of.trait, 
                  DPP = DPP)
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "QQ plot")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "QQ plot")
            print("Manhattan plot (Genomewise)...")
            if (file.output) 
                GAPIT.Manhattan(GI.MP = PWIP$PWIP[, 2:4], name.of.trait = name.of.trait, 
                  DPP = DPP, plot.type = "Genomewise", cutOff = cutOff)
            print("Manhattan plot (Chromosomewise)...")
            if (file.output) 
                GAPIT.Manhattan(GI.MP = PWIP$PWIP[, 2:4], name.of.trait = name.of.trait, 
                  DPP = DPP, plot.type = "Chromosomewise", cutOff = cutOff)
            Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Manhattan plot")
            Memory = GAPIT.Memory(Memory = Memory, Infor = "Manhattan plot")
            print("Association table...")
            GWAS = PWIP$PWIP[PWIP$PWIP[, 9] <= SNP.FDR, ]
            if (file.output) {
                write.table(GWAS, paste("GAPIT.", name.of.trait, 
                  ".GWAS.Results.csv", sep = ""), quote = FALSE, 
                  sep = ",", row.names = FALSE, col.names = TRUE)
                if (!byPass) 
                  write.table(GWAS.2, paste("GAPIT.", name.of.trait, 
                    ".Allelic_Effect_Estimates.csv", sep = ""), 
                    quote = FALSE, sep = ",", row.names = FALSE, 
                    col.names = TRUE)
            }
        }
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Extract GWAS end")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Extract GWAS end")
    }
    if (file.output) 
        log = GAPIT.Log(Y = Y, KI = KI, Z = Z, CV = CV, SNP.P3D = SNP.P3D, 
            group.from = group.from, group.to = group.to, group.by = group.by, 
            kinship.cluster = kinship.cluster, kinship.group = kinship.group, 
            ngrid = ngrid, llin = llin, ulim = ulim, esp = esp, 
            name.of.trait = name.of.trait)
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Report")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Report")
    if (file.output) {
        file = paste("GAPIT.", name.of.trait, ".Timming.csv", 
            sep = "")
        write.table(Timmer, file, quote = FALSE, sep = ",", row.names = FALSE, 
            col.names = TRUE)
        file = paste("GAPIT.", name.of.trait, ".Memory.Stage.csv", 
            sep = "")
        write.table(Memory, file, quote = FALSE, sep = ",", row.names = FALSE, 
            col.names = TRUE)
    }
    print(paste(name.of.trait, "has been analyzed successfully!"))
    print(paste("The results are saved in the directory of ", 
        getwd()))
    print("==========================================================================================")
    if (byPass | Model.selection) 
        Pred <- NA
    return(list(Timmer = Timmer, Compression = Compression, kinship.optimum = theK.return, 
        kinship = KI, PC = PC, GWAS = GWAS, GPS = GPS, Pred = Pred, 
        REMLs = Compression[count, 4], Timmer = Timmer, Memory = Memory))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
