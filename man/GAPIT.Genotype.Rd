\name{GAPIT.Genotype}
\alias{GAPIT.Genotype}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
GAPIT.Genotype(G = NULL, GD = NULL, GM = NULL, KI = NULL, kinship.algorithm = NULL, SNP.effect = "Add", SNP.impute = "Middle", PCA.total = 0, seed = 123, SNP.fraction = 1, file.path = NULL, file.from = NULL, file.to = NULL, file.total = NULL, file.fragment = 1000, SNP.test = TRUE, file.G = NULL, file.Ext.G = NULL, file.GD = NULL, file.Ext.GD = NULL, file.GM = NULL, file.Ext.GM = NULL, SNP.MAF = 0.05, FDR.Rate = 0.05, SNP.FDR = 1, Timmer = NULL, Memory = NULL, LD.chromosome = NULL, LD.location = NULL, LD.range = NULL, SNP.CV = NULL, GP = NULL, GK = NULL, GTindex = NULL, bin.size = 1000, inclosure.size = 100, sangwich.top = NULL, sangwich.bottom = NULL, file.output = TRUE, Create.indicator = FALSE, Major.allele.zero = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
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
  \item{kinship.algorithm}{
%%     ~~Describe \code{kinship.algorithm} here~~
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
  \item{seed}{
%%     ~~Describe \code{seed} here~~
}
  \item{SNP.fraction}{
%%     ~~Describe \code{SNP.fraction} here~~
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
  \item{SNP.test}{
%%     ~~Describe \code{SNP.test} here~~
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
  \item{file.Ext.GD}{
%%     ~~Describe \code{file.Ext.GD} here~~
}
  \item{file.GM}{
%%     ~~Describe \code{file.GM} here~~
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
  \item{Timmer}{
%%     ~~Describe \code{Timmer} here~~
}
  \item{Memory}{
%%     ~~Describe \code{Memory} here~~
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
  \item{SNP.CV}{
%%     ~~Describe \code{SNP.CV} here~~
}
  \item{GP}{
%%     ~~Describe \code{GP} here~~
}
  \item{GK}{
%%     ~~Describe \code{GK} here~~
}
  \item{GTindex}{
%%     ~~Describe \code{GTindex} here~~
}
  \item{bin.size}{
%%     ~~Describe \code{bin.size} here~~
}
  \item{inclosure.size}{
%%     ~~Describe \code{inclosure.size} here~~
}
  \item{sangwich.top}{
%%     ~~Describe \code{sangwich.top} here~~
}
  \item{sangwich.bottom}{
%%     ~~Describe \code{sangwich.bottom} here~~
}
  \item{file.output}{
%%     ~~Describe \code{file.output} here~~
}
  \item{Create.indicator}{
%%     ~~Describe \code{Create.indicator} here~~
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
function (G = NULL, GD = NULL, GM = NULL, KI = NULL, kinship.algorithm = NULL, 
    SNP.effect = "Add", SNP.impute = "Middle", PCA.total = 0, 
    seed = 123, SNP.fraction = 1, file.path = NULL, file.from = NULL, 
    file.to = NULL, file.total = NULL, file.fragment = 1000, 
    SNP.test = TRUE, file.G = NULL, file.Ext.G = NULL, file.GD = NULL, 
    file.Ext.GD = NULL, file.GM = NULL, file.Ext.GM = NULL, SNP.MAF = 0.05, 
    FDR.Rate = 0.05, SNP.FDR = 1, Timmer = NULL, Memory = NULL, 
    LD.chromosome = NULL, LD.location = NULL, LD.range = NULL, 
    SNP.CV = NULL, GP = NULL, GK = NULL, GTindex = NULL, bin.size = 1000, 
    inclosure.size = 100, sangwich.top = NULL, sangwich.bottom = NULL, 
    file.output = TRUE, Create.indicator = FALSE, Major.allele.zero = FALSE) 
{
    print("Genotyping: numericalization, sampling kinship, PCs and much more...")
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Genotype start")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Genotype start")
    byData = !is.null(G) | !is.null(GD)
    byFile = !is.null(file.G) | !is.null(file.GD)
    hasGenotype = (byData | byFile)
    needKinPC = (is.null(KI) | PCA.total > 0 | kinship.algorithm == 
        "Separation")
    if (!is.null(KI) & !byData & !byFile & !SNP.test & kinship.algorithm != 
        "SUPER") {
        print("It return unexpected")
        return(list(GD = NULL, GI = NULL, GT = NULL, hasGenotype = FALSE, 
            genoFormat = NULL, KI = KI, PC = NULL, byFile = FALSE, 
            fullGD = TRUE, Timmer = Timmer, Memory = Memory))
    }
    fullGD = FALSE
    if (byData) 
        fullGD = TRUE
    if (byFile & SNP.fraction == 1 & needKinPC) 
        fullGD = TRUE
    if (!byData & !byFile & is.null(GK) & kinship.algorithm != 
        "SUPER") {
        if (is.null(KI) & is.null(GP) & is.null(GK)) 
            stop("GAPIT says: Kinship has to be provided or estimated from genotype!!!")
        return(list(GD = NULL, GI = NULL, GT = NULL, hasGenotype = FALSE, 
            genoFormat = NULL, KI = KI, PC = NULL, byFile = FALSE, 
            fullGD = TRUE, Timmer = Timmer, Memory = Memory))
    }
    genoFormat = "hapmap"
    if (is.null(G) & is.null(file.G)) 
        genoFormat = "EMMA"
    kinship.algorithm.save = kinship.algorithm
    kinship.algorithm = "SUPER"
    if (is.null(sangwich.top) & is.null(sangwich.bottom)) 
        kinship.algorithm = kinship.algorithm.save
    pass.top = FALSE
    if (!is.null(sangwich.top)) 
        pass.top = !(sangwich.top == "FaST" | sangwich.top == 
            "SUPER" | sangwich.top == "DC")
    pass.bottom = FALSE
    if (!is.null(sangwich.bottom)) 
        pass.bottom = !(sangwich.bottom == "FaST" | sangwich.bottom == 
            "SUPER" | sangwich.bottom == "DC")
    if (pass.top | pass.bottom) 
        kinship.algorithm = kinship.algorithm.save
    if (!is.null(file.from) & !is.null(file.to) & !is.null(file.total)) {
        if (file.total != (file.to - file.from + 1)) 
            stop("GAPIT says: Conflict among file (from, to and total)")
    }
    if (!is.null(file.from) & !is.null(file.to)) {
        if (file.to < file.from) 
            stop("GAPIT says: file.from should smaller than file.to")
    }
    if (is.null(file.from) & !is.null(file.to)) 
        stop("GAPIT says: file.from and file.to must be in pair)")
    if (!is.null(file.from) & is.null(file.to)) 
        stop("GAPIT says: file.from and file.to must be in pair)")
    if (!is.null(file.from) & !is.null(file.to)) 
        file.total = file.to - file.from + 1
    if (byFile & is.null(file.total)) 
        stop("GAPIT says: file.from and file.to must be provided!)")
    if (!is.null(GP) & !is.null(GK)) 
        stop("GAPIT Says: You can not provide GP and GK at same time")
    if (!is.null(GP) & !is.null(KI)) 
        stop("GAPIT Says: You can not provide GP and KI at same time")
    if (!is.null(GK) & !is.null(KI)) 
        stop("GAPIT says: You can not specify GK and KI at same time!!!")
    if (!is.null(GP) & !is.null(sangwich.top)) 
        stop("GAPIT Says: You provided GP. You can not spycify sangwich.top")
    if (!is.null(sangwich.top) & is.null(sangwich.bottom)) 
        stop("GAPIT Says: Top require its Bottom")
    if (is.null(sangwich.top) & !is.null(sangwich.bottom) & (is.null(GP) & 
        is.null(GK))) 
        stop("GAPIT Says: Uncovered Bottom (without TOP) requires GP or GK")
    if (is.null(sangwich.top) & is.null(sangwich.bottom) & (!is.null(GP) | 
        !is.null(GK))) 
        stop("GAPIT Says: You have provide GP or GK, you need to provide Bottom")
    if (!needKinPC & SNP.fraction < 1) 
        stop("GAPIT says: You did not require calculate kinship or PCs. SNP.fraction should not be specified!!!")
    if (!SNP.test & is.null(KI) & !byData & !byFile) 
        stop("GAPIT says: For SNP.test optioin, please input either use KI or use genotype")
    if (is.null(file.total) & !byData & byFile) 
        stop("GAPIT Ssays: Number of file should be provided: >=1")
    if (!is.null(G) & !is.null(GD)) 
        stop("GAPIT Ssays: Both hapmap and EMMA format exist, choose one only.")
    if (!is.null(file.GD) & is.null(file.GM) & (!is.null(GP) | 
        !is.null(GK))) 
        stop("GAPIT Ssays: Genotype data and map files should be in pair")
    if (is.null(file.GD) & !is.null(file.GM) & (!is.null(GP) | 
        !is.null(GK))) 
        stop("GAPIT Ssays: Genotype data and map files should be in pair")
    if (!is.null(GD) & is.null(GM) & (is.null(GP) & is.null(GK)) & 
        kinship.algorithm != "SUPER") 
        stop("GAPIT Says: Genotype data and map files should be in pair")
    if (is.null(GD) & !is.null(GM) & (is.null(GP) & is.null(GK)) & 
        kinship.algorithm != "SUPER") 
        stop("GAPIT Says: Genotype data and map files should be in pair")
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Genotype loaded")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Genotype loaded")
    GLD = NULL
    SNP.QTN = NULL
    GT = NULL
    if (!is.null(GM)) 
        GI = GM
    rm(GM)
    gc()
    if (!is.null(GD)) {
        GT = as.matrix(GD[, 1])
        GD = as.matrix(GD[, -1])
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "GT created from GD)")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "GT created from GD")
    }
    if (!is.null(G)) {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before HapMap")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Before HapMap")
        print(paste("Converting genotype...", sep = ""))
        hm = GAPIT.HapMap(G, SNP.effect = SNP.effect, SNP.impute = SNP.impute, 
            Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "after HapMap")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "after HapMap")
        if (!is.null(LD.chromosome)) {
            chromosome = (G[, 3] == LD.chromosome[1])
            bp = as.numeric(as.vector(G[, 4]))
            deviation = abs(bp - as.numeric(as.vector(LD.location[1])))
            location = deviation < as.numeric(as.vector(LD.range[1]))
            index = chromosome & location
            GLD = G[index, ]
        }
        else {
            GLD = NULL
        }
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "HapMap")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "HapMap")
        print(paste("Converting genotype done.", sep = ""))
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "G removed")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "G removed")
        GT = hm$GT
        GD = hm$GD
        GI = hm$GI
        rm(hm)
        gc()
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "hm removed")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "hm removed")
    }
    if (!byData & byFile) {
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "byFile")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "byFile")
        numFileUsed = file.to
        if (!needKinPC) 
            numFileUsed = file.from
        GD = NULL
        GT = NULL
        GI = NULL
        GLD = NULL
        for (file in file.from:numFileUsed) {
            frag = 1
            numSNP = file.fragment
            myFRG = NULL
            while (numSNP == file.fragment) {
                print(paste("Reading file: ", file, "Fragment: ", 
                  frag))
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before Fragment")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "Before Fragment")
                myFRG = GAPIT.Fragment(file.path = file.path, 
                  file.from = file.from, file.to = file.to, file.total = file.total, 
                  file.G = file.G, file.Ext.G = file.Ext.G, seed = seed, 
                  SNP.fraction = SNP.fraction, SNP.effect = SNP.effect, 
                  SNP.impute = SNP.impute, genoFormat = genoFormat, 
                  file.GD = file.GD, file.Ext.GD = file.Ext.GD, 
                  file.GM = file.GM, file.Ext.GM = file.Ext.GM, 
                  file.fragment = file.fragment, file = file, 
                  frag = frag, LD.chromosome = LD.chromosome, 
                  LD.location = LD.location, LD.range = LD.range, 
                  Create.indicator = Create.indicator, Major.allele.zero = Major.allele.zero)
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "After Fragment")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "After Fragment")
                if (is.null(GT) & !is.null(myFRG$GT)) 
                  GT = as.matrix(myFRG$GT)
                if (is.null(GD)) {
                  GD = myFRG$GD
                }
                else {
                  if (!is.null(myFRG$GD)) {
                    GD = cbind(GD, myFRG$GD)
                  }
                }
                if (is.null(GI)) {
                  GI = myFRG$GI
                }
                else {
                  if (!is.null(myFRG$GI)) {
                    colnames(myFRG$GI) = c("SNP", "Chromosome", 
                      "Position")
                    GI = as.data.frame(rbind(as.matrix(GI), as.matrix(myFRG$GI)))
                  }
                }
                if (is.null(G)) {
                  G = myFRG$G
                }
                else {
                  if (!is.null(myFRG$G)) {
                    G = as.data.frame(rbind(as.matrix(G), as.matrix(myFRG$G[-1, 
                      ])))
                  }
                }
                if (is.null(GLD)) {
                  GLD = myFRG$GLD
                }
                else {
                  if (!is.null(myFRG$GLD)) {
                    if (myFRG$heading) {
                      GLD = as.data.frame(rbind(as.matrix(GLD), 
                        as.matrix(myFRG$GLD[-1, ])))
                    }
                    else {
                      GLD = as.data.frame(rbind(as.matrix(GLD), 
                        as.matrix(myFRG$GLD)))
                    }
                  }
                }
                if (file == file.from & frag == 1) 
                  GT = as.matrix(myFRG$GT)
                frag = frag + 1
                if (!is.null(myFRG$GI)) {
                  numSNP = myFRG$linesRead[1]
                }
                else {
                  numSNP = 0
                }
                if (!needKinPC) 
                  numSNP = 0
                if (is.null(myFRG)) 
                  numSNP = 0
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "END this Fragment")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "END this Fragment")
            }
        }
        print("All files loaded")
    }
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Sampling genotype")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Sampling genotype")
    if (!is.null(KI) & file.output) {
        if (KI != 1) {
            if (nrow(KI) < 1000) {
                print("Plotting Kinship")
                theKin = as.matrix(KI[, -1])
                colnames(theKin) = KI[, 1]
                rownames(theKin) = KI[, 1]
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "set kinship")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "set kinship")
                print("Creating heat map for kinship...")
                pdf(paste("GAPIT.Kin.thirdPart.pdf", sep = ""), 
                  width = 12, height = 12)
                par(mar = c(25, 25, 25, 25))
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "prepare heatmap")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "prepare heatmap")
                heatmap.2(theKin, cexRow = 0.2, cexCol = 0.2, 
                  col = rev(heat.colors(256)), scale = "none", 
                  symkey = FALSE, trace = "none")
                dev.off()
                print("Kinship heat map PDF created!")
                Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "plot heatmap")
                Memory = GAPIT.Memory(Memory = Memory, Infor = "plot heatmap")
            }
        }
    }
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before SUPER")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Before SUPER")
    if (!is.null(GP) & kinship.algorithm == "SUPER" & !is.null(bin.size) & 
        !is.null(inclosure.size)) {
        mySpecify = GAPIT.Specify(GI = GI, GP = GP, bin.size = bin.size, 
            inclosure.size = inclosure.size)
        SNP.QTN = mySpecify$index
        if (!is.null(GD)) {
            GK = GD[, SNP.QTN]
            SNPVar = apply(as.matrix(GK), 2, var)
            GK = GK[, SNPVar > 0]
            GK = cbind(as.data.frame(GT), as.data.frame(GK))
        }
    }
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before PCA")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Before PCA")
    PC = NULL
    thePCA = NULL
    if (PCA.total > 0 | kinship.algorithm == "Separation") {
        thePCA = GAPIT.PCA(X = GD, taxa = GT, PC.number = PCA.total, 
            file.output = file.output)
        PC = thePCA$PCs[, 1:(1 + PCA.total)]
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "PCA")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "PCA")
        print("PC created")
    }
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Before creating kinship")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "Before creating kinship")
    if (is.null(KI) & (!is.null(GD) | !is.null(GK)) & kinship.algorithm != 
        "SUPER") {
        print("Calculating kinship...")
        if (!is.null(GK)) {
            thisGD = GK[, -1]
            myGT = as.matrix(GK[, 1])
            print("GK is used to create KI")
        }
        else {
            thisGD = GD
            myGT = GT
        }
        print(paste("Number of individuals and SNPs are ", nrow(thisGD), 
            " and ", ncol(thisGD)))
        theKin = NULL
        if (kinship.algorithm == "EMMA") {
            half.thisGD = as.matrix(0.5 * thisGD)
            if (length(which(is.na(half.thisGD))) > 0) {
                print("Substituting missing values with heterozygote for kinship matrrix calculation....")
                half.thisGD[which(is.na(half.thisGD))] = 1
            }
            theKin = emma.kinship(snps = t(as.matrix(0.5 * thisGD)), 
                method = "additive", use = "all")
        }
        if (kinship.algorithm == "Loiselle") 
            theKin = GAPIT.kinship.loiselle(snps = t(as.matrix(0.5 * 
                thisGD)), method = "additive", use = "all")
        if (kinship.algorithm == "VanRaden") 
            theKin = GAPIT.kinship.VanRaden(snps = as.matrix(thisGD))
        if (kinship.algorithm == "Separation") 
            theKin = GAPIT.kinship.separation(PCs = thePCA$PCs, 
                EV = thePCA$EV, nPCs = PCA.total)
        if (!is.null(theKin)) {
            colnames(theKin) = myGT
            rownames(theKin) = myGT
            print("kinship calculated")
            if (length(GT) < 1000 & file.output) {
                print("Creating heat map for kinship...")
                pdf(paste("GAPIT.Kin.", kinship.algorithm, ".pdf", 
                  sep = ""), width = 12, height = 12)
                par(mar = c(25, 25, 25, 25))
                heatmap.2(theKin, cexRow = 0.2, cexCol = 0.2, 
                  col = rev(heat.colors(256)), scale = "none", 
                  symkey = FALSE, trace = "none")
                dev.off()
                print("Kinship heat map created")
            }
            print("Adding IDs to kinship...")
            KI = cbind(myGT, as.data.frame(theKin))
            print("Writing kinship to file...")
            if (file.output) 
                write.table(KI, paste("GAPIT.Kin.", kinship.algorithm, 
                  ".csv", sep = ""), quote = FALSE, sep = ",", 
                  row.names = FALSE, col.names = FALSE)
            print("Kinship save as file")
            rm(theKin)
            gc()
        }
        Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "Estimating kinship")
        Memory = GAPIT.Memory(Memory = Memory, Infor = "Estimating kinship")
        print("Kinship created!")
    }
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "after creating kinship")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "after creating kinship")
    if (!is.null(GLD) & file.output) {
        if (nrow(GLD) > 500) {
            GLD = GLD[1, ]
            print("WARNING: The number of SNPs requested is beyond limitation. No LD plot created.")
        }
        if (nrow(GLD) > 1) {
            print("Plot LD...")
            hapmapgeno = data.frame(as.matrix(t(GLD[, -c(1:11)])))
            hapmapgeno[hapmapgeno == "NN"] = NA
            hapmapgeno[hapmapgeno == "XX"] = NA
            hapmapgeno[hapmapgeno == "--"] = NA
            hapmapgeno[hapmapgeno == "++"] = NA
            hapmapgeno[hapmapgeno == "//"] = NA
            LDdist = as.numeric(as.vector(GLD[, 4]))
            LDsnpName = GLD[, 1]
            colnames(hapmapgeno) = LDsnpName
            LDsnpName = LDsnpName[c(1, length(LDsnpName))]
            print("Getting genotype object")
            LDsnp = makeGenotypes(hapmapgeno, sep = "", method = as.genotype)
            print("Caling LDheatmap...")
            pdf(paste("GAPIT.LD.chromosom", LD.chromosome, "(", 
                round(max(0, LD.location - LD.range)/1e+06), 
                "_", round((LD.location + LD.range)/1e+06), "Mb)", 
                ".pdf", sep = ""), width = 12, height = 12)
            par(mar = c(25, 25, 25, 25))
            MyHeatmap <- try(LDheatmap(LDsnp, LDdist, LDmeasure = "r", 
                add.map = TRUE, SNP.name = LDsnpName, color = rev(cm.colors(20)), 
                name = "myLDgrob", add.key = TRUE, geneMapLabelY = 0.1))
            if (!inherits(MyHeatmap, "try-error")) {
                grid.edit(gPath("myLDgrob", "Key", "title"), 
                  gp = gpar(cex = 0.5, col = "blue"))
                grid.edit(gPath("myLDgrob", "geneMap", "title"), 
                  gp = gpar(just = c("center", "bottom"), cex = 0.8, 
                    col = "black"))
                grid.edit(gPath("myLDgrob", "geneMap", "SNPnames"), 
                  gp = gpar(cex = 0.3, col = "black"))
            }
            else {
                print("Warning: error in converting genotype. No LD plot!")
            }
            dev.off()
            print("LD heatmap crated")
        }
        else {
            print("Warning: There are less than two SNPs on the region you sepcified. No LD plot!")
        }
    }
    Timmer = GAPIT.Timmer(Timmer = Timmer, Infor = "after LD plot")
    Memory = GAPIT.Memory(Memory = Memory, Infor = "after LD plot")
    return(list(G = G, GD = GD, GI = GI, GT = GT, hasGenotype = hasGenotype, 
        genoFormat = genoFormat, KI = KI, PC = PC, byFile = byFile, 
        fullGD = fullGD, Timmer = Timmer, Memory = Memory, SNP.QTN = SNP.QTN))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
