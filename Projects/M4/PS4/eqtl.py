from pca import *

## 1. Select most likely SNPs
##  -search radius: only look at SNPs that are close to genes (don't have this data)
##  -minor allele frequency: the frequency at which the minor allele occurs in pop. 
##      can use this as a cut-off (MAF < 0.05 or 5%)
## 2. Covariate adjustments (not done here since no covariate data is provided)
##  -use PCA to separate effects of genotype from confounding factors (biological or
##      population factors)
## 3. Linear regression model
##  -fit data to model Y = B0 + B1*X + E, where E is noise
##  -if you have covariate data, Y = B0 + B1*X + B2*PC1 + B3*PC2 + ... + Bn+1*PCn + E
## 4. Statistical significance
##  -find the p value for each eQTL to determine significance
##  -hypthoses, B0 == 0, B1 != 0

def applyMAF(genotypes, mafthreshold=0.05):

    # assume genotypes is mxn matrix, 
    # m = number patients
    # n = number snps
    
    snpIds = []

    for snp in range(genotypes.shape[1]):
        macount = 0
        for patient in range(genotypes.shape[0]):
            macount += genotypes[patient][snp]
        maf = float(macount) / (2*genotypes.shape[0])
        if maf < mafthreshold:
            snpIds.append(snp)
        
    return snpIds

if __name__ == "__main__":

    snpdata, snps, patientsSnp = readGenotypeData("SnpData.txt")
    exprdata, genes, patientsExpr = readExpressionData("ExpData.txt")

    snpIds = applyMAF(snpdata)

    print("Number of SNPs above threshold: ", len(snpIds))