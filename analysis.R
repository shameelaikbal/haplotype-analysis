#population structure
#convert vcf  to PLINK binary
plinkCommand = "plink --allow-extra-chr --vcf extracted_data.vcf.gz --const-fid --make-bed --out chr"
system(plinkCommand)

#generate eigen values
plink --bfile chr  --pca --out pca_results

#GWAS with rMVP data preparation
library(rMVP)
MVP.Data(fileVCF="data/SNPs_lee.id.biallic_maf_0.05_geno_0.1_916.vcf",
         filePhe="data/pheno.txt",
         fileKin=FALSE,
         filePC=FALSE,
         out="GWAS_lee/lee_mvp") # create a folder such as lee_mvp to store the results

#association analysis
genotype <- attach.big.matrix("lee_mvp.geno.desc")
phenotype <- read.table("lee_mvp.phe",head=TRUE)
map <- read.table("lee_mvp.geno.map" , head = TRUE)


fMVP <- MVP(
    phe=phenotype,
    geno=genotype,
    map=map,
    #K=Kinship,
    #CV.GLM=Covariates,     ##if you have additional covariates, please keep there open.
    #CV.MLM=Covariates,
    #CV.FarmCPU=Covariates,
    nPC.GLM=5,      ##if you have added PC into covariates, please keep there closed.
    nPC.MLM=3,
    nPC.FarmCPU=3,
    priority="memory",       ##for Kinship construction
    #ncpus=10,
    vc.method="BRENT",      ##only works for MLM
    maxLoop=10,
    method.bin="static",      ## "FaST-LMM", "static" (#only works for FarmCPU)
    #permutation.threshold=TRUE,
    #permutation.rep=100,
    threshold=0.05,
    method=c("GLM", "MLM", "FarmCPU"),
    out="GWAS_result" #create this folder to store the results
)

#delimit a region of interest around the significant SNP from GWAS using linkage disequilibrium for local haplotyping
#convert chr 20 VCF to plink binary
plink --allow-extra-chr --vcf SNPs_lee.vcf --const-fid --make-bed --out chr20

#extract a 200kb region around the significant SNP
plink  --bfile chr20 --allow-extra-chr  --snp Gm20_438877528 --window 200  --make-bed --out region_200kb_chr20_43877528

#convert the extracted 200kb region back to VCF
plink --allow-extra-chr --bfile region_200kb_chr20_43877528 --recode vcf-iid --out region_200kb_chr20_43877528

#Visualise the LD in this region using LDBlockShow
./bin/LDBlockShow -InVCF region_200kb_chr20_43877528.vcf -SeleVar 2 -BlockType 2  -Region Gm20:43778602-43977430 -OutPng -OutPut region_200kb_chr20_43877528


#prepare an LD matrix for the final delimited region
plink --allow-extra-chr --r2 square --vcf region_200kb_chr20_43877528.vcf

#inputs to crosshap for haplotype analysis
#the delimited VCF obtained from above step
#make sure VCF has only biallelic SNPs
vcf_1 <- read_vcf ("region_200kb_chr20_43877528.vcf")

#input the phenotype data 
pheno_1 <-  read_pheno ("pheno.txt")

#input metadata
meta_1 <- read_metadata("metadata__crosshap.txt")

# input the prepared LD matrix plink.ld
LD_1 <- read_LD("plink_chr20.ld", vcf = vcf_1)

#Add minimum marker group for SNP count
MGmin_gm <- 25

##Add list of epsilon values to run haplotyping on
epsilon_gm <- c(0.2, 0.4, 0.6, 0.8, 1)


#Run the haplotyping at all provided epsilon value
run_haplo_1 <- run_haplotyping (vcf = vcf_1,
                                 LD = LD_1,
                                 pheno = pheno_1,
                                 metadata = meta_1,
                                 epsilon = epsilon_gm ,
                                 MGmin = MGmin_gm)


#Provide phenotype data and parameters used to create haplotype objects
#Add type = 'MG' to ensure it summarizes Marker Groups rather than haplotypes and type hap to visualise haplotypes
hap_clustree_chr1 <- clustree_viz(HapObject = run_haplo_1,
                                   type = 'hap')

#Add type = 'MG' to ensure it summarizes Marker Groups rather than haplotypes
hap_clustree_chr21_mg <- clustree_viz(HapObject = run_haplo_1,
                                      type = 'MG')

#choosing epsilon 0.6
#Visualize haplotype object created by run_haplotyping() at the chosen epsilon (0.6)
hap_viz_chr1 <- crosshap_viz(HapObject = run_haplo_1, epsilon = 0.6)

#visualization with the position of SNPs in marker groups
hap_viz_chr20_pos <- crosshap_viz(HapObject = run_haplo_1, epsilon = 0.6, plot_left = "pos")

#statistical analysis
#To analyze significant differences in flowering time among distinct subpopulations, Welch’s two-sample t-test was used
# Load Data
data <- read.csv("welshAB.csv") # welshAB.csv has 
coloumn_A <- data$A
coloumn_B <-  data$B

# Perform Welch’s Two Sample t-test
result <- t.test(coloumn_A, coloumn_B, var.equal = FALSE) # var.equal = FALSE implies that it's a Welch’s t-test

