library(dplyr)
library(readr)
library(tibble)
library(edgeR)
library(ggplot2)

#library(org.Sc.sgd.db)
setwd("~/github/2020-pn-dm/")

# determine whether there are gene expression differences based on rootstock


# functions ---------------------------------------------------------------

orf_to_common <- function(keys){
  common <- AnnotationDbi::select(org.Sc.sgd.db::org.Sc.sgd.db,
                                  keys = keys,
                                  columns=c("COMMON"),
                                  keytype="ORF")
  return(common)
}

orf_to_description <- function(keys){
  description <- AnnotationDbi::select(org.Sc.sgd.db::org.Sc.sgd.db,
                                       keys = keys,
                                       columns=c("DESCRIPTION"),
                                       keytype="ORF")
  return(description)
}

dc_contrast <- function(dc_fit, coeff, adj.p.filt = 0.01, filter_logFC = T){
  # Test diffusion component coefficient
  dc_con <- contrasts.fit(dc_fit, coefficients = coeff)   
  dc_b <- eBayes(dc_con)
  dc_top <- topTable(dc_b, sort.by = "logFC", number = Inf)
  dc_top <- dc_top[dc_top$adj.P.Val < adj.p.filt, ]
  
  # Translate ORF to Common
  dc_top$ORF <- rownames(dc_top)
  commons <- orf_to_common(dc_top$ORF)
  # subset to first occurrence of ORF
  commons <- commons[match(unique(commons$ORF), commons$ORF), ]
  commons <- commons[ , -2]
  commons$COMMON <- ifelse(is.na(commons$COMMON), commons$ORF, commons$COMMON)
  commons <- unique(commons[ , ])
  dc_top <- left_join(dc_top, commons, by = c("ORF"))
  
  # Add SGD Description 
  descriptions <- orf_to_description(dc_top$ORF)
  # subset to first occurrence of ORF
  descriptions <- descriptions[match(unique(descriptions$ORF), descriptions$ORF), ]
  descriptions <- descriptions[ , -2] # remove sgd id
  # descriptions$DESCRIPTION <- ifelse(is.na(descriptions$COMMON), descriptions$ORF, descriptions$DESCRIPTION)
  # descriptions <- unique(descriptions[ , ])
  dc_top <- left_join(dc_top, descriptions, by = c("ORF"))
  
  # normalize logfc based on min and max of DC
  normalizer <- max(dc_fit$design[ , coeff]) - min(dc_fit$design[ , coeff])
  dc_top <- dc_top %>%
    mutate(logFC_adj = logFC*normalizer) %>%
    arrange(logFC_adj) %>%
    dplyr::select(logFC, logFC_adj, AveExpr, t, P.Value, adj.P.Val, B, ORF, COMMON, DESCRIPTION)
  if(filter_logFC == F) {
    return(dc_top)
  } else {
    # Join Common to diffex results
    dc_top <- dc_top %>%
      filter(abs(logFC_adj) > 2) 
    return(dc_top)
  }
}

# read in metadata --------------------------------------------------------

info <- read_csv("samples_ava.csv")
info$ava_id <- factor(info$ava_id, levels = c('OR1', 'OR2', 'AV1', 'AV2', 
                                              'SNC1', 'SNC2', 'RRV1', 'RRV2', 'RRV3',
                                              'CRN1', 'AS1', 'AS2', 'SMV1', 
                                              'SMV2','SRH1'))
info$ava <- factor(info$ava, levels = c('OR', 'AV', 'SNC', 'RRV', 'CRN', 
                                        'AS', 'SMV','SRH'))
# table(info$rootstock)
# 101_14          3309C riparia_gloire 
# 196             60             40

# table(info$rootstock, info$ava_id)
#                 OR1 OR2 AV1 AV2 SNC1 SNC2 RRV1 RRV2 RRV3 CRN1 AS1 AS2 SMV1 SMV2 SRH1
# 101_14           0   0  20   0   20   20   20    0    0   20  18  18   20   20   20
# 3309C            0   0   0  20    0    0    0   20   20    0   0   0    0    0    0
# riparia_gloire  20  20   0   0    0    0    0    0    0    0   0   0    0    0    0

rootstock_101_14 <- c("AV1", "SNC1", "SNC2", "RRV1", "CRN1", "AS1", "AS2", "SMV1", "SMV2", "SRH1")
rootstock_other <- c("OR1", "OR2", "AV2", "RRV2", "RRV3")

# read in counts ----------------------------------------------------------

counts19 <- read_tsv("v2019/outputs/counts/raw_counts.tsv") %>%
  dplyr::filter(grepl(pattern = "^Y", x = gene)) %>%
  dplyr::select(-`2019_inoculum_SRH1_readcounts.txt`,
                -`2019_inoculum_AS1_readcounts.txt`, 
                -`2019_inoculum_SMV2_readcounts.txt`,
                -`2019_inoculum_AS2_readcounts.txt`,
                -`2019_inoculum_RRV1_readcounts.txt`,
                -`2019_inoculum_SMV1_readcounts.txt`) %>%
  as.data.frame() %>%
  column_to_rownames("gene")
colnames(counts19) <- gsub("_readcounts.txt", "", colnames(counts19))


# DGE -- brix 2019 --------------------------------------------------------

d0_19 <- DGEList(brix_counts19)                    # initialize the DGE object
d19 <- calcNormFactors(d0_19)
brix19 <- info %>%
  dplyr::filter(year == 2019) %>%
  dplyr::filter(!is.na(brix))
brix_mm19 <- model.matrix(~brix19$brix)
brix_y19 <- voom(d19, brix_mm19, plot = F)
brix_fit19 <- lmFit(brix_y19, brix_mm19)
brix_19 <- dc_contrast(dc_fit = brix_fit19, coeff = 2, filter_logFC = F) # Test brix coefficient

# DGE -- brix 10 2019 ------------------------------------------------------

brix_counts19_10  <- counts19[ , 1:150] %>%
  select(contains(rootstock_101_14))
d0_19_10 <- DGEList(brix_counts19_10)                    # initialize the DGE object
d19_10 <- calcNormFactors(d0_19_10)
brix19_10 <- info %>%
  dplyr::filter(year == 2019) %>%
  dplyr::filter(!is.na(brix)) %>%
  dplyr::filter(id %in% colnames(brix_counts19_10))
brix_mm19_10 <- model.matrix(~brix19_10$brix)
brix_y19_10 <- voom(d19_10, brix_mm19_10, plot = F)
brix_fit19_10 <- lmFit(brix_y19_10, brix_mm19_10)
brix_19_10 <- dc_contrast(dc_fit = brix_fit19_10, coeff = 2, filter_logFC = F) # Test brix coefficient

# DGE -- brix 5 2019 ------------------------------------------------------

brix_counts19_5  <- counts19[ , 1:150] %>%
  select(contains(rootstock_other))
d0_19_5 <- DGEList(brix_counts19_5)                    # initialize the DGE object
d19_5 <- calcNormFactors(d0_19_5)
brix19_5 <- info %>%
  dplyr::filter(year == 2019) %>%
  dplyr::filter(!is.na(brix)) %>%
  dplyr::filter(id %in% colnames(brix_counts19_5))
brix_mm19_5 <- model.matrix(~brix19_5$brix)
brix_y19_5 <- voom(d19_5, brix_mm19_5, plot = F)
brix_fit19_5 <- lmFit(brix_y19_5, brix_mm19_5)
brix_19_5 <- dc_contrast(dc_fit = brix_fit19_5, coeff = 2, filter_logFC = F) # Test brix coefficient

all_brix <- brix_19 %>%
  select(ORF, logFC_adj_19 = logFC_adj) %>%
  left_join(brix_19_10, by = c("ORF")) %>%
  select(ORF, logFC_adj_19, logFC_adj_19_10 = logFC_adj) %>%
  left_join(brix_19_5, by = c("ORF")) %>%
  select(ORF, logFC_adj_19, logFC_adj_19_10, logFC_adj_19_5 = logFC_adj)
View(all_brix)
dim(all_brix)
summary(lm(logFC_adj_19_10 ~ logFC_adj_19_5, data = all_brix))

ggplot(all_brix, aes(x = logFC_adj_19_10, y = logFC_adj_19_5)) +
  geom_point() +
  theme_minimal()
