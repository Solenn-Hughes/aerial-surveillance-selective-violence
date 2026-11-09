# 01: transparent support, correlation, ELI overlap and residualised condition diagnostics.
source('C:/Users/USERNAME/Downloads/JPR_H2_LWJ_00_COMMON.R')
OUT1 <- file.path(OUT,'01_support')
dir.create(OUT1,recursive=TRUE,showWarnings=FALSE)
S <- unique(D[,.(cell_id,district_id_for_h2,eli_final,lwj_contested_share,lwj_taliban_control_share,lwj_goa_control_share,lwj_unconfirmed_share)])
if (nrow(S)!=241L) stop('ELI or LWJ measures not constant within cell; stop and investigate')
S[, group:=fifelse(lwj_contested_share>=0.6,'contested_majority',fifelse(lwj_goa_control_share>=0.6,'government_majority','mixed_or_other'))]
stopifnot(nrow(S)==241L)
q <- S[,.(cells=.N, districts=uniqueN(district_id_for_h2),eli_min=min(eli_final),eli_p10=quantile(eli_final,.1),eli_med=median(eli_final),eli_p90=quantile(eli_final,.9),eli_max=max(eli_final)),by=group]
fwrite(q,file.path(OUT1,'01_groups_eli_support.csv'))
cat('\nELI distribution by pre-Doha majority status:\n'); print(q)
M <- as.matrix(S[,.(eli_final,lwj_contested_share,lwj_taliban_control_share,lwj_goa_control_share,lwj_unconfirmed_share)])
fwrite(as.data.table(cor(M)),file.path(OUT1,'02_cell_correlations.csv'))
cat('\nCell-level correlations:\n'); print(round(cor(M),3))
G <- S[group=='contested_majority',eli_final]; H <- S[group=='government_majority',eli_final]
if (length(G)>=10L && length(H)>=10L) {
  lo <- max(unname(quantile(G,.1)),unname(quantile(H,.1)))
  hi <- min(unname(quantile(G,.9)),unname(quantile(H,.9)))
  S[,overlap_10_90:=eli_final>=lo & eli_final<=hi]
  ov <- data.table(eli_overlap_min=lo,eli_overlap_max=hi, overlap_nonempty=lo<=hi,
                   cells_in_overlap=if(lo<=hi)sum(S$overlap_10_90) else 0L,
                   government_majority_in_overlap=if(lo<=hi)sum(S$overlap_10_90 & S$group=='government_majority') else 0L,
                   contested_majority_in_overlap=if(lo<=hi)sum(S$overlap_10_90 & S$group=='contested_majority') else 0L)
  fwrite(ov,file.path(OUT1,'03_eli_overlap.csv')); print(ov)
} else {cat('Too few majority-status cells to calculate overlap.\n')}
# Compositional shares and ELI have only cell-level variation; no claim of full causal overlap.
fwrite(S,file.path(OUT1,'04_cell_support_records.csv'))
cat('Support files: ',OUT1,'\n')
