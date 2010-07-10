# Transpose table before call to in row major order
FlowSPD.cluster <- function(tbl, k) {
    .Call("FSPD_cluster",t(tbl),as.integer(k))
}
