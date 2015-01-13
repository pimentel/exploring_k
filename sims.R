options(cores = 43)
options(mc.cores = 43)

#' Simulation 'b' in the gap statistic paper
#'
#' By default, simulate three clusters with standard normal error
#'
#' @param counts a vector of counts for each center
#' @param centers a matrix of centers. Each row contains a p-dimensional vector
#' which represents the center
#' @param mu the mean for the standard noise
#' @param sdev the standard deviation for the noise
#' @return a matrix of dimension sum(counts) by ncol(centers)
sim_b <- function(
    counts = c(25, 50, 25),
    centers = matrix(c(0, 0, 0, 5, 5, -3), ncol = 2, byrow = T),
    mu = 0, sdev = 1)
{
    if (length(counts) != nrow(centers))
        stop("Counts must have same length as number of rows.")

    eps <- rnorm(sum(counts) * ncol(centers), mean = mu, sd = sdev)

    mat <- Map(function(row_i, cnt)
        {
            cur_row <- centers[row_i,]
            matrix(rep(cur_row, cnt), ncol = ncol(centers), byrow = T)
        }, 1:nrow(centers), counts)

    mat <- do.call(rbind, mat)
    mat <- mat[sample.int(nrow(mat)),] + eps

    list(data = mat, nclust = nrow(centers), ndim = ncol(centers))
}

#' Replicate centers
#'
#' Replicate a bunch of centers
#'
#' @param centers a matrix where each row represents a center in ncol(centers)
#' dimensional space
#' @param counts a count of each row
#' @return a matrix of \code{sum(counts)} rows and ncol(centers) columns
rep_centers <- function(centers, counts)
{
    mat <- Map(function(row_i, cnt)
        {
            cur_row <- centers[row_i,]
            matrix(rep(cur_row, cnt), ncol = ncol(centers), byrow = T)
        }, 1:nrow(centers), counts)

    do.call(rbind, mat)
}

#' Generate random multivariate normal centers
#'
#' Given some parameters, generate some multivariate normal centers
#'
#' @param nclust then umber of clusters, or centers
#' @param mu a vector of the centers
#' @param sigma the covariance matrix
#' @param min_dist the pairwise euclidean distance of all clusters must be at
#' least this distance. If not, then regenerate the centers.
#' @param max_iter the max number of iterations to try to generate centers. If
#' reaches this max, then fail.
#' @return
gen_mvnorm_centers <- function(nclust, mu, sigma, min_dist, max_iter = 1000)
{
    stopifnot(length(mu) == nrow(sigma) && nrow(sigma) == ncol(sigma))

    cur_iter <- 0
    while(cur_iter < max_iter) {
        centers <- MASS::mvrnorm(nclust, mu, sigma)

        if (all(dist(centers) >= min_dist))
            return(centers)

        cur_iter <- cur_iter + 1
    }

    stop("Couldn't generate centers with min_dist by max_iter.")
}

#' Simulation "c" in the gap statistic paper
#'
#' @param nclust number of clusters
#' @param nobs the number of possible observations in each cluster
#' @param mu a vector of centers
#' @param sigma covariance matrix
#' @return a matrix of
sim_c <- function(
    nclust = 4,
    nobs = c(25, 50),
    mu = rep(0, 4),
    sigma = diag(5, 4)
    )
{
    clust_size <- sample(nobs, nclust, replace = T)
    centers <- gen_mvnorm_centers(nclust, mu, sigma, 1)
    sim_vals <- rep_centers(centers, clust_size)
    eps <- rnorm(length(sim_vals))
    sim_vals <- sim_vals + eps
    sim_vals <- sim_vals[sample.int(nrow(sim_vals)),]

    list(data = sim_vals, clust_size = clust_size, nclust = length(clust_size),
        ndim = ncol(sim_vals))
}

sim_d <- function(nobs = c(25, 50))
{
    sim_c(nclust = 4,
        nobs = nobs,
        mu = rep(0, 10),
        sigma = diag(1.9, 10))
}


sim_e <- function(nsamp, rng = c(-0.5, 0.5))
{
    clust_1 <- matrix(rep(seq(from = rng[1], to = rng[2], length.out = nsamp), 3),
        ncol = 3)
    clust_2 <- matrix(rep(seq(from = rng[1], to = rng[2], length.out = nsamp), 3),
        ncol = 3) + 10

    clust_1 <- clust_1 + rnorm(length(clust_1), mean = 0, sd = 0.1)
    clust_2 <- clust_2 + rnorm(length(clust_2), mean = 0, sd = 0.1)

    res <- rbind(clust_1, clust_2)
    res <- res[sample.int(nrow(res)),]

    list(data = res, nclust = 2, ndim = ncol(clust_1))
}

overlap_clust <- function(delta, nsamp)
{
    clust_1 <- matrix(rnorm(nsamp[1] * 2), ncol = 2)
    clust_2 <- matrix(rnorm(nsamp[2] * 2), ncol = 2)
    clust_2[,1] <- clust_2[,1] + delta

    res <- rbind(clust_1, clust_2)
    res <- res[sample.int(nrow(res)),]

    list(data = , nclust = 2, ndim = ncol(clust_1))
}

#' @param n_sim how many simulations to run per combination
#' @param dims the number of dimensions
#' @param ... additional parameters to sim_fun()
benchmark_all_generic <- function(n_sim, sim_fun, ...)
{
    data <- lapply(1:n_sim, function(x) sim_fun(...)$data)

    res <- list()
    cat("Pham\n")
    res$pham <- pham_batch_generic(data)
    cat("Gap\n")
    res$gap <- gap_batch_generic(data)
    cat("Pred\n")
    res$pred_strength <- pred_batch_generic(data)

    res
}

#' @param n_samps the number of samples (observations)
#' @param n_sim the number of simulations to run
#' @param p required parameter sent to uniform_sample
pham_batch_generic <- function(X, ...)
{
    res <- mclapply(X, function(tmp_x)
        {
            time_res <- system.time(ks <- kselection(tmp_x,
                    max_centers = 10,
                    iter.max = 1000,
                    nstart = 50, ...))

            list(time_res = as.data.frame(as.list(time_res)), ks = ks)
        })

    res_df <- rbindlist(lapply(1:length(res), function(i)
            {
                t_res <- res[[i]]$time_res
                ks <- res[[i]]$ks

                tmp <- data.frame(
                    f_k = ks$f_k,
                    sim_num = i,
                    k = 1:length(ks$f_k),
                    k_opt = ks$k)

                data.frame(list(tmp, t_res))
            }))

    as.data.frame(res_df)
}


gap_batch_generic <- function(X, ...)
{
    res <- mclapply(X, function(tmp_x)
    # res <- lapply(X, function(tmp_x)
        {
            time_res <- system.time(cg <- clusGap(tmp_x,
                    FUNcluster = kmeans,
                    K.max = 10,
                    iter.max = 1000,
                    nstart = 50))

            list(time_res = as.data.frame(as.list(time_res)), cg = cg)
        })

    res_df <- rbindlist(lapply(1:length(res), function(i)
            {
                t_res <- res[[i]]$time_res
                cg <- res[[i]]$cg

                gap <- cg$Tab[,"gap"]
                gapSE <- cg$Tab[,"SE.sim"]

                k_opt <- maxSE(gap, gapSE)

                tmp <- data.frame(
                    gap = gap,
                    gapSE = gapSE,
                    sim_num = i,
                    k = 1:length(gap),
                    k_opt = k_opt
                    )

                data.frame(list(tmp, t_res))
            }
        ))

    as.data.frame(res_df)
}

pred_batch_generic <- function(X, ...)
{
    res <- mclapply(X, function(tmp_x)
        {
            time_res <- system.time(ps <- prediction.strength(tmp_x,
                    Gmax = 10,
                    runs = 50,
                    iter.max = 1000
                    ))

            list(time_res = as.data.frame(as.list(time_res)), ps = ps)
        })

    res_df <- rbindlist(lapply(1:length(res), function(i)
            {
                t_res <- res[[i]]$time_res
                ps <- res[[i]]$ps

                k_opt <- ps$optimalk

                tmp <- data.frame(
                    mean_pred = ps$mean.pred,
                    sim_num = i,
                    k = 1:length(ps$mean.pred),
                    k_opt = k_opt
                    )

                data.frame(list(tmp, t_res))
            }
        ))

    as.data.frame(res_df)
}

################################################################################
# testing benchmark
################################################################################

b_counts <- lapply(c(1, 2, 5, 100), function(mult) mult * c(25, 50, 25))
set.seed(42)
sim_b_res <- lapply(b_counts, function(cnts)
    {
        cat("Counts: ", cnts, "\n")
        benchmark_all_generic(100, sim_b, counts = cnts)
    })
save(sim_b_res, b_counts, file = "sim_res.RData")

c_counts <- lapply(c(1, 2, 5, 100), function(mult) mult * c(25, 50))
set.seed(42)
sim_c_res <- lapply(c_counts, function(cnts)
    {
        cat("Counts: ", cnts, "\n")
        benchmark_all_generic(100, sim_c, nobs = cnts)
    })
save(sim_c_res, c_counts,
    sim_b_res, b_counts,
    file = "sim_res.RData")

d_counts <- lapply(c(1, 2, 5, 100), function(mult) mult * c(25, 50))
set.seed(42)
sim_d_res <- lapply(d_counts, function(cnts)
    {
        cat("Counts: ", cnts, "\n")
        benchmark_all_generic(100, sim_d, nobs = cnts)
    })

save(sim_d_res, d_counts,
    sim_c_res, c_counts,
    sim_b_res, b_counts,
    file = "sim_res.RData")

e_counts <- ceiling(c(200, 500, 1000, 15000, 30000) / 2)
set.seed(42)
sim_e_res <- lapply(e_counts, function(cnts)
    {
        cat("Counts: ", cnts, "\n")
        benchmark_all_generic(100, sim_e, nsamp = cnts)
    })

save(
    sim_e_res, e_counts,
    sim_d_res, d_counts,
    sim_c_res, c_counts,
    sim_b_res, b_counts,
    file = "sim_res.RData")
