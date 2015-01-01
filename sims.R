
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
    mat + eps
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

gen_mvnorm_centers(4, c(0, 0), diag(5, 2), 1)

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

    sim_vals
}

sim_d <- function()
{
    sim_c(nclust = 4,
        nobs = c(25, 50),
        mu = rep(0, 10),
        sigma = diag(1.9, 10))
}

sim_d()

sim_e <- function(nsamp, rng = c(-0.5, 0.5))
{
    clust_1 <- matrix(rep(seq(from = rng[1], to = rng[2], length.out = nsamp), 3),
        ncol = 3)
    clust_2 <- matrix(rep(seq(from = rng[1], to = rng[2], length.out = nsamp), 3),
        ncol = 3) + 10

    clust_1 <- clust_1 + rnorm(length(clust_1), mean = 0, sd = 0.1)
    clust_2 <- clust_2 + rnorm(length(clust_2), mean = 0, sd = 0.1)
    rbind(clust_1, clust_2)
}

plot(sim_e(100))

overlap_clust <- function(delta, nsamp)
{
    clust_1 <- matrix(rnorm(nsamp[1] * 2), ncol = 2)
    clust_2 <- matrix(rnorm(nsamp[2] * 2), ncol = 2)
    clust_2[,1] <- clust_2[,1] + delta
    rbind(clust_1, clust_2)
}

plot(overlap_clust(5, c(50, 50)))
