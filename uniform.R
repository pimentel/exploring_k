library(cluster)
library(data.table)
library(dplyr)
library(kselection)
library(multicore)

# XXX: remember to set the number of core to use
option(cores = 30)

#' @param n sample size
#' @param p the number of dimensions

#' @param rng a matrix with each row describing lower and upper bound for the
#' uniform distribution
#' @return a matrix with n uniform samples over rng in p dimensions. If p is NULL, then uses rng.
uniform_sample <- function(n, p = NULL, rng = matrix(c(0, 1), ncol = 2))
{
    if (!is(rng, "matrix"))
        stop("'rng' must be a matrix!")

    if ( ncol(rng) != 2 )
        stop("'rng' must have two columns (min, max) for each dimension")

    if (is.null(p) && nrow(rng) < 2)
        stop("Must specify 'p' or have more than one row for 'rng'")
    else if (!is.null(p) && nrow(rng) == 1)
        rng <- matrix(rep(rng, p), ncol = 2, byrow = T)
    else
        warning("'p' and 'rng' specified -- ignoring 'p'")


    apply(rng, 1, function(row)
        {
            runif(n, min = row[1], max = row[2])
        })
}

#' @param n_samps the number of samples (observations)
#' @param n_sim the number of simulations to run
#' @param p required parameter sent to uniform_sample
pham_batch <- function(X, n_samps, n_sim,  p, ...)
{
    if (length(n_samps) != 1L || length(n_sim) != 1L)
        stop("Only takes one argument for value 'n_samps'")

    res <- lapply(X, function(tmp_x)
        {
            time_res <- system.time(ks <- kselection(tmp_x, ...))

            list(time_res = as.data.frame(as.list(time_res)), ks = ks)
        })

    res_df <- rbindlist(lapply(1:length(res), function(i)
            {
                t_res <- res[[i]]$time_res
                ks <- res[[i]]$ks

                tmp <- data.frame(n_samps = n_samps,
                    n_sim = n_sim,
                    p = p,
                    f_k = ks$f_k,
                    sim_num = i,
                    k = 1:length(ks$f_k),
                    k_opt = ks$k)

                data.frame(list(tmp, t_res))
            }))

    as.data.frame(res_df)
}


#' @param n_samps_per_sim number of samples per simulation (e.g. generate n_samps_per_sim vectors of dimension dims)
#' @param n_sim how many simulations to run per combination
#' @param dims the number of dimensions
#' @param ... additional parameters to uniforom_sample()
benchmark_pham_gap <- function(n_samps_per_sim, n_sim, dims, ...)
{
    all_comb <- expand.grid(n_samp = n_samps_per_sim, n_dim = dims)

    res <- mclapply(1:nrow(all_comb), function(i)
        {
            cat(sprintf("num samples:\t%d, dim\t%d\n", all_comb$n_samp[i],
                    all_comb$n_dim[i]))

            cur_data <- lapply(1:n_sim, function(x)
                uniform_sample(all_comb$n_samp[i], all_comb$n_dim[i], ...)
                )
            pham_res <- pham_batch(cur_data, all_comb$n_samp[i], n_sim, all_comb$n_dim[i])
            gap_res <- gap_batch(cur_data, all_comb$n_samp[i], n_sim, all_comb$n_dim[i])

            list(pham = pham_res, gap = gap_res)
        })

    pham <- rbindlist(lapply(res, function(x) x$pham))
    gap <- rbindlist(lapply(res, function(x) x$gap))

    list(pham = as.data.frame(pham), gap = as.data.frame(gap))
}

gap_batch <- function(X, n_samps, n_sim,  p, ...)
{
    if (length(n_samps) != 1L || length(n_sim) != 1L)
        stop("Only takes one argument for value 'n_samps' and 'n_sim'")

    res <- lapply(X, clusGap, FUNcluster = kmeans,
        K.max = 15, nstart = 25)

    res <- lapply(X, function(tmp_x)
        {
            time_res <- system.time(cg <- clusGap(tmp_x, FUNcluster = kmeans, K.max = 15, nstart = 25))

            list(time_res = as.data.frame(as.list(time_res)), cg = cg)
        })

    res_df <- rbindlist(lapply(1:length(res), function(i)
            {
                t_res <- res[[i]]$time_res
                cg <- res[[i]]$cg

                gap <- cg$Tab[,"gap"]
                gapSE <- cg$Tab[,"SE.sim"]

                k_opt <- maxSE(gap, gapSE)

                tmp <- data.frame(n_samps = n_samps,
                    n_sim = n_sim,
                    p = p,
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


set.seed(52)
bench <- benchmark_pham_gap(c(50, 100, 250, 500, 1000), 100, c(2, 3, 10, 20))
save(bench, file = "~/pham.RData")


load("~/Dropbox/blog/phamKmeans/pham.RData", verbose = T)

pham_distinct <- bench$pham %>%
    select(-c(f_k, k)) %>%
    distinct() %>%
    mutate(method = "pham")

gap_distinct <- bench$gap %>%
    select(-c(gap, gapSE, k)) %>%
    distinct() %>%
    mutate(method = "gap")

pham_gap_distinct <- rbind_list(pham_distinct, gap_distinct)

ggplot(pham_gap_distinct, aes(k_opt, fill = method), group = method) +
    geom_histogram(position = "dodge") +
    scale_x_continuous(breaks = 1:15) +
    facet_wrap(~ p + n_samps) +
    theme_bw() +
    xlab("k selected") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("~/Dropbox/blog/phamKmeans/img/optimal_k.png")

ggsave("~/Dropbox/blog/phamKmeans/img/optimal_k.pdf")

pham_bench <- pham_distinct %>%
    select(-c(method)) %>%
    group_by(n_samps, p, k_opt) %>%
    summarise(count_pham = n())

gap_bench <- gap_distinct %>%
    select(-c(method)) %>%
    group_by(n_samps, p, k_opt) %>%
    summarise(count_gap = n())

bench_table <- left_join(pham_bench, gap_bench,
    by = c("n_samps", "p", "k_opt"))
bench_table <- left_join(gap_bench, pham_bench, by = c("n_samps", "p", "k_opt")) %>%
    union(bench_table)

bench_table[is.na(bench_table)] <- 0

################################################################################
# looking more closely at pham
################################################################################

neg_profiles <- pham_distinct %>%
    filter(n_samps == 50 & k_opt != 1) %>%
    semi_join(bench$pham, ., by = c("sim_num", "n_samps", "p", "n_sim"))

ggplot(neg_profiles, aes(k, f_k, group = factor(sim_num), colour = factor(sim_num))) +
    geom_point(alpha = 0.3, shape = 3) +
    geom_line(alpha = 0.3, size = 0.7) +
    geom_abline(slope = 0, intercept = 0.85, colour = "red", size = 1, linetype = 3) +
    ylim(0, 1.2) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("f(k)") +
    scale_x_continuous(breaks = 1:15) +
    ggtitle("Incorrectly classified simulations (n = 50, p = 2)")
ggsave("~/Dropbox/blog/phamKmeans/img/neg_profiles.png")

pos_profiles <- pham_distinct %>%
    filter(n_samps == 50 & k_opt == 1) %>%
    slice(11:20) %>%
    semi_join(bench$pham, ., by = c("sim_num", "n_samps", "p", "n_sim"))

ggplot(pos_profiles, aes(k, f_k, group = factor(sim_num), colour = factor(sim_num))) +
    geom_point(alpha = 0.4, shape = 3) +
    geom_line(alpha = 0.4, size = 0.7) +
    geom_abline(slope = 0, intercept = 0.85, colour = "red", size = 1, linetype = 3) +
    ylim(0, 1.2) +
    theme_bw() +
    theme(legend.position = "none") +
    ylab("f(k)") +
    scale_x_continuous(breaks = 1:15) +
    ggtitle("Correctly classified simulations (n = 50, p = 2)")
ggsave("~/Dropbox/blog/phamKmeans/img/pos_profiles.png")

################################################################################
# larger space -- not reported on blog as looks very similar to
################################################################################

set.seed(42)
system.time(
bench_lrg_rng <- benchmark_pham_gap(
    n_samps_per_sim = c(50, 100, 250, 500, 1000),
    n_sim = 100,
    dims = c(2, 3, 10, 20),
    rng = matrix(c(0, 10000), ncol = 2))
)

save(bench_lrg_rng, file = "~/bench_lrg_rng.RData")

load("~/Dropbox/blog/phamKmeans/bench_lrg_rng.RData")

pham_lrg_distinct <- bench_lrg_rng$pham %>%
    select(-c(f_k, k)) %>%
    distinct() %>%
    mutate(method = "pham")

gap_lrg_distinct <- bench_lrg_rng$gap %>%
    select(-c(gap, gapSE, k)) %>%
    distinct() %>%
    mutate(method = "gap")

pham_gap_lrg_distinct <- rbind_list(pham_lrg_distinct, gap_lrg_distinct)

ggplot(pham_gap_lrg_distinct, aes(k_opt, fill = method), group = method) +
    geom_histogram(position = "dodge") +
    facet_wrap(~ p + n_samps) +
    xlab("Optimal k")

################################################################################
# look at timings
################################################################################

pham_timing <- bench$pham %>%
    filter(k == 1) %>%
    select(n_samps, p, k_opt, elapsed) %>%
    mutate(method = "pham")
gap_timing <- bench$gap %>%
    filter(k == 1) %>%
    select(n_samps, p, k_opt, elapsed) %>%
    mutate(method = "gap")

all_timing <- rbind_list(pham_timing, gap_timing)

# this is the figure used in the blog
ggplot(all_timing, aes(elapsed, fill = method), group = method) +
    geom_histogram(position = "stack", binwidth = 3.5) +
    facet_wrap(~ p + n_samps) +
    theme_bw() +
    xlab("Elapsed time (s)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("~/Dropbox/blog/phamKmeans/img/time_hist.png")
ggsave("~/Dropbox/blog/phamKmeans/img/time_hist.pdf")

ggplot(all_timing, aes(method, elapsed)) +
    geom_boxplot() +
    facet_wrap(~ p + n_samps) +
    theme_bw() +
    xlab("Method") +
    ylab("Elapsed time (s)")
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("~/Dropbox/blog/phamKmeans/img/time_boxplot.pdf")


ggplot(all_timing, aes(elapsed, fill = method), group = method) +
    geom_density() +
    facet_wrap(~ p + n_samps) +
    theme_bw() +
    xlab("Elapsed time (s)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave("~/Dropbox/blog/phamKmeans/img/time_density.pdf")
