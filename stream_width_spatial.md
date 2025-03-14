Width class (spatial model)
================

``` r
library(sf)
```

    ## Linking to GEOS 3.12.1, GDAL 3.8.4, PROJ 9.3.1; sf_use_s2() is TRUE

``` r
library(sfnetworks)
library(SSN2)
```

    ## Warning: package 'SSN2' was built under R version 4.4.1

``` r
library(ggplot2)
```

Here is an idea about how we might be able to model width and width
class where we have both point observations of width and depth, and
reach-level observations of width class. Below, I will simulate what
this data might look like and fit a model to it.

We model log(width) $log(w_{s,r})$ at location $s$ in reach $r$ as a
spatially continuous variable. Width class $C_r$ is a function of all
the $w_{s,r}$ over reach $r$.

\$\$

\$\$

Load stream network. Use a small area so that computations aren’t too
slow.

``` r
# Fit ssn model to get covariance matrix 
library(data.table)
if (!file.exists("streams_clipped/nodes_of_interest.shp")) {
  nodes <- read_sf("streams_clipped/nodes_testc.shp")
  st_crs(nodes) <- "EPSG:26911"
  aoi <- read_sf("ssn_files_keep/ssn_test_region.shp")
  aoi <- st_transform(aoi, st_crs(nodes))
  n <- st_intersection(nodes, aoi)
  write_sf(n, "streams_clipped/nodes_of_interest.shp")
  rm(nodes)
} else {
  n <- read_sf("streams_clipped/nodes_of_interest.shp")
}


setDT(n)
n$nid <- 1:nrow(n)
n$to_nid <- 0
# make table

# This always hangs... 
nothing <- sapply(1:nrow(n), \(x) {
  t = n[NodeNum == n[nid == x]$ToNode]$nid
  if (length(t) == 0) {t = NA}
  n[nid == x, to_nid := t]
  return(NULL)
})
# edges
e <- n[, .(from = nid, to = to_nid)]
e <- e[!is.na(to)]

# For now just let length = 1 (evenly spaced nodes)
e$length <- 1
```

Create covariance matrix based only on downstream distance.

``` r
library(Matrix)
# set up

N <- nrow(n)
# Use contributing area for weights

## Parameters
# sigma <- 1.3 # variance parameter for x
# theta <- .04 # autocorrelation parameter
# alpha = 0 # intercept
# beta <- .1 # effect of x on width
# phi <- 20 # variance parameter for y

# simulate multivariate normal from precision matrix
rmvnorm_prec <-
function( mu, # estimated fixed and random effects
          prec, # estimated joint precision
          n.sims = 1) {

  require(Matrix)
  # Simulate values
  z0 = matrix(rnorm(length(mu) * n.sims), ncol=n.sims)
  # Q = t(P) * L * t(L) * P
  L = Cholesky(prec, super=TRUE)
  # Calcualte t(P) * solve(t(L)) * z0 in two steps
  z = solve(L, z0, system = "Lt") # z = Lt^-1 * z
  z = solve(L, z, system = "Pt") # z = Pt    * z
  return(mu + as.matrix(z))
}

make_stream_cov <- function(edges, 
                            nodes, 
                            theta = .04, # autocorrelation parameter
                            prec = FALSE) {
  
  N <- nrow(nodes)
  # expected covariance
  # weights matrix
  A <- sparseMatrix(i=edges$from, 
                    j=edges$to, 
                    x=1, 
                    dims=c(N,N) )
  sources <- which(colSums(A) < 1)
  confluences <- which(colSums(A) > 1)
  weights_mx <- A
  for (c in confluences) {
    # from 
    f <- which(A[, c] == 1)
    w <- nodes[f, AREA_SQKM]/(sum(nodes[f, AREA_SQKM]))
    weights_mx[f, c] <- w
  }
  # autocorrelation for each path
  rho_mx <- sparseMatrix(i=edges$from, 
                         j=edges$to, 
                         x=1, dims=c(N,N) )
  
  # variance contribution to each x
  var_mx <- sparseMatrix(i=edges$from, 
                         j=edges$to, 
                         x=(1-exp(-2*theta*edges$length)), dims=c(N,N) )
  
  # path matrix: weights times autocorrelation (element-wise, not mx multiplication)
  Gamma <- rho_mx * weights_mx
  
  # variance: weighted average for each x
  v <- apply((weights_mx * var_mx), 2, sum)
  # set variance for initial conditions
  v[sources] <- 1
  V <- diag(x=v)
  v_inv <- diag(x = sqrt(v))
  
  # construct precision matrix
  
  I <- diag(N)
  if (prec) {
    return(Matrix((I-Gamma) %*% v_inv %*% t(I-Gamma), sparse = TRUE))
  } else {
    return(Matrix(t(solve(I-Gamma)) %*% V %*% (solve(I-Gamma)), sparse = TRUE))
  }
}

# Precision (inverse covariance matrix)
# which is much sparser and faster.
Q <- make_stream_cov(e, n, prec = TRUE)
image(Q)
```

![](stream_width_spatial_files/figure-gfm/cov-mx-1.png)<!-- -->

``` r
# covariance matrix will also be sparse because it is only nonzero for 
# nodes which are flow-connected. 
```

Simulate width at each point based on covariance matrix.

``` r
set.seed(10161994)

# Two spatial random effects. One which changes the coefficient for 
# contributing area and one which changes the overall mean. 
# Each can have different levels of spatial autocorrelation and variance. 
# It will likely be hard to estimate both but let's see!

# random spatial effect which changes the coefficient for area. 
psi_area <- rmvnorm_prec(rep(0, dim(Q)[1]), Q)
# random spatial effect which changes the mean. 
psi_mean <- rmvnorm_prec(rep(0, dim(Q)[1]), Q)
n$psi_area <- psi_area
n$psi_mean <- psi_mean
n$mu <- -2 + (1.3 + .02*psi_area)*n$AREA_SQKM + .01*psi_mean

phi <- .1
n$w <- rgamma(length(n$mu), 1/phi, , exp(n$mu)*phi)
plot(st_as_sf(n)[, 'AREA_SQKM'], pch = 3, cex = .5)
```

![](stream_width_spatial_files/figure-gfm/simulate-width-1.png)<!-- -->

``` r
plot(st_as_sf(n)[, 'psi_area'], pch = 3, cex = .5)
```

![](stream_width_spatial_files/figure-gfm/simulate-width-2.png)<!-- -->

``` r
plot(st_as_sf(n)[, 'psi_mean'], pch = 3, cex = .5)
```

![](stream_width_spatial_files/figure-gfm/simulate-width-3.png)<!-- -->

``` r
plot(st_as_sf(n)[, 'w'], pch = 3, cex = .5)
```

![](stream_width_spatial_files/figure-gfm/simulate-width-4.png)<!-- -->

Simulate width class at each reach based on simulated widths.  
First get reaches. For now, just use links.

``` r
A <- sparseMatrix(i=e$from, 
                    j=e$to, 
                    x=1, 
                    dims=c(N,N) )

# start with only sources, outlets, and confluences
sources <- which(colSums(A) == 0)
outlets <- which(rowSums(A) == 0)
confluences <- which(colSums(A) == 2)

new_n <- data.table(old_id = c(sources, confluences))
setorder(new_n, old_id)
new_n$new_id <- 1:nrow(new_n)
new_e <- data.table(from = new_n$new_id)
new_n <- rbind(new_n, data.table(old_id = outlets, new_id = 1:length(outlets) + nrow(new_n)))

get_downstream <- function(id, i = 1, n_list = id) {
  to = n[id, to_nid]
  n_list <- c(n_list, to)
  i = i+1 
  # id <- to
  # to %in% c(confluences, outlets)
  if (to %in% c(confluences, outlets)) {
    return(list(to = new_n[old_id == to, new_id], i = i, n_list = n_list))
  } else {
    get_downstream(to, i, n_list)
  }
}
sets <- lapply(new_n[new_e$from, old_id], get_downstream)

n_sf <- st_as_sf(n)
new_e$geometry <- lapply(sets, \(s) {st_linestring(st_coordinates(n_sf)[s$n_list, ])}) |> 
  st_as_sfc()
new_e$to <- sapply(sets, \(s) s$to)
new_e$rid <- 1:nrow(new_e)
new_e <- st_as_sf(new_e)

# Add new edge identifier to old nodes
```
