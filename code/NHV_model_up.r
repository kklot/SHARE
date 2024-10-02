sessionInfo()

here::i_am("code/NHV_model_up.r")

library(INLA)
library(dplyr)
library(tidyr)
library(knitr)
# devtools::install('~/Code/R/ktools/') # remotes::install_github('kklot/ktools')
library(ktools)
library(here)
library(ggplot2)
library(haven)
library(spdep)
library(sf)
library(rnaturalearth)
library(countrycode)
sf::sf_use_s2(FALSE)

d <- readRDS(here("data/dNi_up.rds")) |>
    filter(country != "Ireland") |>
    filter(country != "Israel") |>
    select(-ph006dno, -ph006dot) |>
    select(-paste0("ph006d", c(15, 16, 18, 19, 20, 21))) |> # later collected diseases
    mutate(
        across(where(is.labelled), as.double),
        across(starts_with("ph006d"), as.numeric),
        age = as.numeric(age),
        year = as.numeric(int_year),
        yob = year - age,
        sex = case_when(gender.x == "Male" ~ 0, gender.x == "Female" ~ 1, TRUE ~ NA_real_)
    ) |>
    filter(newedu <= 25) |>
    filter(age >= 50, age <= 100) |>
    select(-int_year, -gender.x)

d <- d |> filter(!is.na(ph006d1))

pm_data <- readRDS(here("data/pm_data.rds")) |> select(-iso)
head(pm_data)


d_ps <- d |>
    ungroup() |>
    filter(!wave %in% c("COVID-I", "COVID-II", "W9")) |>
    filter(acc >= 0) |> # some wrong data
    rename(Y1 = acc) |>
    mutate(pm = 12, covid = 0) |>
    drop_na(Y1)

d_cp <- d |>
    ungroup() |>
    filter(wave %in% c("COVID-I", "COVID-II")) |>
    rename(Y2 = acc) |>
    left_join(pm_data, c("wave", "mergeid")) |>
    mutate(Y2 = if_else(Y2 >= 0, 0, 1), covid = 1) |>
    drop_na(Y2)

d_cc <- d |>
    ungroup() |>
    filter(wave %in% c("W9")) |>
    filter(acc >= 0) |> # some wrong data
    rename(Y3 = acc) |>
    mutate(pm = 12, covid = 2) |>
    drop_na(Y3)

d_aug <-
    crossing(
        country = d |> pull(country) |> as.character() |> unique(),
        year = unique(d$year),
        sex = 0:1,
    ) |>
    mutate(
        mergeid = "K2907M2908",
        newedu = mean(d_ps$newedu),
        yedu = mean(d_ps$yedu),
        age = mean(d_ps$age),
        doctor = mean(d_ps$doctor),
        chronic = mean(d_ps$chronic),
        pm = 12,
        covid = 0,
        Y1 = NA,
    ) |>
    bind_rows(d_ps[0, ]) |>
    mutate(across(starts_with("ph006"), function(x) 0)) |>
    as.data.frame()

d_full <- bind_rows(
    d_ps, # data before covid
    d_aug, # augmented extra years for prediction up to wave 9
    d_cp, # covid waves
    d_cc, # carried over to wave 9
    .id = "d_src"
)

# country id
c_meta <- tibble(country = d$country |> unique()) |>
    arrange() |>
    mutate(
        cid = 1:n(),
        iso = countrycode::countrycode(country, "country.name", "iso3c"),
        country = if_else(iso == "CZE", "Czechia", country)
    )
write.csv(c_meta, here("data/c_meta.csv"), row.names = F)
c_meta


malta <- ne_countries(country = "Malta", type = "tiny_countries", returnclass = "sf") |>
    select(country = name_long) |>
    st_buffer(0.02)

country.sf <-
    ne_countries(
        geounit = c_meta$country,
        type = "map_units", returnclass = "sf"
    ) |>
    select(country = name) |>
    bind_rows(malta) |>
    mutate(
        iso = countrycode(country, "country.name", "iso3c")
    ) |>
    left_join(c_meta, c("iso", "country")) |>
    arrange(cid)
country.sf |> st_drop_geometry()

country.cen <- st_centroid(st_geometry(country.sf), of_largest_polygon = TRUE)

queen.nb <- country.sf |>
    st_geometry() |>
    poly2nb()

# Connect Malta and Italy
queen.nb[[20]] <- 16L
queen.nb[[16]] <- c(queen.nb[[16]], 20L)

# Connect Cyprus and Greece
queen.nb[[4]] <- 13L
queen.nb[[13]] <- c(queen.nb[[13]], 4L)

# Connect Finland and Estonia
queen.nb[[12]] <- c(queen.nb[[12]], 9L)
queen.nb[[9]] <- c(queen.nb[[9]], 12L)

# Connect Sweden and Denmark
queen.nb[[25]] <- c(queen.nb[[25]], 8L)
queen.nb[[8]] <- c(queen.nb[[8]], 25L)

queen.sf <- as(nb2lines(queen.nb, coords = country.cen, proj4string = CRS("+proj=longlat +datum=WGS84")), "sf")

nb2INLA(queen.nb, file = "eu.graph")

eu.graph <- inla.read.graph(filename = "eu.graph")

# prep for inla
d_inla <- d_full |>
    mutate(
        iso = countrycode(country, "country.name", "iso3c")
    ) |>
    select(
        starts_with("Y", F), pm,
        mergeid, iso, covid, age, newedu, year, sex,
        matches("ph006d[0-9]{1,2}")
    ) |>
    left_join(c_meta) |>
    mutate(
        iid = as.numeric(factor(mergeid)),
        # separate effects for covid changes - spatial?
        i0 = if_else(covid == 0, cid, NA_real_), # intercept
        i1 = if_else(covid == 1, cid, NA_real_), # intercept copied down
        i2 = if_else(covid == 2, cid, NA_real_), # intercept copied down
        # epidemic effects
        epi = if_else(covid == 1, cid, NA_real_),
        post = if_else(covid == 2, cid, NA_real_),
        age1 = round(age - min(age) + 1),
        age2 = age1,
        edu1 = newedu,
        edu2 = edu1,
        year = year - min(year, na.rm = T) + 1,
        t1id = if_else(covid == 0, 1L, NA_integer_),
        t2id = if_else(covid == 1, 1L, NA_integer_),
        t3id = if_else(covid == 2, 1L, NA_integer_),

        # index for diseases
        d1_0 = if_else(covid == 0, 1L, NA_integer_),
        d1_1 = if_else(covid == 1, 1L, NA_integer_),
        d1_2 = if_else(covid == 2, 1L, NA_integer_),
        d2_0 = d1_0, d2_1 = d1_1, d2_2 = d1_2,
        d3_0 = d1_0, d3_1 = d1_1, d3_2 = d1_2,
        d4_0 = d1_0, d4_1 = d1_1, d4_2 = d1_2,
        d5_0 = d1_0, d5_1 = d1_1, d5_2 = d1_2,
        d6_0 = d1_0, d6_1 = d1_1, d6_2 = d1_2,
        d10_0 = d1_0, d10_1 = d1_1, d10_2 = d1_2,
        d11_0 = d1_0, d11_1 = d1_1, d11_2 = d1_2,
        d12_0 = d1_0, d12_1 = d1_1, d12_2 = d1_2,
        d13_0 = d1_0, d13_1 = d1_1, d13_2 = d1_2,
        d14_0 = d1_0, d14_1 = d1_1, d14_2 = d1_2,
        cid1 = if_else(covid == 0, cid, NA_integer_),
        cid2 = if_else(covid == 1, cid, NA_integer_),
        cid3 = if_else(covid == 2, cid, NA_integer_),
        covid.factor = factor(covid),
        covid = covid + 1, # for `replicate`,
        iii = 1:n(),
        log_pm = log(pm)
    )

male.data <- d_inla |>
    filter(sex == 0) |>
    select(-c(1:3)) |>
    as.list()
male.data$Y <- d_inla |>
    filter(sex == 0) |>
    select(c(1:3)) |>
    as.matrix()
female.data <- d_inla |>
    filter(sex == 1) |>
    select(-c(1:3)) |>
    as.list()
female.data$Y <- d_inla |>
    filter(sex == 1) |>
    select(c(1:3)) |>
    as.matrix()

(NHV.interval <- c(1, quantile(d_inla$Y1, na.rm = T, probs = 0.975)))

# base model
hp_rw1 <- list(theta = list(prior = "pc.prec", param = c(1, 0.01)))
hp_rw2 <- list(prec = list(param = c(1, 0.00005)))
hp_iid <- list(theta = list(prior = "loggamma", param = c(1, 0.01)))
hp_fix <- list(prec = list(initial = -4, fixed = T))
hp_pcp <- list(prec = list(prior = "pcprec", param = c(1, 0.1)))
hp_bym <- list(
    prec = list(prior = "pc.prec", param = c(1, 0.01)),
    phi = list(prior = "pc", param = c(0.5, 0.5))
)
hp_beta_fix <- list(prec = list(initial = 1, fixed = TRUE))
hp_prec_fix <- list(prec = list(initial = 1, fixed = TRUE))
hp_beta_cp <- list(beta = list(fixed = FALSE)) # ~ interaction
hp_beta_cp_free <- list(beta = list(fixed = FALSE))

f0 <- Y ~
    -1 +
    f(i0) +
    f(i1, copy = "i0", precision = 1e14) +
    f(i2, copy = "i0", precision = 1e14) +

    # linear trend
    f(cid1, year) +
    f(cid2, year, copy = "cid1", precision = 1e14) +
    f(cid3, year, copy = "cid1", precision = 1e14) +

    # demographics
    f(age1, replicate = covid, model = "rw2", scale.model = TRUE) +
    f(age2, replicate = covid) +

    # covid effect separated by countries
    f(epi) +
    f(post) +

    # copy and estimate scale
    f(d1_0, ph006d1, hyper = hp_beta_fix) +
    f(d1_1, ph006d1, copy = "d1_0", hyper = hp_beta_cp_free) +
    f(d1_2, ph006d1, copy = "d1_0", hyper = hp_beta_cp_free) +
    f(d2_0, ph006d2, hyper = hp_beta_fix) +
    f(d2_1, ph006d2, copy = "d2_0", hyper = hp_beta_cp_free) +
    f(d2_2, ph006d2, copy = "d2_0", hyper = hp_beta_cp_free) +
    f(d3_0, ph006d3, hyper = hp_beta_fix) +
    f(d3_1, ph006d3, copy = "d3_0", hyper = hp_beta_cp_free) +
    f(d3_2, ph006d3, copy = "d3_0", hyper = hp_beta_cp_free) +
    f(d4_0, ph006d4, hyper = hp_beta_fix) +
    f(d4_1, ph006d4, copy = "d4_0", hyper = hp_beta_cp_free) +
    f(d4_2, ph006d4, copy = "d4_0", hyper = hp_beta_cp_free) +
    f(d5_0, ph006d5, hyper = hp_beta_fix) +
    f(d5_1, ph006d5, copy = "d5_0", hyper = hp_beta_cp_free) +
    f(d5_2, ph006d5, copy = "d5_0", hyper = hp_beta_cp_free) +
    f(d6_0, ph006d6, hyper = hp_beta_fix) +
    f(d6_1, ph006d6, copy = "d6_0", hyper = hp_beta_cp_free) +
    f(d6_2, ph006d6, copy = "d6_0", hyper = hp_beta_cp_free) +
    f(d10_0, ph006d10, hyper = hp_beta_fix) +
    f(d10_1, ph006d10, copy = "d10_0", hyper = hp_beta_cp_free) +
    f(d10_2, ph006d10, copy = "d10_0", hyper = hp_beta_cp_free) +
    f(d11_0, ph006d11, hyper = hp_beta_fix) +
    f(d11_1, ph006d11, copy = "d11_0", hyper = hp_beta_cp_free) +
    f(d11_2, ph006d11, copy = "d11_0", hyper = hp_beta_cp_free) +
    f(d12_0, ph006d13, hyper = hp_beta_fix) +
    f(d12_1, ph006d13, copy = "d12_0", hyper = hp_beta_cp_free) +
    f(d12_2, ph006d13, copy = "d12_0", hyper = hp_beta_cp_free) +
    f(d13_0, ph006d13, hyper = hp_beta_fix) +
    f(d13_1, ph006d13, copy = "d13_0", hyper = hp_beta_cp_free) +
    f(d13_2, ph006d13, copy = "d13_0", hyper = hp_beta_cp_free) +
    f(d14_0, ph006d14, hyper = hp_beta_fix) +
    f(d14_1, ph006d14, copy = "d14_0", hyper = hp_beta_cp_free) +
    f(d14_2, ph006d14, copy = "d14_0", hyper = hp_beta_cp_free) +
    f(iid) +
    f(iii)


# md.male = readRDS("/scratch/fuchs/fias/knguyen/SHARE/fitted/md.male.upd.dis.rds")
md.male <- readRDS("/scratch/fuchs/fias/knguyen/SHARE/fitted/md.male.free.rds")

md.male <- inla(f0,
    E = log_pm, data = male.data,
    family = c("poisson", "cenpoisson", "poisson"),
    control.family = list(list(), list(cenpoisson.I = NHV.interval), list()),
    control.compute = list(config = TRUE),
    control.inla = list(int.strategy = "eb"),
    control.mode = list(x = md.male$mode$x, restart = TRUE),
    control.predictor = list(link = 1),
    verbose = FALSE
)

saveRDS(md.male, "/scratch/fuchs/fias/knguyen/SHARE/fitted/md.male.free.update2.rds")
