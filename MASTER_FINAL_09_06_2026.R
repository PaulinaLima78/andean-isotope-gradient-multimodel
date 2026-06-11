# =============================================================================
# MASTER SCRIPT FINAL — Artículo isótopos MHP & CCSHP
# Autor: Paulina Lima | UCE | 2026
# Versión: 09/06/2026 — Revisión mayor EJRH
#
# CAMBIOS v09_06_2026:
#   - Sección 1: ANCOVA formal para test de diferencia entre LMWLs
#   - Figure 3A: make_biplot_panel mejorada (sin OLS por tipo, leyenda unificada,
#                "Surface water" en lugar de "Runoff", escala fill/shape separadas)
#
# OUTPUTS:
#   Table 1   → LMWL + LC-excess + ANCOVA
#   Table 2   → Medias ponderadas mensuales IAEA-INAMHI
#   Table 3   → ΔδD vs LMWL por sitio
#   Table 4   → Lapse rates δ18O y δD vs H50
#   Table 5A  → Definición de end-members e intakes
#   Table 5B  → Comparación MODEL vs MEASURED (después de sección 6.7)
#   Table 5   → Mezcla isotópica + caudal mensual
#   Figure 2  → LMM 2×2
#   Figure 3  → δ18O y δD vs H50
#   Figure 3A → δ2H vs δ18O con GMWL + LMWL por proyecto (paneles Manduriacu | Coca)
#   Figure 4  → Mezcla mensual en intakes
# =============================================================================

rm(list = ls())
gc()

# =============================================================================
# 0) CONFIGURACIÓN
# =============================================================================

BASE <- getwd()

RAW     <- file.path(BASE, "data")
PROC    <- file.path(BASE, "output")
FIG_PNG <- file.path(BASE, "output/figures/png")
FIG_TIF <- file.path(BASE, "output/figures/tiff")
TAB_CSV <- file.path(BASE, "output/tables")
TAB_SUP <- file.path(BASE, "output/supplementary")

for (d in c(PROC, FIG_PNG, FIG_TIF, TAB_CSV, TAB_SUP)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

PATH_RAIN_IAEA <- file.path(RAW, "TABLE_ISO_RAIN.csv")
PATH_RUN_O18   <- file.path(RAW, "Data_O18.csv")
PATH_RUN_D2H   <- file.path(RAW, "Data_2H.csv")
PATH_RAIN_O18  <- file.path(RAW, "DATA_O18_RAIN_HP.csv")
PATH_RAIN_D2H  <- file.path(RAW, "DATA_2H_RAIN_HP.csv")

PATH_CUENCAS_ACC <- file.path(BASE, "gis", "cuencas_acumuladas_runoff_v2_20260529_082734.shp")

cat("\n========== VERIFICACIÓN DE ARCHIVOS ==========\n")

files_check <- c(
  PATH_RAIN_IAEA,
  PATH_RUN_O18,
  PATH_RUN_D2H,
  PATH_RAIN_O18,
  PATH_RAIN_D2H,
  PATH_CUENCAS_ACC
)

all_ok <- TRUE

for (f in files_check) {
  ok <- file.exists(f)
  cat(ifelse(ok, "✅", "❌"), f, "\n")
  if (!ok) all_ok <- FALSE
}

if (!all_ok) stop("❌ Faltan archivos requeridos.")
cat("✅ Todos los archivos encontrados.\n")

# =============================================================================
# 1) PAQUETES
# =============================================================================

pkgs <- c(
  "readr", "dplyr", "tidyr", "stringr", "lubridate",
  "ggplot2", "patchwork", "ggrepel", "scales", "forcats",
  "lme4", "lmerTest", "mgcv", "nlme",
  "broom", "broom.mixed", "performance", "ggeffects", "tibble",
  "sf", "janitor"
)

miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]

if (length(miss)) install.packages(miss)

suppressPackageStartupMessages(
  invisible(lapply(pkgs, library, character.only = TRUE))
)

options(dplyr.summarise.inform = FALSE)

cat("✅ Paquetes cargados.\n")

# =============================================================================
# 2) PARÁMETROS Y HELPERS
# =============================================================================

LMWL_MHP_a <- 7.9878
LMWL_MHP_b <- 10.4301

LMWL_CCS_a <- 8.3862
LMWL_CCS_b <- 14.5124

GMWL_a <- 8
GMWL_b <- 10

MHP_SITES <- c(
  "ESMERALDAS",
  "ALLURIQUIN",
  "IZOBAMBA",
  "LA CONCORDIA",
  "QUITO-INAMHI"
)

CCSHP_SITES <- c(
  "LAGO AGRIO",
  "CUYUJA",
  "PAPALLACTA",
  "BAEZA",
  "EL CHACO"
)

# Usar "serif" evita errores de Times New Roman en Windows
FONT <- "serif"

theme_pub <- function(base_size = 12) {
  theme_bw(base_size = base_size, base_family = FONT) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold")
    )
}

save_fig <- function(p, name, w = 9, h = 6) {
  ggsave(
    file.path(FIG_PNG, paste0(name, ".png")),
    p,
    width = w,
    height = h,
    dpi = 300,
    bg = "white"
  )
  
  ggsave(
    file.path(FIG_TIF, paste0(name, ".tiff")),
    p,
    width = w,
    height = h,
    dpi = 600,
    bg = "white",
    compression = "lzw"
  )
  
  cat("✅ Figura guardada:", name, "\n")
}

wmean <- function(x, w) {
  i <- is.finite(x) & is.finite(w) & w > 0
  if (!any(i)) return(NA_real_)
  sum(x[i] * w[i]) / sum(w[i])
}

safe_mean <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else mean(x)
}

safe_sd <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <= 1) NA_real_ else sd(x)
}

safe_min <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else min(x)
}

safe_max <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) NA_real_ else max(x)
}

MONTH_MAP <- c(
  ene = "01", jan = "01",
  feb = "02",
  mar = "03",
  abr = "04", apr = "04",
  may = "05",
  jun = "06",
  jul = "07",
  ago = "08", aug = "08",
  sep = "09", sept = "09",
  oct = "10",
  nov = "11",
  dic = "12", dec = "12"
)

wide_long <- function(path, valname) {
  df <- read_csv(path, show_col_types = FALSE)
  
  id_c <- intersect(c("Id_site", "Name_site", "ID", "IdN"), names(df))
  
  if (!length(id_c)) {
    stop(paste("Sin columna ID en", basename(path)))
  }
  
  df %>%
    pivot_longer(
      -all_of(id_c),
      names_to = "mon_raw",
      values_to = valname,
      values_transform = setNames(list(as.numeric), valname)
    ) %>%
    mutate(
      mon_raw  = tolower(mon_raw),
      mon_abbr = str_extract(mon_raw, "^[a-z]+"),
      yy       = str_extract(mon_raw, "\\d{2}"),
      mm       = unname(MONTH_MAP[mon_abbr]),
      date     = suppressWarnings(ymd(paste0("20", yy, "-", mm, "-01"))),
      Id_site_u = toupper(trimws(.data[[id_c[1]]])),
      "{valname}" := suppressWarnings(as.numeric(.data[[valname]]))
    ) %>%
    filter(!is.na(date), !is.na(mm)) %>%
    select(Id_site_u, date, all_of(valname))
}

read_rain_iaea <- function(path) {
  x <- read_csv(path, show_col_types = FALSE)
  nm <- names(x)
  
  if ("Name of site" %in% nm) {
    x <- rename(x, Name_of_site = `Name of site`)
  }
  
  if (!"VO18" %in% names(x) && "O18" %in% nm) {
    x <- rename(x, VO18 = O18)
  }
  
  x %>%
    mutate(
      Name_of_site  = trimws(Name_of_site),
      VO18          = suppressWarnings(as.numeric(VO18)),
      H2            = suppressWarnings(as.numeric(H2)),
      Precipitation = suppressWarnings(as.numeric(Precipitation)),
      Altitude      = suppressWarnings(as.numeric(Altitude)),
      month         = as.integer(month)
    )
}

resid_lmwl <- function(d18O, dD, basin) {
  case_when(
    basin == "MHP"   ~ dD - (LMWL_MHP_a * d18O + LMWL_MHP_b),
    basin == "CCSHP" ~ dD - (LMWL_CCS_a * d18O + LMWL_CCS_b),
    TRUE ~ NA_real_
  )
}

assign_basin_from_id <- function(id) {
  u <- toupper(trimws(id))
  
  case_when(
    str_detect(u, "MAND|MASH|PITA|CARI|GUAY|GRAN|GUAC") ~ "MHP",
    str_detect(u, "COCA|ANTI|MALO|QUIJ|SALA|TRIB|SAN|RESR|TURB|DSCH") ~ "CCSHP",
    TRUE ~ NA_character_
  )
}

rmse_vec <- function(obs, pred) {
  i <- is.finite(obs) & is.finite(pred)
  if (!any(i)) return(NA_real_)
  sqrt(mean((obs[i] - pred[i])^2))
}

extract_model_metrics <- function(model, data, response, model_name, dataset_name) {
  pred <- tryCatch(
    predict(model, newdata = data, type = "response", allow.new.levels = TRUE),
    error = function(e) rep(NA_real_, nrow(data))
  )
  
  out <- tibble(
    Dataset = dataset_name,
    Model = model_name,
    n = nrow(data),
    AIC = suppressWarnings(tryCatch(AIC(model), error = function(e) NA_real_)),
    BIC = suppressWarnings(tryCatch(BIC(model), error = function(e) NA_real_)),
    RMSE = rmse_vec(data[[response]], pred)
  )
  
  if (inherits(model, "lm")) {
    gl <- glance(model)
    
    out <- mutate(
      out,
      R2 = gl$r.squared,
      Adj_R2 = gl$adj.r.squared,
      Dev_expl = NA_real_
    )
    
  } else if (inherits(model, "lmerMod")) {
    r2v <- tryCatch(performance::r2(model), error = function(e) NULL)
    
    out <- mutate(
      out,
      R2 = if (!is.null(r2v)) r2v$R2_marginal else NA_real_,
      Adj_R2 = if (!is.null(r2v)) r2v$R2_conditional else NA_real_,
      Dev_expl = NA_real_
    )
    
  } else if (inherits(model, "gam")) {
    sm <- summary(model)
    
    out <- mutate(
      out,
      R2 = ifelse(
        is.numeric(sm$r.sq) && is.finite(sm$r.sq),
        sm$r.sq,
        NA_real_
      ),
      Adj_R2 = NA_real_,
      Dev_expl = ifelse(
        is.numeric(sm$dev.expl) && is.finite(sm$dev.expl),
        round(sm$dev.expl, 4),
        NA_real_
      )
    )
    
  } else {
    out <- mutate(
      out,
      R2 = NA_real_,
      Adj_R2 = NA_real_,
      Dev_expl = NA_real_
    )
  }
  
  out
}

# =============================================================================
# 3) NUEVA GEOMORFOLOGÍA DESDE SHAPEFILE DE CUENCAS ACUMULADAS
# =============================================================================

cat("\n========== NUEVA GEOMORFOLOGÍA ==========\n")

cuencas_acc <- st_read(PATH_CUENCAS_ACC, quiet = FALSE)

cat("\nColumnas del shapefile:\n")
print(names(cuencas_acc))

geom_new <- cuencas_acc %>%
  st_drop_geometry() %>%
  mutate(
    Punto = toupper(trimws(Nombre_P)),
    
    Id_site_u = case_when(
      Punto == "P2A" ~ "1MAND",
      Punto == "P2"  ~ "2MAND",
      Punto == "P1"  ~ "3MASH",
      Punto == "P4"  ~ "4GUAY",
      Punto == "P5"  ~ "5PITA",
      Punto == "P15" ~ "6MICA",
      Punto == "P9"  ~ "7ANTI",
      Punto == "P8"  ~ "8GUAC",
      Punto == "P7"  ~ "9GRAN",
      Punto == "P16" ~ "10TURB",
      Punto == "P17" ~ "11TURB",
      Punto == "P18" ~ "12RESR",
      Punto == "P14" ~ "13MALO",
      Punto == "P12A" ~ "14COCA",
      Punto == "P12" ~ "15COCA",
      Punto == "P13" ~ "16SALA",
      Punto == "P11" ~ "17QUIJ",
      Punto == "P10" ~ "18TRIB",
      Punto == "P19" ~ "19RESR",
      Punto == "P3"  ~ "20SAN",
      Punto == "P20" ~ "21DSCH",
      Punto == "P6"  ~ "22CARI",
      TRUE ~ Punto
    ),
    
    Area_km2   = as.numeric(Area_km2),
    Perim_km   = as.numeric(Perim_km),
    Kc         = as.numeric(Kc),
    
    Altit_min  = as.numeric(Elev_min),
    Altit_max  = as.numeric(Elev_max),
    Altit_mean = as.numeric(Elev_mean),
    
    H50_m      = Altit_mean,
    H50_km     = H50_m / 1000,
    
    Relieve_m  = as.numeric(Relieve_m),
    Slope_mean = as.numeric(Pend_mean),
    Slope_max  = as.numeric(Pend_max),
    
    LenRio_km  = as.numeric(LenRio_km),
    Dren_dens  = as.numeric(Dren_dens),
    
    HI = ifelse(
      is.finite(Altit_min) &
        is.finite(Altit_max) &
        is.finite(H50_m) &
        Altit_max > Altit_min,
      (H50_m - Altit_min) / (Altit_max - Altit_min),
      NA_real_
    ),
    
    Project = case_when(
      Id_site_u == "6MICA" ~ "Interbasin transfer / Quito water supply",
      str_detect(Id_site_u, "MAND|MASH|PITA|CARI|GUAY|GRAN|GUAC") ~ "Manduriacu (a)",
      str_detect(Id_site_u, "COCA|ANTI|MALO|QUIJ|SALA|TRIB|SAN|RESR|TURB|DSCH") ~ "Coca Codo Sinclair (b)",
      TRUE ~ NA_character_
    ),
    
    Basin = case_when(
      Project == "Manduriacu (a)" ~ "MHP",
      Project == "Coca Codo Sinclair (b)" ~ "CCSHP",
      Project == "Interbasin transfer / Quito water supply" ~ "Interbasin transfer",
      TRUE ~ NA_character_
    ),
    
    Drainage_System = Project,
    
    Source_zone = case_when(
      H50_m >= 3000 ~ "High mountain / headwater",
      H50_m < 3000  ~ "Low-elevation / local",
      TRUE ~ NA_character_
    )
  ) %>%
  select(
    Id_site_u,
    Nombre_P,
    Basin,
    Project,
    Drainage_System,
    Area_km2,
    Perim_km,
    Kc,
    Altit_min,
    Altit_mean,
    Altit_max,
    H50_m,
    H50_km,
    Relieve_m,
    HI,
    Slope_mean,
    Slope_max,
    LenRio_km,
    Dren_dens,
    Source_zone
  )

write_csv(
  geom_new,
  file.path(PROC, "Geomorf_actualizada_cuencas_acumuladas_20260529.csv")
)

cat("\n✅ geom_new creado correctamente desde el shapefile actualizado.\n")

print(
  geom_new %>%
    select(
      Id_site_u,
      Nombre_P,
      Basin,
      Project,
      Area_km2,
      H50_m,
      HI,
      Source_zone
    )
)

# =============================================================================
# SECCIÓN 1 — TABLE 1
# =============================================================================

cat("\n========== SECCIÓN 1: TABLE 1 — LMWL + ANCOVA ==========\n")

rain_iaea <- read_rain_iaea(PATH_RAIN_IAEA) %>%
  filter(Name_of_site %in% c(MHP_SITES, CCSHP_SITES)) %>%
  mutate(
    Basin = case_when(
      Name_of_site %in% MHP_SITES ~ "MHP",
      Name_of_site %in% CCSHP_SITES ~ "CCSHP",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(VO18),
    !is.na(H2),
    month %in% 1:12,
    Basin %in% c("MHP", "CCSHP")
  ) %>%
  mutate(
    LC_excess = resid_lmwl(VO18, H2, Basin)
  )

# ---- ANCOVA: test formal de diferencia entre LMWLs ----
# H0: las dos LMWLs tienen el mismo intercepto Y la misma pendiente
# Modelo: δ2H ~ δ18O * Basin (interacción = test de pendiente diferente)
lmwl_test_data <- rain_iaea %>%
  filter(
    Basin %in% c("MHP", "CCSHP"),
    is.finite(VO18),
    is.finite(H2)
  ) %>%
  mutate(
    Basin = factor(Basin, levels = c("CCSHP", "MHP"))
  )

ancova_lmwl    <- lm(H2 ~ VO18 * Basin, data = lmwl_test_data)
anova_lmwl     <- anova(ancova_lmwl)

ancova_summary <- anova_lmwl %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Term") %>%
  rename(
    Df      = Df,
    Sum_Sq  = `Sum Sq`,
    Mean_Sq = `Mean Sq`,
    F_value = `F value`,
    p_value = `Pr(>F)`
  )

format_p <- function(p) {
  case_when(
    is.na(p)   ~ NA_character_,
    p < 0.001  ~ "<0.001",
    TRUE       ~ sprintf("%.3f", p)
  )
}

p_intercept <- ancova_summary %>% filter(Term == "Basin")      %>% pull(p_value)
F_intercept <- ancova_summary %>% filter(Term == "Basin")      %>% pull(F_value)
p_slope     <- ancova_summary %>% filter(Term == "VO18:Basin") %>% pull(p_value)
F_slope     <- ancova_summary %>% filter(Term == "VO18:Basin") %>% pull(F_value)

cat("\n--- ANCOVA LMWL: test de diferencia entre cuencas ---\n")
cat(sprintf("  Intercepto (Basin):      F = %.2f, p = %s\n",
            F_intercept, format_p(p_intercept)))
cat(sprintf("  Pendiente (VO18:Basin):  F = %.2f, p = %s\n",
            F_slope, format_p(p_slope)))

# ---- Table 1 integrada ----
table1 <- rain_iaea %>%
  group_by(Basin) %>%
  summarise(
    Slope_a        = ifelse(first(Basin) == "MHP", LMWL_MHP_a, LMWL_CCS_a),
    Intercept_b    = ifelse(first(Basin) == "MHP", LMWL_MHP_b, LMWL_CCS_b),
    mean_LC_excess = round(mean(LC_excess, na.rm = TRUE), 3),
    sd_LC_excess   = round(sd(LC_excess,   na.rm = TRUE), 2),
    n              = n(),
    .groups        = "drop"
  ) %>%
  mutate(
    LMWL_eq = sprintf("\u03b4\u00b2H = %.4f\u00b7\u03b4\u00b9\u2078O + %.4f",
                      Slope_a, Intercept_b),
    ANCOVA_intercept_test = sprintf("F = %.2f, p = %s",
                                    F_intercept, format_p(p_intercept)),
    ANCOVA_slope_test     = sprintf("F = %.2f, p = %s",
                                    F_slope,     format_p(p_slope))
  )

print(table1)

write_csv(table1,
          file.path(TAB_CSV, "Table1_LMWL_LCexcess_ANCOVA.csv"))
write_csv(ancova_summary,
          file.path(TAB_CSV, "Table1_LMWL_ANCOVA_details.csv"))
write_csv(broom::tidy(ancova_lmwl),
          file.path(TAB_CSV, "Table1_LMWL_ANCOVA_coefficients.csv"))

cat("✅ Table 1 guardada con LMWL, LC-excess y ANCOVA.\n")

# =============================================================================
# SECCIÓN 2 — TABLE 2
# =============================================================================

cat("\n========== SECCIÓN 2: TABLE 2 — Medias ponderadas MHP + CCSHP ==========\n")

make_table2 <- function(basin_name) {
  rain_iaea %>%
    filter(Basin == basin_name, is.finite(Precipitation), Precipitation > 0) %>%
    group_by(month) %>%
    summarise(
      Basin = basin_name,
      d18O_w = round(wmean(VO18, Precipitation), 2),
      dD_w = round(wmean(H2, Precipitation), 2),
      Precip_mm = round(mean(Precipitation, na.rm = TRUE), 1),
      n_d18O = sum(is.finite(VO18)),
      n_dD = sum(is.finite(H2)),
      d_excess = round(wmean(H2, Precipitation) - 8 * wmean(VO18, Precipitation), 2),
      .groups = "drop"
    ) %>%
    arrange(month) %>%
    mutate(Month = month.abb[month]) %>%
    select(Basin, Month, d18O_w, dD_w, Precip_mm, n_d18O, n_dD, d_excess)
}

table2 <- bind_rows(
  make_table2("MHP"),
  make_table2("CCSHP")
)

cat("\n--- MHP ---\n")
print(table2 %>% filter(Basin == "MHP"))

cat("\n--- CCSHP ---\n")
print(table2 %>% filter(Basin == "CCSHP"))

write_csv(
  table2,
  file.path(TAB_CSV, "Table2_monthly_weighted_MHP_CCSHP.csv")
)

cat("✅ Table 2 guardada.\n")

# =============================================================================
# SECCIÓN 3 — FIGURE 2: LMM 2×2
# =============================================================================

cat("\n========== SECCIÓN 3: FIGURE 2 — LMM ==========\n")

dat_lmm <- rain_iaea %>%
  filter(
    !is.na(VO18),
    !is.na(Altitude),
    month %in% 1:12
  ) %>%
  mutate(
    Name_of_site = as.factor(Name_of_site)
  ) %>%
  group_by(Name_of_site) %>%
  mutate(
    Altitude = median(Altitude, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  mutate(
    month_f = factor(month, levels = 1:12, labels = month.abb),
    Alt_c = as.numeric(scale(Altitude))
  )

m_lmm_int <- lmer(
  VO18 ~ Alt_c + month_f + (1 | Name_of_site),
  data = dat_lmm
)

r2_vals <- performance::r2(m_lmm_int)

cat(
  sprintf(
    "R² marginal = %.3f | R² conditional = %.3f\n",
    r2_vals$R2_marginal,
    r2_vals$R2_conditional
  )
)

fixef_tbl <- broom.mixed::tidy(
  m_lmm_int,
  effects = "fixed",
  conf.int = TRUE,
  conf.method = "Wald"
) %>%
  mutate(
    p.value = 2 * (1 - pnorm(abs(statistic)))
  )

forest_df <- fixef_tbl %>%
  filter(term != "(Intercept)") %>%
  mutate(
    term_label = case_when(
      term == "Alt_c" ~ "Altitude (z-score)",
      grepl("^month_f", term) ~ gsub("month_f", "", term),
      TRUE ~ term
    ),
    sig = if_else(p.value < 0.05, "p < 0.05", "n.s."),
    term_label = factor(term_label, levels = c("Altitude (z-score)", month.abb))
  ) %>%
  arrange(term_label)

p_forest <- ggplot(
  forest_df,
  aes(
    estimate,
    term_label,
    xmin = conf.low,
    xmax = conf.high,
    colour = sig
  )
) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(size = 0.4) +
  scale_colour_manual(values = c("p < 0.05" = "#009E73", "n.s." = "#999999")) +
  labs(
    title = "A) Fixed effects (95% CI)",
    x = expression("Change in " * delta^{18} * "O (‰)"),
    y = NULL
  ) +
  theme_pub() +
  theme(legend.position = "bottom")

eff_alt <- ggpredict(m_lmm_int, terms = "Alt_c [all]")

p_alt <- ggplot(eff_alt, aes(x, predicted)) +
  geom_line() +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(
    title = "B) Altitude effect",
    x = "Altitude (z-score)",
    y = expression(delta^{18} * "O (‰)")
  ) +
  theme_pub()

eff_month <- ggpredict(m_lmm_int, terms = "month_f")

p_month <- ggplot(eff_month, aes(x, predicted, group = 1)) +
  geom_line() +
  geom_point(size = 1.8) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), alpha = 0.2) +
  labs(
    title = "C) Month effect",
    x = "Month",
    y = expression(delta^{18} * "O (‰)")
  ) +
  theme_pub()

re_list <- ranef(m_lmm_int, condVar = TRUE)
ri_mat <- re_list$Name_of_site
postVar <- attr(ri_mat, "postVar")

se_vec <- sqrt(
  sapply(
    1:dim(postVar)[3],
    function(i) postVar[1, 1, i]
  )
)

df_ri <- tibble(
  Site = rownames(ri_mat),
  intercept = ri_mat[, "(Intercept)"],
  se = se_vec,
  lo = intercept - 1.96 * se,
  hi = intercept + 1.96 * se
) %>%
  arrange(intercept) %>%
  mutate(Site = fct_reorder(Site, intercept))

p_cater <- ggplot(df_ri, aes(intercept, Site)) +
  geom_vline(xintercept = 0, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.4) +
  labs(
    title = "D) Random intercepts by site (95% CI)",
    x = expression("Intercept deviation (‰)"),
    y = "Site"
  ) +
  theme_pub()

fig2 <- (p_forest | p_alt) / (p_month | p_cater) +
  plot_annotation(
    title = expression("Mixed-effects model: spatio-temporal variability of " * delta^{18} * "O"),
    subtitle = sprintf(
      "R² marginal = %.3f | R² conditional = %.3f",
      r2_vals$R2_marginal,
      r2_vals$R2_conditional
    ),
    theme = theme(
      plot.title = element_text(face = "bold", family = FONT)
    )
  )

print(fig2)

save_fig(
  fig2,
  "Figure2_LMM_d18O_2x2",
  w = 14,
  h = 10
)

cat("✅ Figure 2 guardada.\n")
# =============================================================================
# SECCIÓN 4 — TABLE 3: ΔδD vs LMWL
# =============================================================================

cat("\n========== SECCIÓN 4: TABLE 3 ==========\n")

geom_t3 <- geom_new %>%
  mutate(
    Id_site_u = toupper(trimws(Id_site_u)),
    Basin = case_when(
      Basin %in% c("MHP", "CCSHP") ~ Basin,
      str_detect(toupper(Drainage_System), "MHP|MAND") ~ "MHP",
      str_detect(toupper(Drainage_System), "CCSHP|COCA") ~ "CCSHP",
      TRUE ~ assign_basin_from_id(Id_site_u)
    )
  ) %>%
  filter(Basin %in% c("MHP", "CCSHP")) %>%
  select(Id_site_u, Basin) %>%
  distinct()

rain18_l <- wide_long(PATH_RAIN_O18, "d18O")
rainD_l  <- wide_long(PATH_RAIN_D2H, "dD")

rain_field <- full_join(
  rain18_l, rainD_l,
  by = c("Id_site_u", "date"),
  relationship = "many-to-many"
) %>%
  left_join(geom_t3, by = "Id_site_u") %>%
  mutate(
    Basin = ifelse(is.na(Basin), assign_basin_from_id(Id_site_u), Basin),
    SampleType = "Rainfall",
    deltaD_resid = resid_lmwl(d18O, dD, Basin)
  ) %>%
  filter(Basin %in% c("MHP", "CCSHP"), is.finite(d18O), is.finite(dD))

run18_t3 <- wide_long(PATH_RUN_O18, "d18O")
runD_t3  <- wide_long(PATH_RUN_D2H, "dD")

runoff_t3 <- full_join(
  run18_t3, runD_t3,
  by = c("Id_site_u", "date"),
  relationship = "many-to-many"
) %>%
  left_join(geom_t3, by = "Id_site_u") %>%
  mutate(
    Basin = ifelse(is.na(Basin), assign_basin_from_id(Id_site_u), Basin),
    SampleType = "Runoff",
    deltaD_resid = resid_lmwl(d18O, dD, Basin),
    slope = ifelse(Basin == "MHP", LMWL_MHP_a, LMWL_CCS_a),
    dist_LMWL = abs(deltaD_resid) / sqrt(1 + slope^2),
    class = case_when(
      abs(deltaD_resid) < 1.0 ~ "on",
      deltaD_resid < 0 ~ "below",
      TRUE ~ "above"
    )
  ) %>%
  filter(Basin %in% c("MHP", "CCSHP"), is.finite(d18O), is.finite(dD))

summ_t3 <- function(df) {
  df %>%
    group_by(SampleType, Basin, Id_site_u) %>%
    summarise(
      n = n(),
      d18O_mean = round(safe_mean(d18O), 2),
      d18O_sd   = round(safe_sd(d18O), 2),
      dD_mean   = round(safe_mean(dD), 2),
      dD_sd     = round(safe_sd(dD), 2),
      Delta_dD_mean = round(safe_mean(deltaD_resid), 2),
      Delta_dD_sd   = round(safe_sd(deltaD_resid), 2),
      .groups = "drop"
    )
}

table3 <- bind_rows(
  summ_t3(rain_field),
  summ_t3(runoff_t3)
) %>%
  arrange(SampleType, Basin, Id_site_u)

table3_class <- runoff_t3 %>%
  group_by(Basin, Id_site_u) %>%
  summarise(
    n = n(),
    n_on = sum(class == "on"),
    n_below = sum(class == "below"),
    n_above = sum(class == "above"),
    pct_on = round(100 * n_on / n, 1),
    Delta_dD_med = round(median(deltaD_resid, na.rm = TRUE), 2),
    Med_dist_LMWL = round(median(dist_LMWL, na.rm = TRUE), 3),
    .groups = "drop"
  )

print(table3)
print(table3_class)

write_csv(table3, file.path(TAB_CSV, "Table3_deltaD_by_site.csv"))
write_csv(table3_class, file.path(TAB_CSV, "Table3_runoff_LMWL_class.csv"))

cat("✅ Table 3 guardada.\n")

# =============================================================================
# FIGURE 3A — δ2H vs δ18O: GMWL + basin-specific LMWL
# Reviewer request: show rainfall and surface-water samples by project
# =============================================================================

cat("\n========== FIGURE 3A: δ2H vs δ18O — GMWL + LMWL ==========\n")

# -----------------------------------------------------------------------------
# DIAGNÓSTICO DE OUTLIERS EN LLUVIA DE CAMPO
# -----------------------------------------------------------------------------

rain_field_diag <- rain_field %>%
  group_by(Basin) %>%
  mutate(
    resid_lmwl_val = resid_lmwl(d18O, dD, Basin),
    q1  = quantile(resid_lmwl_val, 0.25, na.rm = TRUE),
    q3  = quantile(resid_lmwl_val, 0.75, na.rm = TRUE),
    iqr = q3 - q1,
    outlier_flag = abs(resid_lmwl_val - median(resid_lmwl_val, na.rm = TRUE)) >
      3 * iqr
  ) %>%
  ungroup()

rain_outliers <- rain_field_diag %>%
  filter(outlier_flag) %>%
  select(Basin, Id_site_u, date, d18O, dD, resid_lmwl_val)

print(rain_outliers)

write_csv(
  rain_outliers,
  file.path(TAB_CSV, "Figure3A_rainfall_outliers_flagged.csv")
)

# -----------------------------------------------------------------------------
# MAPA DE ID A PUNTO
# -----------------------------------------------------------------------------

site_label_map <- tibble::tribble(
  ~Id_site_u, ~Punto,
  "1MAND",  "P2a",
  "2MAND",  "P2",
  "3MASH",  "P1",
  "4GUAY",  "P4",
  "5PITA",  "P5",
  "6MICA",  "P15",
  "7ANTI",  "P9",
  "8GUAC",  "P8",
  "9GRAN",  "P7",
  "10TURB", "P16",
  "11TURB", "P17",
  "12RESR", "P18",
  "13MALO", "P14",
  "14COCA", "P12a",
  "15COCA", "P12",
  "16SALA", "P13",
  "17QUIJ", "P11",
  "18TRIB", "P10",
  "19RESR", "P19",
  "20SAN",  "P3",
  "21DSCH", "P20",
  "22CARI", "P6"
)

# -----------------------------------------------------------------------------
# DATASET PARA FIGURA
# -----------------------------------------------------------------------------

iso_biplot_data <- bind_rows(
  rain_field_diag %>%
    filter(!outlier_flag) %>%
    select(Basin, Id_site_u, date, d18O, dD, SampleType),
  runoff_t3 %>%
    select(Basin, Id_site_u, date, d18O, dD, SampleType)
) %>%
  filter(
    Basin %in% c("MHP", "CCSHP"),
    is.finite(d18O),
    is.finite(dD)
  ) %>%
  left_join(site_label_map, by = "Id_site_u") %>%
  mutate(
    Project = case_when(
      Basin == "MHP"   ~ "Manduriacu (a)",
      Basin == "CCSHP" ~ "Coca Codo Sinclair (b)"
    ),
    Project = factor(
      Project,
      levels = c("Manduriacu (a)", "Coca Codo Sinclair (b)")
    ),
    SampleType = factor(
      SampleType,
      levels = c("Rainfall", "Runoff"),
      labels = c("Rainfall", "Surface water")
    ),
    Site_label = paste0(Punto, " (", SampleType, ")")
  ) %>%
  filter(!is.na(Punto))

iso_biplot_all <- bind_rows(
  rain_field %>%
    select(Basin, Id_site_u, date, d18O, dD, SampleType),
  runoff_t3 %>%
    select(Basin, Id_site_u, date, d18O, dD, SampleType)
) %>%
  filter(
    Basin %in% c("MHP", "CCSHP"),
    is.finite(d18O),
    is.finite(dD)
  ) %>%
  left_join(site_label_map, by = "Id_site_u") %>%
  mutate(
    Project = case_when(
      Basin == "MHP"   ~ "Manduriacu (a)",
      Basin == "CCSHP" ~ "Coca Codo Sinclair (b)"
    ),
    Project = factor(
      Project,
      levels = c("Manduriacu (a)", "Coca Codo Sinclair (b)")
    ),
    SampleType = factor(
      SampleType,
      levels = c("Rainfall", "Runoff"),
      labels = c("Rainfall", "Surface water")
    )
  ) %>%
  filter(!is.na(Punto))

# =============================================================================
# PALETA Y ROLES HIDROLÓGICOS PARA FIGURE 3A
# Consistente con Figure 4:
# Headwater = orange | Local = blue | Intake = teal | Downstream = dark grey
# =============================================================================

COL_HEAD  <- "#F4A62A"
COL_LOCAL <- "#2E86DE"
COL_INTAK <- "#2A9D8F"
COL_DOWN  <- "grey40"
COL_OTHER <- "grey70"

assign_plot_role <- function(Basin, Punto) {
  case_when(
    Basin == "MHP"   & Punto %in% c("P5", "P6")    ~ "Headwater reference",
    Basin == "MHP"   & Punto == "P1"               ~ "Local reference",
    Basin == "MHP"   & Punto %in% c("P2", "P2a")   ~ "Intake",
    
    Basin == "CCSHP" & Punto == "P9"               ~ "Headwater reference",
    Basin == "CCSHP" & Punto == "P14"              ~ "Local reference",
    Basin == "CCSHP" & Punto %in% c("P12", "P12a") ~ "Intake",
    Basin == "CCSHP" & Punto %in% c("P16", "P17", "P18", "P19", "P20") ~
      "Downstream / regulated reach",
    
    TRUE ~ "Other tributary / monitoring site"
  )
}

role_levels <- c(
  "Headwater reference",
  "Local reference",
  "Intake",
  "Downstream / regulated reach",
  "Other tributary / monitoring site"
)

role_colors <- c(
  "Headwater reference" = COL_HEAD,
  "Local reference" = COL_LOCAL,
  "Intake" = COL_INTAK,
  "Downstream / regulated reach" = COL_DOWN,
  "Other tributary / monitoring site" = COL_OTHER
)

iso_biplot_data <- iso_biplot_data %>%
  mutate(
    Plot_role = assign_plot_role(Basin, Punto),
    Plot_role = factor(Plot_role, levels = role_levels)
  )

iso_biplot_all <- iso_biplot_all %>%
  mutate(
    Plot_role = assign_plot_role(Basin, Punto),
    Plot_role = factor(Plot_role, levels = role_levels)
  )

cat("\n--- Check Plot_role Figure 3A ---\n")
print(table(iso_biplot_data$Plot_role, useNA = "ifany"))

# =============================================================================
# FUNCIÓN FIGURE 3A — SOLO GMWL + LMWL
# =============================================================================

make_biplot_panel <- function(basin_code, basin_title,
                              lmwl_a, lmwl_b,
                              data_all,
                              lmwl_label_x = NULL) {
  
  d <- data_all %>%
    filter(
      Basin == basin_code,
      is.finite(d18O),
      is.finite(dD),
      !is.na(Punto)
    ) %>%
    mutate(
      SampleType = factor(SampleType, levels = c("Rainfall", "Surface water")),
      Plot_role = factor(Plot_role, levels = role_levels)
    )
  
  if (nrow(d) == 0) stop(paste("No hay datos para", basin_code))
  
  x_range <- range(d$d18O, na.rm = TRUE)
  y_range <- range(d$dD, na.rm = TRUE)
  
  x_pad <- diff(x_range) * 0.08
  y_pad <- diff(y_range) * 0.08
  
  x_lo <- x_range[1] - x_pad
  x_hi <- x_range[2] + x_pad
  
  line_df <- bind_rows(
    tibble(
      d18O = c(x_lo, x_hi),
      dD   = c(GMWL_a * x_lo + GMWL_b, GMWL_a * x_hi + GMWL_b),
      Line = "GMWL"
    ),
    tibble(
      d18O = c(x_lo, x_hi),
      dD   = c(lmwl_a * x_lo + lmwl_b, lmwl_a * x_hi + lmwl_b),
      Line = "Basin-specific LMWL"
    )
  ) %>%
    mutate(Line = factor(Line, levels = c("GMWL", "Basin-specific LMWL")))
  
  lmwl_lx <- if (is.null(lmwl_label_x)) x_lo + diff(x_range) * 0.05 else lmwl_label_x
  lmwl_ly <- lmwl_a * lmwl_lx + lmwl_b + diff(y_range) * 0.06
  
  lmwl_eq <- sprintf(
    "\u03b4\u00b2H = %.4f\u00b7\u03b4\u00b9\u2078O %+.4f",
    lmwl_a, lmwl_b
  )
  
  label_df <- d %>%
    group_by(Punto, Plot_role) %>%
    summarise(
      d18O = mean(d18O, na.rm = TRUE),
      dD   = mean(dD, na.rm = TRUE),
      .groups = "drop"
    )
  
  ggplot(d, aes(x = d18O, y = dD)) +
    geom_line(
      data = line_df,
      aes(x = d18O, y = dD, linetype = Line, colour = Line),
      linewidth = 1.0,
      inherit.aes = FALSE
    ) +
    annotate(
      "text",
      x = lmwl_lx,
      y = lmwl_ly,
      label = lmwl_eq,
      colour = "#185FA5",
      size = 3,
      hjust = 0,
      fontface = "italic",
      family = FONT
    ) +
    geom_point(
      aes(fill = Plot_role, shape = SampleType),
      colour = "grey20",
      size = 2.8,
      alpha = 0.90,
      stroke = 0.35
    ) +
    ggrepel::geom_text_repel(
      data = label_df,
      aes(x = d18O, y = dD, label = Punto),
      inherit.aes = FALSE,
      size = 3.2,
      family = FONT,
      colour = "black",
      fontface = "bold",
      
      nudge_y = 4,          # mueve etiquetas hacia arriba
      direction = "y",      # solo movimiento vertical
      
      box.padding = 0.6,
      point.padding = 0.5,
      min.segment.length = 0,
      
      segment.color = "grey40",
      segment.size = 0.3,
      
      max.overlaps = Inf,
      seed = 123
    ) +
    scale_linetype_manual(
      name = "Reference line",
      values = c("GMWL" = "dashed", "Basin-specific LMWL" = "solid")
    ) +
    scale_colour_manual(
      name = "Reference line",
      values = c("GMWL" = "black", "Basin-specific LMWL" = "#185FA5")
    ) +
    scale_shape_manual(
      name = "Sample type",
      values = c("Rainfall" = 24, "Surface water" = 21)
    ) +
    scale_fill_manual(
      name = "Hydrologic role",
      values = role_colors,
      drop = FALSE
    ) +
    coord_cartesian(
      xlim = c(x_lo, x_hi),
      ylim = c(y_range[1] - y_pad, y_range[2] + y_pad)
    ) +
    labs(
      title = basin_title,
      x = expression(delta^{18}*"O (\u2030)"),
      y = expression(delta^2*"H (\u2030)")
    ) +
    theme_pub(base_size = 11) +
    theme(
      legend.position  = "right",
      legend.key.size  = unit(0.45, "cm"),
      legend.text      = element_text(size = 8),
      legend.title     = element_text(face = "bold", size = 9),
      plot.title       = element_text(face = "bold", hjust = 0.5, size = 12),
      
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      
      panel.border = element_rect(
        colour = "black",
        fill = NA,
        linewidth = 0.5
      )
    ) +
    guides(
      fill = guide_legend(
        order = 2,
        override.aes = list(shape = 21, size = 3.5)
      ),
      shape = guide_legend(
        order = 3,
        override.aes = list(size = 3.5, fill = "grey70")
      ),
      colour = guide_legend(order = 1),
      linetype = guide_legend(order = 1)
    )
}

# =============================================================================
# PANEL A — MANDURIACU
# =============================================================================

fig3A_MHP <- make_biplot_panel(
  basin_code  = "MHP",
  basin_title = "(a) Manduriacu",
  lmwl_a      = LMWL_MHP_a,
  lmwl_b      = LMWL_MHP_b,
  data_all    = iso_biplot_data
)

save_fig(
  fig3A_MHP,
  "Figure3A_MHP_d2H_d18O_GMWL_LMWL",
  w = 9,
  h = 7
)

# =============================================================================
# PANEL B — COCA CODO SINCLAIR
# =============================================================================

fig3A_CCSHP <- make_biplot_panel(
  basin_code  = "CCSHP",
  basin_title = "(b) Coca Codo Sinclair",
  lmwl_a      = LMWL_CCS_a,
  lmwl_b      = LMWL_CCS_b,
  data_all    = iso_biplot_data
)

save_fig(
  fig3A_CCSHP,
  "Figure3A_CCSHP_d2H_d18O_GMWL_LMWL",
  w = 9,
  h = 7
)

# =============================================================================
# FIGURA COMBINADA — SIN TÍTULO GENERAL
# =============================================================================

fig3A_combined <- fig3A_MHP + fig3A_CCSHP

save_fig(
  fig3A_combined,
  "Figure3A_d2H_d18O_GMWL_LMWL_by_project",
  w = 17,
  h = 8
)

write_csv(
  iso_biplot_data,
  file.path(TAB_CSV, "Figure3A_biplot_data_GMWL_LMWL_clean.csv")
)

write_csv(
  iso_biplot_all,
  file.path(TAB_CSV, "Figure3A_biplot_data_GMWL_LMWL_all.csv")
)

cat("✅ Figure 3A guardada sin título general, con GMWL + LMWL y roles hidrológicos claros.\n")
# =============================================================================
# CHECK SAMPLE COUNTS — RUNOFF + RAINFALL
# =============================================================================

cat("\n========== SAMPLE COUNT SUMMARY ==========\n")

runoff_samples <- runoff_t3 %>%
  filter(is.finite(d18O), is.finite(dD))

rain_samples <- rain_field %>%
  filter(is.finite(d18O), is.finite(dD))

n_runoff_total <- nrow(runoff_samples)
n_runoff_sites <- n_distinct(runoff_samples$Id_site_u)

n_rain_total <- nrow(rain_samples)
n_rain_sites <- n_distinct(rain_samples$Id_site_u)

runoff_dates <- range(runoff_samples$date, na.rm = TRUE)
rain_dates   <- range(rain_samples$date, na.rm = TRUE)

n_months_runoff <- length(unique(format(runoff_samples$date, "%Y-%m")))
n_months_rain   <- length(unique(format(rain_samples$date, "%Y-%m")))

cat("\n--- RUNOFF ---\n")
cat("Samples:", n_runoff_total, "\n")
cat("Sites:", n_runoff_sites, "\n")
cat("Period:", as.character(runoff_dates[1]), "to", as.character(runoff_dates[2]), "\n")
cat("Months covered:", n_months_runoff, "\n")

cat("\n--- RAINFALL (FIELD) ---\n")
cat("Samples:", n_rain_total, "\n")
cat("Sites:", n_rain_sites, "\n")
cat("Period:", as.character(rain_dates[1]), "to", as.character(rain_dates[2]), "\n")
cat("Months covered:", n_months_rain, "\n")

cat("\n--- TOTAL ---\n")
cat("Total samples:", n_runoff_total + n_rain_total, "\n")

sample_summary <- tibble(
  Dataset = c("Runoff", "Rainfall"),
  Samples = c(n_runoff_total, n_rain_total),
  Sites   = c(n_runoff_sites, n_rain_sites),
  Months  = c(n_months_runoff, n_months_rain),
  Start   = c(runoff_dates[1], rain_dates[1]),
  End     = c(runoff_dates[2], rain_dates[2])
)

write_csv(sample_summary, file.path(TAB_CSV, "Table_sample_summary.csv"))
write_csv(runoff_samples, file.path(TAB_CSV, "Runoff_samples_used_in_analysis.csv"))
write_csv(rain_samples, file.path(TAB_CSV, "Rainfall_samples_used_in_analysis.csv"))

cat("\n✅ Sample summary table saved.\n")

# =============================================================================
# SECCIÓN 5 — TABLE 4 + FIGURE 3
# =============================================================================

cat("\n========== SECCIÓN 5: TABLE 4 + FIGURE 3 ==========\n")

geom <- geom_new %>%
  mutate(
    Id_site_u = toupper(trimws(Id_site_u)),
    H50_km = H50_m / 1000,
    Project = case_when(
      Project %in% c("Manduriacu (a)", "Coca Codo Sinclair (b)") ~ Project,
      Basin == "MHP" ~ "Manduriacu (a)",
      Basin == "CCSHP" ~ "Coca Codo Sinclair (b)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    is.finite(H50_m),
    !is.na(Project),
    Project %in% c("Manduriacu (a)", "Coca Codo Sinclair (b)")
  )

run18_l <- wide_long(PATH_RUN_O18, "d18O") %>%
  group_by(Id_site_u, date) %>%
  summarise(d18O = mean(d18O, na.rm = TRUE), .groups = "drop")

runD_l <- wide_long(PATH_RUN_D2H, "dD") %>%
  group_by(Id_site_u, date) %>%
  summarise(dD = mean(dD, na.rm = TRUE), .groups = "drop")

runoff <- full_join(
  run18_l, runD_l,
  by = c("Id_site_u", "date"),
  relationship = "many-to-many"
) %>%
  left_join(
    geom %>%
      select(
        Id_site_u, H50_m, H50_km, Project,
        Area_km2, Altit_min, Altit_mean, Altit_max,
        HI, Slope_mean, Source_zone
      ),
    by = "Id_site_u"
  ) %>%
  filter(is.finite(H50_m), !is.na(Project))

# Sitios excluidos del gradiente principal de runoff:
# 6MICA = transferencia intercuencas Quito
# 18TRIB, 19RESR, 17TURB, 20TURB, 16DSCH = infraestructura / reservorio / descarga
excl_ids <- c("6MICA", "18TRIB", "19RESR", "17TURB", "20TURB", "16DSCH")

runoff_clean <- runoff %>%
  filter(!Id_site_u %in% excl_ids)

# -----------------------------------------------------------------------------
# Usar medias por sitio para evitar pseudorreplicación temporal
# -----------------------------------------------------------------------------

runoff_site_mean <- runoff_clean %>%
  group_by(
    Project, Id_site_u,
    H50_m, H50_km,
    Area_km2, Altit_min, Altit_mean, Altit_max,
    HI, Slope_mean, Source_zone
  ) %>%
  summarise(
    n_months = n_distinct(date),
    d18O = mean(d18O, na.rm = TRUE),
    dD   = mean(dD, na.rm = TRUE),
    d18O_sd = sd(d18O, na.rm = TRUE),
    dD_sd   = sd(dD, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(is.finite(d18O), is.finite(dD))

write_csv(
  runoff_site_mean,
  file.path(TAB_CSV, "Table4_site_mean_runoff_isotopes_H50.csv")
)

fit_lapse <- function(df, yvar, proj) {
  d <- df %>%
    filter(
      Project == proj,
      is.finite(.data[[yvar]]),
      is.finite(H50_km)
    )
  
  if (nrow(d) < 3) return(NULL)
  
  fit <- lm(reformulate("H50_km", response = yvar), data = d)
  td <- broom::tidy(fit, conf.int = TRUE)
  gl <- broom::glance(fit)
  
  b  <- td[td$term == "H50_km", ]
  ic <- td[td$term == "(Intercept)", ]
  
  sw <- tryCatch(
    shapiro.test(residuals(fit))$p.value,
    error = function(e) NA_real_
  )
  
  tibble(
    Project = proj,
    Isotope = yvar,
    n_sites = nrow(d),
    R2 = round(gl$r.squared, 3),
    p_model = signif(gl$p.value, 3),
    Slope_per_km = round(b$estimate, 3),
    CI_low = round(b$conf.low, 3),
    CI_high = round(b$conf.high, 3),
    Intercept = round(ic$estimate, 3),
    Pearson_r = round(cor(d$H50_km, d[[yvar]], use = "complete.obs"), 3),
    Shapiro_p = round(sw, 3),
    Equation = sprintf(
      "%s = %.3f %+.3f·H50_km",
      yvar, ic$estimate, b$estimate
    )
  )
}

projs <- c("Manduriacu (a)", "Coca Codo Sinclair (b)")

table4 <- bind_rows(lapply(projs, function(p) {
  bind_rows(
    fit_lapse(runoff_site_mean, "d18O", p),
    fit_lapse(runoff_site_mean, "dD", p)
  )
}))

print(table4)

write_csv(
  table4,
  file.path(TAB_CSV, "Table4_lapse_rates_site_mean.csv")
)

cat("✅ Table 4 guardada usando medias por sitio.\n")

# =============================================================================
# SUPPLEMENTARY / DIAGNOSTIC MODEL COMPARISON — LM, LMM, GAM
# Evaluates whether isotope–elevation relationships are linear or require
# site/project-level structure or nonlinear terms
# =============================================================================

cat("\n========== MODEL DIAGNOSTICS: LM vs LMM vs GAM ==========\n")

diag_data <- runoff_clean %>%
  filter(
    is.finite(d18O),
    is.finite(H50_km),
    !is.na(Project),
    !is.na(Id_site_u),
    !is.na(date)
  ) %>%
  mutate(
    Project = factor(Project),
    Site = factor(Id_site_u),
    Month_num = lubridate::month(date),
    Month = factor(Month_num, levels = 1:12, labels = month.abb),
    H50_z = as.numeric(scale(H50_km))
  ) %>%
  filter(!is.na(Month), is.finite(H50_z))

if (nrow(diag_data) < 10) {
  stop("❌ diag_data tiene muy pocas observaciones para comparar modelos.")
}

# -----------------------------
# 1) Linear Model
# -----------------------------
m_lm <- lm(
  d18O ~ H50_z * Project + Month,
  data = diag_data
)

# -----------------------------
# 2) Linear Mixed-Effects Model
# -----------------------------
m_lmm <- lme4::lmer(
  d18O ~ H50_z * Project + Month + (1 | Site),
  data = diag_data,
  REML = FALSE
)

# -----------------------------
# 3) Generalized Additive Model
# -----------------------------
# Nota: para el suavizado cíclico mensual se usa Month_num, no factor Month.
m_gam <- mgcv::gam(
  d18O ~ Project +
    s(H50_km, by = Project, k = 4) +
    s(Month_num, bs = "cc", k = 6) +
    s(Site, bs = "re"),
  data = diag_data,
  method = "REML",
  knots = list(Month_num = c(0.5, 12.5))
)

# -----------------------------
# Helper metrics robusto
# -----------------------------
get_model_fit_value <- function(model) {
  out <- NA_real_
  
  if (inherits(model, "gam")) {
    sm <- tryCatch(summary(model), error = function(e) NULL)
    if (!is.null(sm) && !is.null(sm$dev.expl)) {
      out <- as.numeric(sm$dev.expl)
    }
    
  } else if (inherits(model, "merMod")) {
    r2_obj <- tryCatch(performance::r2(model), error = function(e) NULL)
    if (!is.null(r2_obj) && "R2_conditional" %in% names(r2_obj)) {
      out <- as.numeric(r2_obj$R2_conditional)
    }
    
  } else if (inherits(model, "lm")) {
    sm <- tryCatch(summary(model), error = function(e) NULL)
    if (!is.null(sm) && !is.null(sm$r.squared)) {
      out <- as.numeric(sm$r.squared)
    }
  }
  
  out
}

model_metrics <- function(model, data, model_name) {
  pred <- tryCatch(
    predict(model, newdata = data, type = "response", allow.new.levels = TRUE),
    error = function(e) {
      tryCatch(
        predict(model, newdata = data, type = "response"),
        error = function(e2) rep(NA_real_, nrow(data))
      )
    }
  )
  
  obs <- data$d18O
  ok <- is.finite(obs) & is.finite(pred)
  
  tibble::tibble(
    Model = model_name,
    n = sum(ok),
    AIC = suppressWarnings(AIC(model)),
    BIC = suppressWarnings(BIC(model)),
    RMSE = ifelse(any(ok), sqrt(mean((obs[ok] - pred[ok])^2)), NA_real_),
    MAE = ifelse(any(ok), mean(abs(obs[ok] - pred[ok])), NA_real_),
    R2_or_Deviance_explained = get_model_fit_value(model)
  )
}

model_comparison <- dplyr::bind_rows(
  model_metrics(m_lm, diag_data, "LM"),
  model_metrics(m_lmm, diag_data, "LMM"),
  model_metrics(m_gam, diag_data, "GAM")
) %>%
  mutate(
    delta_AIC = AIC - min(AIC, na.rm = TRUE),
    AIC_weight = exp(-0.5 * delta_AIC) /
      sum(exp(-0.5 * delta_AIC), na.rm = TRUE)
  ) %>%
  arrange(AIC)

print(model_comparison)

readr::write_csv(
  model_comparison,
  file.path(TAB_CSV, "Table_model_diagnostics_LM_LMM_GAM.csv")
)

# -----------------------------
# ICC for LMM
# -----------------------------
icc_lmm <- tryCatch(
  performance::icc(m_lmm),
  error = function(e) NULL
)

icc_table <- tibble::tibble(
  Model = "LMM",
  ICC_adjusted = if (!is.null(icc_lmm) && "ICC_adjusted" %in% names(icc_lmm)) {
    as.numeric(icc_lmm$ICC_adjusted)
  } else {
    NA_real_
  },
  ICC_unadjusted = if (!is.null(icc_lmm) && "ICC_unadjusted" %in% names(icc_lmm)) {
    as.numeric(icc_lmm$ICC_unadjusted)
  } else {
    NA_real_
  }
)

print(icc_table)

readr::write_csv(
  icc_table,
  file.path(TAB_CSV, "Table_LMM_ICC_diagnostics.csv")
)

# -----------------------------
# Residual diagnostics
# -----------------------------
diag_residuals <- diag_data %>%
  mutate(
    pred_LM = as.numeric(predict(m_lm, newdata = diag_data)),
    resid_LM = d18O - pred_LM,
    
    pred_LMM = as.numeric(predict(m_lmm, newdata = diag_data, allow.new.levels = TRUE)),
    resid_LMM = d18O - pred_LMM,
    
    pred_GAM = as.numeric(predict(m_gam, newdata = diag_data, type = "response")),
    resid_GAM = d18O - pred_GAM
  )

readr::write_csv(
  diag_residuals,
  file.path(TAB_CSV, "Model_diagnostic_residuals_LM_LMM_GAM.csv")
)

# -----------------------------
# Plot observed vs predicted
# -----------------------------
obs_pred_long <- diag_residuals %>%
  select(Project, Site, date, d18O, pred_LM, pred_LMM, pred_GAM) %>%
  tidyr::pivot_longer(
    cols = starts_with("pred_"),
    names_to = "Model",
    values_to = "Predicted"
  ) %>%
  mutate(
    Model = dplyr::recode(
      Model,
      pred_LM = "LM",
      pred_LMM = "LMM",
      pred_GAM = "GAM"
    )
  )

p_obs_pred <- ggplot(obs_pred_long, aes(d18O, Predicted)) +
  geom_point(alpha = 0.65) +
  geom_abline(intercept = 0, slope = 1, linetype = 2) +
  facet_wrap(~Model) +
  labs(
    title = "Observed vs predicted δ¹⁸O for diagnostic models",
    x = expression("Observed " * delta^{18} * "O (‰)"),
    y = expression("Predicted " * delta^{18} * "O (‰)")
  ) +
  theme_pub()

save_fig(
  p_obs_pred,
  "Figure_model_diagnostics_observed_vs_predicted",
  w = 10,
  h = 5
)

cat("✅ Model diagnostics saved: LM, LMM, GAM, AIC, BIC, RMSE, MAE, ICC.\n")
# =============================================================================
# COEFICIENTES PARA END-MEMBERS GEOMORFOLÓGICOS DE FIGURE 4
# δ = intercept + slope * H50_km
# =============================================================================

lapse_coef <- table4 %>%
  select(Project, Isotope, Intercept, Slope_per_km, R2, n_sites)

write_csv(
  lapse_coef,
  file.path(TAB_CSV, "Table4_lapse_coefficients_for_geomorphic_endmembers.csv")
)

cat("\n✅ Coeficientes guardados para Figure 4 geomorfológica.\n")
print(lapse_coef)

predict_iso_from_h50 <- function(project, isotope, h50_km) {
  cf <- lapse_coef %>%
    filter(Project == project, Isotope == isotope)
  
  if (nrow(cf) != 1) {
    stop(paste("No se encontró coeficiente único para:", project, isotope))
  }
  
  cf$Intercept + cf$Slope_per_km * h50_km
}

plot_iso_h50 <- function(df, yvar, ylab_expr, panel_lab) {
  df_use <- df %>%
    filter(is.finite(.data[[yvar]]), is.finite(H50_m), !is.na(Project))
  
  ggplot(df_use, aes(H50_m, .data[[yvar]])) +
    geom_point(
      aes(size = n_months),
      shape = 21,
      fill = "grey70",
      color = "grey30",
      alpha = 0.85
    ) +
    geom_smooth(
      method = "lm",
      formula = y ~ x,
      se = TRUE,
      linewidth = 1.1,
      color = "#185FA5"
    ) +
    ggrepel::geom_text_repel(
      aes(label = Id_site_u),
      size = 3,
      min.segment.length = 0
    ) +
    facet_wrap(~Project, scales = "free_x") +
    scale_x_continuous(labels = label_number(big.mark = ",")) +
    scale_size_continuous(name = "Months", range = c(2, 5)) +
    labs(
      title = panel_lab,
      x = "Median elevation H50 (m a.s.l.)",
      y = ylab_expr
    ) +
    theme_pub()
}

p3_d18O <- plot_iso_h50(
  runoff_site_mean,
  "d18O",
  expression(delta^{18}*"O (‰)"),
  expression("(a) "*delta^{18}*"O–H"[50]*" relationship")
)

p3_dD <- plot_iso_h50(
  runoff_site_mean,
  "dD",
  expression(delta*"D (‰)"),
  expression("(b) "*delta*"D–H"[50]*" relationship")
)

fig3 <- p3_d18O / p3_dD +
  plot_annotation(
    title = "Isotope–elevation relationships by hydropower basin",
    subtitle = "Site-mean runoff isotopes using accumulated contributing catchment H50",
    theme = theme(plot.title = element_text(face = "bold", family = FONT))
  )

print(fig3)

save_fig(
  fig3,
  "Figure3_isotope_H50_by_basin_site_mean",
  w = 11,
  h = 10
)

# Captura de warnings del modelo LMM (si los hay)
w <- warnings()
if (length(w) > 0) {
  cat('\nWarnings del modelo LMM:\n')
  print(w[seq_len(min(5, length(w)))])
} else {
  cat('\nSin warnings en el modelo.\n')
}
# =============================================================================
# SECCIÓN 6 — TABLE 5 + FIGURE 4
# End-member transparency for Reviewer 2
#
# MODEL:
#   End-members predicted from basin-specific δ18O–H50 lapse-rate regressions
#   using the H50 of the reference catchments.
#
# MEASURED:
#   Monthly measured end-members from reference catchments:
#     MHP   : headwater = P5 Pita;      local = P1 Mashpi
#     CCSHP : headwater = P9 Antisana;  local = P14 Malo
#
# NOTE:
#   Mashpi and Malo are not assumed to be physical tributaries to the intakes.
#   They are adjacent low-elevation reference catchments used to represent
#   isotopically local runoff without high-elevation paramo/glacier influence.
# =============================================================================

cat("\n========== SECCIÓN 6: TABLE 5 + FIGURE 4 ==========\n")

# -----------------------------------------------------------------------------
# 6.1 Definición explícita de sitios usados en el modelo de mezcla
# -----------------------------------------------------------------------------

endmember_sites <- c(
  "5PITA",   # P5  Pita River — MHP high-elevation headwater
  "3MASH",   # P1  Mashpi River — MHP low-elevation local reference
  "7ANTI",   # P9  Antisana Stream — CCSHP high-elevation headwater
  "13MALO"   # P14 Malo River — CCSHP low-elevation local reference
)

intake_sites <- c(
  "2MAND",   # P2  Manduriacu intake
  "15COCA"   # P12 Coca intake
)

site_roles <- tibble::tibble(
  Project = c(
    "Manduriacu (a)",
    "Manduriacu (a)",
    "Manduriacu (a)",
    "Coca Codo Sinclair (b)",
    "Coca Codo Sinclair (b)",
    "Coca Codo Sinclair (b)"
  ),
  Id_site_u = c(
    "5PITA",
    "3MASH",
    "2MAND",
    "7ANTI",
    "13MALO",
    "15COCA"
  ),
  Point_ID = c(
    "P5",
    "P1",
    "P2",
    "P9",
    "P14",
    "P12"
  ),
  Site_label = c(
    "Pita River",
    "Mashpi River",
    "Manduriacu intake",
    "Antisana Stream",
    "Malo River",
    "Coca intake"
  ),
  Role = c(
    "headwater",
    "local",
    "intake",
    "headwater",
    "local",
    "intake"
  ),
  Endmember_type = c(
    "Measured and H50-modelled headwater end-member",
    "Measured and H50-modelled local reference end-member",
    "Mixture / hydropower intake",
    "Measured and H50-modelled headwater end-member",
    "Measured and H50-modelled local reference end-member",
    "Mixture / hydropower intake"
  ),
  Justification = c(
    "High-elevation Pita catchment representing paramo/glacier-influenced runoff in the Manduriacu system.",
    "Adjacent low-elevation Mashpi catchment used as local runoff reference without high-elevation paramo input.",
    "Integrated isotope signal at the Manduriacu hydropower intake.",
    "High-elevation Antisana stream representing paramo/glacier-influenced runoff in the Coca Codo Sinclair system.",
    "Adjacent low-elevation Malo catchment used as local runoff reference without high-elevation paramo input.",
    "Integrated isotope signal at the Coca Codo Sinclair hydropower intake."
  )
)

# Guardar tabla explícita para responder al Reviewer 2:
# "what precise values were used for the endmembers and how were they determined?"
table5A_endmember_definition <- site_roles %>%
  left_join(
    geom_new %>%
      select(
        Id_site_u, Basin, Area_km2, H50_m, H50_km,
        Altit_min, Altit_max, HI, Source_zone
      ),
    by = "Id_site_u"
  ) %>%
  arrange(Project, Role)

print(table5A_endmember_definition)

write_csv(
  table5A_endmember_definition,
  file.path(TAB_CSV, "Table5A_endmember_site_definition.csv")
)

cat("✅ Table 5A guardada: definición explícita de end-members e intakes.\n")
# =============================================================================
# 6.2 End-members geomorfológicos H50
#     MODEL = H50 de cuencas de referencia individuales
# =============================================================================

geom_endmembers <- geom_new %>%
  filter(Id_site_u %in% endmember_sites) %>%
  # site_roles tambien tiene columna Project -> se renombra antes del join
  # para evitar conflicto Project.x / Project.y que rompe mapply()
  left_join(
    site_roles %>% rename(Project_role = Project),
    by = "Id_site_u"
  ) %>%
  mutate(
    d18O_geom = mapply(
      function(project, h50_km) {
        predict_iso_from_h50(project, "d18O", h50_km)
      },
      Project, H50_km
    ),
    dD_geom = mapply(
      function(project, h50_km) {
        predict_iso_from_h50(project, "dD", h50_km)
      },
      Project, H50_km
    )
  ) %>%
  select(
    Project, Role, Id_site_u, Site_label,
    Area_km2, H50_m, H50_km, HI, Slope_mean,
    d18O_geom, dD_geom
  ) %>%
  arrange(Project, Role)

print(geom_endmembers)

write_csv(
  geom_endmembers,
  file.path(TAB_CSV, "Table5_geomorphic_H50_endmembers.csv")
)

get_geom_em <- function(proj, role, iso = "d18O_geom") {
  val <- geom_endmembers %>%
    filter(Project == proj, Role == role) %>%
    pull(.data[[iso]])
  
  if (length(val) != 1) {
    stop(paste("End-member geomorfológico no único:", proj, role, iso))
  }
  
  val
}
# =============================================================================
# 6.3 End-members medidos mensuales
# =============================================================================

em_monthly_measured <- runoff %>%
  filter(Id_site_u %in% endmember_sites) %>%
  # runoff ya tiene columna Project -> se renombra en site_roles para evitar .x/.y
  left_join(
    site_roles %>% rename(Project_role = Project),
    by = "Id_site_u"
  ) %>%
  mutate(Month = format(date, "%b-%y")) %>%
  select(Project, Month, date, Role, Id_site_u, Site_label, d18O, dD) %>%
  pivot_wider(
    names_from = Role,
    values_from = c(d18O, dD, Id_site_u, Site_label),
    names_sep = "_"
  )

write_csv(
  em_monthly_measured,
  file.path(TAB_CSV, "Table5_measured_monthly_endmembers.csv")
)

# =============================================================================
# 6.4 Intake mensual
# =============================================================================

intake_for_fig <- runoff %>%
  filter(Id_site_u %in% intake_sites, is.finite(d18O)) %>%
  mutate(Month = format(date, "%b-%y")) %>%
  group_by(Project, date, Month) %>%
  summarise(
    d18O_intake = mean(d18O, na.rm = TRUE),
    dD_intake   = mean(dD, na.rm = TRUE),
    .groups = "drop"
  )

# =============================================================================
# 6.5 MODEL
# =============================================================================

mix_model <- intake_for_fig %>%
  rowwise() %>%
  mutate(
    d18O_head  = get_geom_em(Project, "headwater", "d18O_geom"),
    d18O_local = get_geom_em(Project, "local", "d18O_geom"),
    denom      = d18O_head - d18O_local,
    f_head_raw = ifelse(
      is.finite(denom) & denom != 0,
      (d18O_intake - d18O_local) / denom,
      NA_real_
    ),
    inside_interval = d18O_intake >= min(d18O_head, d18O_local) &
      d18O_intake <= max(d18O_head, d18O_local),
    f_head  = pmin(pmax(f_head_raw, 0), 1),
    f_local = 1 - f_head,
    Type = "MODEL"
  ) %>%
  ungroup()

# =============================================================================
# 6.6 MEASURED
# =============================================================================

mix_measured <- intake_for_fig %>%
  left_join(em_monthly_measured, by = c("Project", "Month", "date")) %>%
  mutate(
    d18O_head  = d18O_headwater,
    d18O_local = d18O_local,
    denom      = d18O_head - d18O_local,
    f_head_raw = ifelse(
      is.finite(denom) & denom != 0,
      (d18O_intake - d18O_local) / denom,
      NA_real_
    ),
    inside_interval = d18O_intake >= pmin(d18O_head, d18O_local, na.rm = TRUE) &
      d18O_intake <= pmax(d18O_head, d18O_local, na.rm = TRUE),
    f_head  = pmin(pmax(f_head_raw, 0), 1),
    f_local = 1 - f_head,
    Type = "MEASURED"
  )

# =============================================================================
# 6.7 Comparación MODEL vs MEASURED
# =============================================================================

mix_compare <- bind_rows(
  mix_model %>%
    select(Project, date, Month, Type, d18O_intake, d18O_head, d18O_local,
           f_head_raw, f_head, f_local, inside_interval),
  mix_measured %>%
    select(Project, date, Month, Type, d18O_intake, d18O_head, d18O_local,
           f_head_raw, f_head, f_local, inside_interval)
) %>%
  arrange(Project, Type, date) %>%
  mutate(
    pct_head  = round(100 * f_head, 1),
    pct_local = round(100 * f_local, 1)
  )

print(mix_compare)

write_csv(
  mix_compare,
  file.path(TAB_CSV, "Table5_model_H50_vs_measured_mixing.csv")
)

outside_tbl <- mix_compare %>%
  filter(!inside_interval | f_head_raw < 0 | f_head_raw > 1)

if (nrow(outside_tbl) > 0) {
  print(outside_tbl)
  write_csv(
    outside_tbl,
    file.path(TAB_CSV, "Table5_outside_endmember_interval.csv")
  )
}

# =============================================================================
# 6.7B — TABLE 5B: COMPARACIÓN MODEL vs MEASURED
# =============================================================================

cat("\n========== TABLE 5B: MODEL vs MEASURED END-MEMBERS ==========\n")

table5B_model_vs_measured <- mix_compare %>%
  mutate(
    Type = as.character(Type)
  ) %>%
  group_by(Project, Type) %>%
  summarise(
    n_months = n(),

    d18O_head_mean  = round(mean(d18O_head, na.rm = TRUE), 2),
    d18O_head_sd    = round(sd(d18O_head, na.rm = TRUE), 2),

    d18O_local_mean = round(mean(d18O_local, na.rm = TRUE), 2),
    d18O_local_sd   = round(sd(d18O_local, na.rm = TRUE), 2),

    d18O_intake_mean = round(mean(d18O_intake, na.rm = TRUE), 2),
    d18O_intake_sd   = round(sd(d18O_intake, na.rm = TRUE), 2),

    mean_headwater_fraction = round(mean(f_head, na.rm = TRUE), 3),
    mean_headwater_pct      = round(100 * mean(f_head, na.rm = TRUE), 1),
    mean_local_pct          = round(100 * mean(f_local, na.rm = TRUE), 1),

    min_headwater_pct = round(100 * min(f_head, na.rm = TRUE), 1),
    max_headwater_pct = round(100 * max(f_head, na.rm = TRUE), 1),

    n_outside_interval = sum(
      !inside_interval | f_head_raw < 0 | f_head_raw > 1,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    Interpretation = case_when(
      Project == "Manduriacu (a)" & mean_headwater_pct >= 50 ~
        "Headwater-influenced intake composition",
      Project == "Manduriacu (a)" & mean_headwater_pct < 50 ~
        "Mixed intake composition",
      Project == "Coca Codo Sinclair (b)" & mean_headwater_pct < 50 ~
        "Local/intermediate-influenced intake composition",
      TRUE ~ "Mixed intake composition"
    )
  ) %>%
  arrange(Project, Type)

print(table5B_model_vs_measured)

write_csv(
  table5B_model_vs_measured,
  file.path(TAB_CSV, "Table5B_model_vs_measured_endmembers.csv")
)

cat("✅ Table 5B guardada: comparación MODEL vs MEASURED.\n")
table5B_model_vs_measured %>%
  select(
    Project,
    Type,
    mean_headwater_pct,
    mean_local_pct
  )

# =============================================================================
# 6.8 Caudales operacionales
# =============================================================================

q_tbl <- tribble(
  ~Project,                 ~date,                    ~Q_turb,
  "Manduriacu (a)",         as.Date("2017-09-01"),  90.0,
  "Manduriacu (a)",         as.Date("2017-10-01"),  90.1,
  "Manduriacu (a)",         as.Date("2017-11-01"),  73.2,
  "Manduriacu (a)",         as.Date("2017-12-01"), 107.0,
  "Manduriacu (a)",         as.Date("2018-01-01"), 181.4,
  "Coca Codo Sinclair (b)", as.Date("2017-09-01"), 149.2,
  "Coca Codo Sinclair (b)", as.Date("2017-10-01"), 142.2,
  "Coca Codo Sinclair (b)", as.Date("2017-11-01"), 140.3,
  "Coca Codo Sinclair (b)", as.Date("2017-12-01"), 155.4,
  "Coca Codo Sinclair (b)", as.Date("2018-01-01"), 149.7
)

table5 <- mix_model %>%
  left_join(q_tbl, by = c("Project", "date")) %>%
  filter(!is.na(Q_turb)) %>%
  mutate(
    Q_head    = round(f_head  * Q_turb, 1),
    Q_local   = round(f_local * Q_turb, 1),
    pct_head  = paste0(round(100 * f_head), "%"),
    pct_local = paste0(round(100 * f_local), "%")
  ) %>%
  arrange(Project, date)

print(table5)

write_csv(
  table5,
  file.path(TAB_CSV, "Table5_monthly_mixing_discharge_H50_model.csv")
)

# =============================================================================
# 6.9 Figure 4
# =============================================================================

month_order <- c(
  "sept-17", "oct-17", "nov-17", "dic-17", "ene-18", "feb-18",
  "mar-18", "abr-18", "may-18", "jun-18", "jul-18", "ago-18"
)

month_labels_en <- c(
  "Sep-17", "Oct-17", "Nov-17", "Dec-17", "Jan-18", "Feb-18",
  "Mar-18", "Apr-18", "May-18", "Jun-18", "Jul-18", "Aug-18"
)

fig4_data <- mix_compare %>%
  mutate(
    Month = factor(Month, levels = month_order, labels = month_labels_en),
    Project_lab = factor(
      Project,
      levels = c("Manduriacu (a)", "Coca Codo Sinclair (b)")
    ),
    Type = factor(
      Type,
      levels = c("MODEL", "MEASURED"),
      labels = c(
        "MODEL\narea-weighted H50",
        "MEASURED\nreference catchments"
      )
    )
  ) %>%
  filter(!is.na(Month))

fig4 <- ggplot(fig4_data, aes(x = Month)) +
  geom_linerange(
    aes(
      ymin = pmin(d18O_local, d18O_head, na.rm = TRUE),
      ymax = pmax(d18O_local, d18O_head, na.rm = TRUE)
    ),
    linewidth = 2.5,
    alpha = 0.55,
    colour = "grey55"
  ) +
  geom_point(
    aes(y = d18O_local),
    shape = 21,
    size = 2.9,
    fill = "#378ADD",
    colour = "black",
    stroke = 0.4
  ) +
  geom_point(
    aes(y = d18O_head),
    shape = 21,
    size = 2.9,
    fill = "#EF9F27",
    colour = "black",
    stroke = 0.4
  ) +
  geom_point(
    aes(y = d18O_intake, fill = f_head * 100),
    shape = 21,
    size = 4.0,
    colour = "black",
    stroke = 0.5
  ) +
  scale_fill_gradientn(
    colours = c("#440154", "#3B528B", "#21908C", "#5DC863", "#FDE725"),
    name = "Headwater\ncontribution (%)",
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  scale_y_reverse(
    breaks = seq(-16, -4, by = 2)
  ) +
  facet_grid(Project_lab ~ Type, scales = "free_x") +
  labs(
    title = expression("Monthly isotope mixing at hydropower intakes ("*delta^{18}*"O)"),
    subtitle = paste(
      "Orange = headwater end-member; blue = local end-member;",
      "colored intake point = estimated headwater contribution;",
      "grey bar = isotopic mixing interval"
    ),
    x = NULL,
    y = expression(delta^{18}*"O (‰)")
  ) +
  theme_pub() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    strip.text = element_text(face = "bold", size = 10),
    strip.background = element_rect(fill = "grey92", colour = "grey70"),
    legend.position = "bottom",
    plot.title = element_text(face = "bold", size = 12, hjust = 0),
    plot.subtitle = element_text(size = 9, colour = "grey35")
  )

print(fig4)

save_fig(
  fig4,
  "Figure4_monthly_runoff_mixing_area_weighted_H50_model_vs_measured",
  w = 12,
  h = 8
)

cat("\n✅ Figure 4 actualizada guardada.\n")

geom_endmembers

mix_compare %>%
  group_by(Project, Type) %>%
  summarise(
    mean_head  = mean(d18O_head, na.rm = TRUE),
    mean_local = mean(d18O_local, na.rm = TRUE),
    mean_fhead = mean(f_head, na.rm = TRUE),
    .groups = "drop"
  )
# =============================================================================
# SECCIÓN 7 — MODEL COMPARISON LM/LMM/GAM
# =============================================================================

cat("\n========== SECCIÓN 7: MODEL COMPARISON ==========\n")

# 7A — Lluvia
rain_mod <- rain_iaea %>%
  filter(is.finite(VO18),is.finite(Altitude),month %in% 1:12) %>%
  group_by(Name_of_site) %>% mutate(Altitude=median(Altitude,na.rm=TRUE)) %>% ungroup() %>%
  mutate(Alt_z=as.numeric(scale(Altitude)),
         month_f=factor(month,levels=1:12,labels=month.abb),
         month_num=month, Site=factor(Name_of_site))

m_rain_lm     <- lm(VO18~Alt_z+month_f+Site,data=rain_mod)
m_rain_lmm_ri <- lmer(VO18~Alt_z+month_f+(1|Site),data=rain_mod,REML=TRUE)
m_rain_lmm_rs <- tryCatch(
  lmer(VO18~Alt_z+month_f+(1+Alt_z|Site),data=rain_mod,REML=TRUE,
       control=lmerControl(check.conv.singular="ignore",optimizer="bobyqa")),
  error=function(e) NULL)
m_rain_gam <- gam(VO18~s(Alt_z,k=5)+s(month_num,bs="cc",k=12)+s(Site,bs="re"),
                  data=rain_mod,method="REML",knots=list(month_num=c(0.5,12.5)))

rain_models_tbl <- bind_rows(
  extract_model_metrics(m_rain_lm,    rain_mod,"VO18","LM",                        "Rainfall δ18O"),
  extract_model_metrics(m_rain_lmm_ri,rain_mod,"VO18","LMM random intercept",      "Rainfall δ18O"),
  if(!is.null(m_rain_lmm_rs))
    extract_model_metrics(m_rain_lmm_rs,rain_mod,"VO18","LMM random intercept+slope","Rainfall δ18O"),
  extract_model_metrics(m_rain_gam,   rain_mod,"VO18","GAM cyclic month",           "Rainfall δ18O")
) %>% arrange(AIC)
print(rain_models_tbl)
write_csv(rain_models_tbl,file.path(TAB_CSV,"TableS_model_comparison_rainfall.csv"))

# 7B — Escorrentía
runoff_mod <- runoff %>%
  filter(is.finite(d18O),is.finite(H50_km)) %>%
  mutate(month_num=month(date),
         month_f  =factor(month_num,levels=1:12,labels=month.abb),
         H50_z    =as.numeric(scale(H50_km)),
         Subset   =case_when(
           Id_site_u %in% c("1MAND","2MAND","3MASH","4GUAY","5PITA","8GUAC") ~"MHP + Mashpi",
           Id_site_u %in% c("7ANTI","13MALO","14COCA","15COCA","16SALA",
                             "17QUIJ","18TRIB","20SAN","22CARI")              ~"CCSHP + Malo",
           TRUE~NA_character_)) %>%
  filter(!is.na(Subset))

fit_suite_runoff <- function(df,subset_name) {
  d <- df %>% filter(Subset==subset_name) %>% mutate(Site=factor(Id_site_u))
  if (nrow(d)<8||n_distinct(d$Site)<2) return(NULL)
  m_lm     <- lm(d18O~H50_z+month_f+Site,data=d)
  m_lmm_ri <- lmer(d18O~H50_z+month_f+(1|Site),data=d,REML=TRUE)
  m_lmm_rs <- tryCatch(
    lmer(d18O~H50_z+month_f+(1+H50_z|Site),data=d,REML=TRUE,
         control=lmerControl(check.conv.singular="ignore",optimizer="bobyqa")),
    error=function(e) NULL)
  k_alt <- min(5,max(3,length(unique(d$H50_km))-1))
  m_gam <- gam(d18O~s(H50_km,k=k_alt)+s(month_num,bs="cc",k=12)+s(Site,bs="re"),
               data=d,method="REML",knots=list(month_num=c(0.5,12.5)))
  bind_rows(
    extract_model_metrics(m_lm,    d,"d18O","LM",                        subset_name),
    extract_model_metrics(m_lmm_ri,d,"d18O","LMM random intercept",      subset_name),
    if(!is.null(m_lmm_rs))
      extract_model_metrics(m_lmm_rs,d,"d18O","LMM random intercept+slope",subset_name),
    extract_model_metrics(m_gam,   d,"d18O","GAM cyclic month",           subset_name))
}

runoff_models_tbl <- bind_rows(
  fit_suite_runoff(runoff_mod,"MHP + Mashpi"),
  fit_suite_runoff(runoff_mod,"CCSHP + Malo")) %>% arrange(Dataset,AIC)
print(runoff_models_tbl)
write_csv(runoff_models_tbl,file.path(TAB_CSV,"TableS_model_comparison_runoff.csv"))

# 7C — Figura suplementaria GAM
eff_gam_alt <- as.data.frame(ggpredict(m_rain_gam,terms="Alt_z [all]"))
eff_gam_mon <- as.data.frame(ggpredict(m_rain_gam,terms="month_num [1:12]")) %>%
  mutate(x_lab=factor(round(x),levels=1:12,labels=month.abb))

p_gam_alt <- ggplot(eff_gam_alt,aes(x,predicted)) +
  geom_line() + geom_ribbon(aes(ymin=conf.low,ymax=conf.high),alpha=.2) +
  labs(title="A) Rainfall GAM: altitude smooth",x="Altitude (z-score)",
       y=expression(delta^{18}*"O (‰)")) + theme_pub()
p_gam_mon <- ggplot(eff_gam_mon,aes(x_lab,predicted,group=1)) +
  geom_line() + geom_point(size=1.6) +
  geom_ribbon(aes(ymin=conf.low,ymax=conf.high,group=1),alpha=.2) +
  labs(title="B) Rainfall GAM: cyclic month smooth",x="Month",
       y=expression(delta^{18}*"O (‰)")) + theme_pub()

figS <- p_gam_alt/p_gam_mon +
  plot_annotation(title=expression("Supplementary: GAM effects on rainfall "*delta^{18}*"O"),
                  theme=theme(plot.title=element_text(face="bold",family=FONT)))
print(figS)
save_fig(figS,"FigureS_rainfall_GAM_effects",w=10,h=8)

# 7D — Interpretación
cat("\n--- Interpretación modelos lluvia ---\n")
best_rain <- rain_models_tbl %>% slice_min(AIC,n=1,with_ties=FALSE)
cat("Mejor modelo (AIC):",best_rain$Model,"| AIC=",round(best_rain$AIC,2),
    "| R²=",round(best_rain$R2,3),"\n")

cat("\n--- Interpretación modelos escorrentía ---\n")
runoff_models_tbl %>% group_by(Dataset) %>%
  slice_min(AIC,n=1,with_ties=FALSE) %>% ungroup() %>%
  { for(i in seq_len(nrow(.)))
    cat("Dataset:",.[[i,"Dataset"]],"| Mejor:",.[[i,"Model"]],
        "| AIC=",round(.[[i,"AIC"]],2),"\n"); . } %>% invisible()

# =============================================================================
# RESUMEN FINAL
# =============================================================================

cat("\n╔══════════════════════════════════════════════════════╗\n")
cat("║         SCRIPT COMPLETADO EXITOSAMENTE              ║\n")
cat("╠══════════════════════════════════════════════════════╣\n")
cat("║  TABLAS en 03_outputs/tables/CSV/                   ║\n")
cat("║    ✅ Table1_LMWL_LCexcess.csv                      ║\n")
cat("║    ✅ Table2_monthly_weighted_MHP_CCSHP.csv          ║\n")
cat("║    ✅ Table3_deltaD_by_site.csv                      ║\n")
cat("║    ✅ Table3_runoff_LMWL_class.csv                   ║\n")
cat("║    ✅ Table4_lapse_rates.csv                         ║\n")
cat("║    ✅ Table5_monthly_mixing_discharge.csv            ║\n")
cat("║    ✅ TableS_model_comparison_rainfall.csv           ║\n")
cat("║    ✅ TableS_model_comparison_runoff.csv             ║\n")
cat("║  FIGURAS en PNG 300dpi y TIFF 600dpi                ║\n")
cat("║    ✅ Figure2  ✅ Figure3  ✅ Figure4  ✅ FigureS    ║\n")
cat("╚══════════════════════════════════════════════════════╝\n")
################
# GRAPH METHODOLOGY
# =============================================================================
# METHODOLOGICAL FLOWCHART — Isotopic analysis MHP & CCSHP
# Author: Paulina Lima
# Associated manuscript:
# Isotope-Based Source Attribution of Hydropower Inflows
# in Two Contrasting Ecuadorian Andean Watersheds
# Produces a publication-quality flowchart of the methodological workflow
# Output: PNG 300 dpi + TIFF 600 dpi
# Packages required: ggplot2, ggforce, patchwork, grid, gridExtra
# =============================================================================

# ---- Packages ----------------------------------------------------------------
pkgs <- c("ggplot2","ggforce","dplyr","tibble","grid","gridExtra","patchwork")
miss <- pkgs[!sapply(pkgs, requireNamespace, quietly = TRUE)]
if (length(miss)) install.packages(miss)
suppressPackageStartupMessages(invisible(lapply(pkgs, library, character.only = TRUE)))

# ---- Output paths ------------------------------------------------------------
BASE    <- "C:/Users/Pauli/Desktop/UCE/MANUSCRIPT_ISOTOPES/MANUSCRIPT_24_03_2026"
FIG_PNG <- file.path(BASE, "03_outputs/figures/PNG_300dpi")
FIG_TIF <- file.path(BASE, "03_outputs/figures/TIFF_600dpi")
for (d in c(FIG_PNG, FIG_TIF)) if (!dir.exists(d)) dir.create(d, recursive = TRUE)

FONT <- "serif"
# =============================================================================
# LAYOUT DESIGN
# =============================================================================
# The flowchart has 5 vertical columns (phases) and flows top-to-bottom
# with lateral connections between parallel processes.
#
# Column positions (x center):
#   Col 1 (x=1.5): Data inputs
#   Col 2 (x=4.0): Processing / analysis
#   Col 3 (x=6.5): Statistical models
#   Col 4 (x=9.0): Outputs
#   Col 5 (x=11.5): Interpretation
#
# =============================================================================

# ---- Color palette (publication-safe) ----------------------------------------
COL <- list(
  phase_data   = "#1F4E79",   # dark blue   — data inputs
  phase_proc   = "#2E75B6",   # mid blue    — processing
  phase_model  = "#0F6E56",   # dark green  — statistical models
  phase_output = "#7B3F00",   # dark brown  — outputs
  phase_interp = "#4A235A",   # dark purple — interpretation
  arrow        = "#444444",
  bg           = "white",
  box_border   = "#CCCCCC",
  phase_header = "#F2F2F2"
)

# ---- Node definitions --------------------------------------------------------
# Each node: id, label (text shown), x, y, width, height, fill, text_color
# y increases downward; width/height in data units

nodes <- tribble(
  ~id,           ~label,                                         ~x,    ~y,    ~w,    ~h,    ~fill,               ~tcol,    ~shape,
  # ── PHASE HEADERS ──
  "ph_data",     "DATA COLLECTION\n& PREPARATION",               1.5,   0.5,   2.4,   0.5,   COL$phase_data,     "white",  "rect",
  "ph_proc",     "ISOTOPIC\nANALYSIS",                           4.5,   0.5,   2.4,   0.5,   COL$phase_proc,     "white",  "rect",
  "ph_model",    "STATISTICAL\nMODELING",                        7.5,   0.5,   2.4,   0.5,   COL$phase_model,    "white",  "rect",
  "ph_output",   "RESULTS &\nOUTPUTS",                          10.5,   0.5,   2.4,   0.5,   COL$phase_output,   "white",  "rect",
  
  # ── DATA INPUTS ──
  "iaea",        "IAEA–INAMHI\nRainfall Network\n(1968–2020)\n10 stations | 538 records\nδ¹⁸O, δD, Precipitation",
  1.5,   2.0,   2.2,   1.4,   "#D6E4F0",          "#1F4E79", "rect",
  "field",       "Field Sampling\n(Sept 2017–Aug 2018)\n22 runoff sites\n5 rainfall stations\n281 surface + 47 rain samples",
  1.5,   4.2,   2.2,   1.4,   "#D6E4F0",          "#1F4E79", "rect",
  "geomorf",     "GIS / DEM Analysis\nWatershed delineation\nH₅₀, H_min, H_max\nArea, Stream order",
  1.5,   6.4,   2.2,   1.2,   "#D6E4F0",          "#1F4E79", "rect",
  "lab",         "Laboratory Analysis\nIRMS (δ¹⁸O ±0.2‰, δD ±0.2‰)\nUNM cross-check\nPicarro L1102-I",
  1.5,   8.3,   2.2,   1.2,   "#D6E4F0",          "#1F4E79", "rect",
  
  # ── ISOTOPIC ANALYSIS ──
  "lmwl",        "LMWL Calculation\n(Eq. 1–2)\nOLS: δD = a·δ¹⁸O + b\nMHP: a=7.99, b=10.43\nCCSHP: a=8.39, b=14.51",
  4.5,   2.0,   2.2,   1.4,   "#C8E6C9",          "#0F6E56", "rect",
  "lcexcess",    "LC-excess\n(Eq. 3)\nLC = δD − a·δ¹⁸O − b\nDeviation from LMWL",
  4.5,   4.0,   2.2,   1.1,   "#C8E6C9",          "#0F6E56", "rect",
  "awmean",      "Amount-weighted\nMeans (Eq. 4–5)\nδ_w = Σ(δᵢ·Pᵢ)/ΣPᵢ\nd-excess by month",
  4.5,   5.8,   2.2,   1.1,   "#C8E6C9",          "#0F6E56", "rect",
  "hyps",        "Hypsometric Analysis\n(Eq. 6)\nHI = (H₅₀−H_min)/(H_max−H_min)\nClassification: High / Low",
  4.5,   7.6,   2.2,   1.1,   "#C8E6C9",          "#0F6E56", "rect",
  "lapse",       "Isotope–Elevation\nRegressions (Eq. 7)\nδ = a + b·H₅₀(km)\nLapse rates ‰ km⁻¹",
  4.5,   9.4,   2.2,   1.1,   "#C8E6C9",          "#0F6E56", "rect",
  
  # ── STATISTICAL MODELS ──
  "lm",          "Linear Model (LM)\n(Eq. 10)\nδ¹⁸O ~ Alt_z + month\n+ site (fixed effects)",
  7.5,   2.5,   2.2,   1.1,   "#FFF9C4",          "#5D4037", "rect",
  "lmm",         "Linear Mixed Model\n(LMM) (Eq. 11)\nδ¹⁸O ~ Alt_z + month\n+ (1|site)",
  7.5,   4.3,   2.2,   1.1,   "#FFF9C4",          "#5D4037", "rect",
  "gam",         "Generalized Additive\nModel (GAM) (Eq. 12)\ns(Alt) + s(month,cc)\n+ s(site,re)  AIC-best",
  7.5,   6.1,   2.2,   1.1,   "#FFF9C4",          "#5D4037", "rect",
  "mix",         "Two-Component\nMixing Model (Eq. 8)\nf_head = (δ_intake−δ_local)\n         /(δ_head−δ_local)",
  7.5,   8.1,   2.2,   1.1,   "#FFF9C4",          "#5D4037", "rect",
  "qdisch",      "Discharge Partitioning\n(Eq. 9)\nQ_local = [(δ_hydro−δ_head)·Q]\n         /(δ_local−δ_head)",
  7.5,   9.9,   2.2,   1.1,   "#FFF9C4",          "#5D4037", "rect",
  
  # ── OUTPUTS ──
  "t1",          "Table 1\nLMWL coefficients\n+ LC-excess statistics",
  10.5,   2.0,   2.2,   0.9,   "#F3E5F5",          "#4A235A", "rect",
  "t2",          "Table 2\nMonthly amount-weighted\nδ¹⁸O, δD, d-excess",
  10.5,   3.5,   2.2,   0.9,   "#F3E5F5",          "#4A235A", "rect",
  "t3",          "Table 3\nΔδD deviation\nfrom LMWL by site",
  10.5,   5.0,   2.2,   0.9,   "#F3E5F5",          "#4A235A", "rect",
  "t4",          "Table 4\nIsotope lapse rates\nvs H₅₀ (R², slope, CI)",
  10.5,   6.5,   2.2,   0.9,   "#F3E5F5",          "#4A235A", "rect",
  "t5",          "Table 5\nMonthly mixing +\ndischarge partition",
  10.5,   8.0,   2.2,   0.9,   "#F3E5F5",          "#4A235A", "rect",
  "f2",          "Figure 2\nLMM 2×2: fixed effects\naltitude · month · sites",
  10.5,   9.3,   2.2,   0.9,   "#FCE4EC",          "#880E4F", "rect",
  "f3",          "Figure 3\nδ¹⁸O & δD vs H₅₀\nOLS ± 95% CI by basin",
  10.5,  10.6,   2.2,   0.9,   "#FCE4EC",          "#880E4F", "rect",
  "f4",          "Figure 4\nMixing model vs measured\n2×2 panels MHP/CCSHP",
  10.5,  11.9,   2.2,   0.9,   "#FCE4EC",          "#880E4F", "rect",
)

# ---- Arrow definitions -------------------------------------------------------
# Each arrow: from_id, to_id, label (optional)
arrows <- tribble(
  ~from,      ~to,        ~lbl,
  # Data → Analysis
  "iaea",     "lmwl",     "",
  "iaea",     "awmean",   "",
  "field",    "lcexcess", "",
  "field",    "lapse",    "",
  "geomorf",  "hyps",     "",
  "geomorf",  "lapse",    "",
  "lab",      "lapse",    "",
  
  # Analysis → Models
  "lmwl",     "lcexcess", "",
  "lmwl",     "lm",       "",
  "lmwl",     "lmm",      "",
  "lmwl",     "gam",      "",
  "awmean",   "lm",       "",
  "hyps",     "lapse",    "",
  "lapse",    "mix",      "",
  "lapse",    "qdisch",   "",
  
  # Models → Outputs
  "lmwl",     "t1",       "",
  "awmean",   "t2",       "",
  "lcexcess", "t3",       "",
  "lapse",    "t4",       "",
  "mix",      "t5",       "",
  "qdisch",   "t5",       "",
  "lm",       "f2",       "",
  "lmm",      "f2",       "",
  "gam",      "f2",       "",
  "lapse",    "f3",       "",
  "mix",      "f4",       "",
)

# =============================================================================
# BUILD PLOT
# =============================================================================

# Helper: get node center coordinates
node_xy <- function(nid) {
  r <- nodes[nodes$id == nid, ]
  list(x = r$x, y = r$y, w = r$w, h = r$h)
}

# Helper: compute arrow endpoints (edge of boxes)
get_arrow <- function(from_id, to_id) {
  f <- node_xy(from_id)
  t <- node_xy(to_id)
  # Determine direction and pick appropriate edge
  dx <- t$x - f$x
  dy <- t$y - f$y
  if (abs(dx) >= abs(dy)) {
    # horizontal dominant
    if (dx > 0) {
      x0 <- f$x + f$w/2; y0 <- f$y
      x1 <- t$x - t$w/2; y1 <- t$y
    } else {
      x0 <- f$x - f$w/2; y0 <- f$y
      x1 <- t$x + t$w/2; y1 <- t$y
    }
  } else {
    # vertical dominant
    if (dy > 0) {
      x0 <- f$x; y0 <- f$y + f$h/2
      x1 <- t$x; y1 <- t$y - t$h/2
    } else {
      x0 <- f$x; y0 <- f$y - f$h/2
      x1 <- t$x; y1 <- t$y + t$h/2
    }
  }
  data.frame(x=x0, y=y0, xend=x1, yend=y1)
}

# Build arrow dataframe
arrow_df <- do.call(rbind, lapply(seq_len(nrow(arrows)), function(i) {
  tryCatch(get_arrow(arrows$from[i], arrows$to[i]), error = function(e) NULL)
}))

# =============================================================================
# GGPLOT
# =============================================================================

p <- ggplot() +
  
  # ── Phase background strips ──
  annotate("rect", xmin=0.3,  xmax=2.7,  ymin=-0.1, ymax=12.5,
           fill="#EBF3FB", alpha=0.35, color=NA) +
  annotate("rect", xmin=3.3,  xmax=5.7,  ymin=-0.1, ymax=12.5,
           fill="#E8F5E9", alpha=0.35, color=NA) +
  annotate("rect", xmin=6.3,  xmax=8.7,  ymin=-0.1, ymax=12.5,
           fill="#FFFDE7", alpha=0.35, color=NA) +
  annotate("rect", xmin=9.3,  xmax=11.7, ymin=-0.1, ymax=12.5,
           fill="#F3E5F5", alpha=0.35, color=NA) +
  
  # ── Arrows ──
  geom_segment(data = arrow_df,
               aes(x=x, y=y, xend=xend, yend=yend),
               arrow = arrow(length=unit(0.18,"cm"), type="closed"),
               color = COL$arrow, linewidth = 0.45, alpha = 0.7) +
  
  # ── Boxes ──
  geom_tile(data = nodes,
            aes(x=x, y=y, width=w, height=h, fill=fill),
            color = "#888888", linewidth = 0.4) +
  scale_fill_identity() +
  
  # ── Text labels ──
  geom_text(data = nodes,
            aes(x=x, y=y, label=label, color=tcol),
            size = 2.5, family = FONT,
            lineheight = 1.15, fontface = "plain") +
  scale_color_identity() +
  
  # ── Phase header labels (above strips) ──
  annotate("text", x=1.5,  y=-0.25, label="DATA\nINPUTS",
           size=3.2, fontface="bold", family=FONT, color=COL$phase_data,   hjust=0.5) +
  annotate("text", x=4.5,  y=-0.25, label="ISOTOPIC\nANALYSIS",
           size=3.2, fontface="bold", family=FONT, color=COL$phase_proc,   hjust=0.5) +
  annotate("text", x=7.5,  y=-0.25, label="STATISTICAL\nMODELS",
           size=3.2, fontface="bold", family=FONT, color=COL$phase_model,  hjust=0.5) +
  annotate("text", x=10.5, y=-0.25, label="RESULTS &\nOUTPUTS",
           size=3.2, fontface="bold", family=FONT, color=COL$phase_output, hjust=0.5) +
  
  # ── Main title ──
  annotate("text", x=6.0, y=12.3,
           label="Methodological workflow — Isotopic analysis of hydropower runoff sources",
           size=4.5, fontface="bold", family=FONT, color="#1F4E79", hjust=0.5) +
  annotate("text", x=6.0, y=12.0,
           label="Manduriacu (MHP) and Coca Codo Sinclair (CCSHP) — Northern Andes, Ecuador",
           size=3.2, fontface="italic", family=FONT, color="#444444", hjust=0.5) +
  
  # ── Legend ──
  annotate("rect",  xmin=0.4, xmax=0.9,  ymin=11.2, ymax=11.5, fill="#D6E4F0", color="#888888", linewidth=0.3) +
  annotate("text",  x=1.0, y=11.35, label="Data input",   size=2.4, hjust=0, family=FONT, color="#1F4E79") +
  annotate("rect",  xmin=2.2, xmax=2.7,  ymin=11.2, ymax=11.5, fill="#C8E6C9", color="#888888", linewidth=0.3) +
  annotate("text",  x=2.8, y=11.35, label="Analysis",     size=2.4, hjust=0, family=FONT, color="#0F6E56") +
  annotate("rect",  xmin=4.0, xmax=4.5,  ymin=11.2, ymax=11.5, fill="#FFF9C4", color="#888888", linewidth=0.3) +
  annotate("text",  x=4.6, y=11.35, label="Model",        size=2.4, hjust=0, family=FONT, color="#5D4037") +
  annotate("rect",  xmin=5.8, xmax=6.3,  ymin=11.2, ymax=11.5, fill="#F3E5F5", color="#888888", linewidth=0.3) +
  annotate("text",  x=6.4, y=11.35, label="Table output", size=2.4, hjust=0, family=FONT, color="#4A235A") +
  annotate("rect",  xmin=7.6, xmax=8.1,  ymin=11.2, ymax=11.5, fill="#FCE4EC", color="#888888", linewidth=0.3) +
  annotate("text",  x=8.2, y=11.35, label="Figure output",size=2.4, hjust=0, family=FONT, color="#880E4F") +
  
  # ── Axes and theme ──
  scale_x_continuous(limits=c(0.2, 12.0)) +
  scale_y_continuous(limits=c(-0.6, 12.6), trans="reverse") +
  coord_cartesian(clip="off") +
  theme_void() +
  theme(
    plot.background  = element_rect(fill="white", color=NA),
    panel.background = element_rect(fill="white", color=NA),
    plot.margin      = margin(15, 15, 15, 15)
  )

# =============================================================================
# SAVE
# =============================================================================

ggsave(file.path(FIG_PNG, "Flowchart_methodology.png"),
       p, width=16, height=14, dpi=300, bg="white")
ggsave(file.path(FIG_TIF, "Flowchart_methodology.tiff"),
       p, width=16, height=14, dpi=600, bg="white",
       compression="lzw")

cat("✅ Flowchart saved:\n")
cat("   PNG 300dpi →", file.path(FIG_PNG, "Flowchart_methodology.png"), "\n")
cat("   TIFF 600dpi →", file.path(FIG_TIF, "Flowchart_methodology.tiff"), "\n")
###############
# =============================================================================
# SECCIÓN X — COMPARISON OF MODELLED VS MEASURED END-MEMBERS
# Residual analysis + improved Table 5 + Figure 6
# =============================================================================

library(dplyr)
library(readr)
library(ggplot2)
library(mgcv)
library(patchwork)
library(broom)
library(scales)

# -----------------------------------------------------------------------------
# 1) BUILD RESIDUAL DATASET: MODEL vs MEASURED
# -----------------------------------------------------------------------------

# mix_model and mix_measured are assumed to exist from the previous section
# mix_model    -> end-members from annual/regression-derived values
# mix_measured -> end-members from monthly measured values

df_residuals <- mix_model %>%
  select(
    Project, date, Month,
    f_head_model = f_head,
    d18O_intake_model = d18O_intake,
    d18O_head_model   = d18O_head,
    d18O_local_model  = d18O_local
  ) %>%
  left_join(
    mix_measured %>%
      select(
        Project, date, Month,
        f_head_measured = f_head,
        d18O_head_measured = d18O_head,
        d18O_local_measured = d18O_local
      ),
    by = c("Project", "date", "Month")
  ) %>%
  mutate(
    residual = f_head_model - f_head_measured,
    project = factor(Project),
    month_num = as.integer(format(date, "%m")),
    Basin = case_when(
      Project == "Manduriacu (a)" ~ "MHP",
      Project == "Coca Codo Sinclair (b)" ~ "CCSHP",
      TRUE ~ NA_character_
    )
  )

# Atmospheric covariates from Table 2
table2_cov <- table2 %>%
  mutate(month_num = match(Month, month.abb)) %>%
  select(Basin, month_num, d18O_w, Precip_mm, d_excess)

df_residuals <- df_residuals %>%
  left_join(table2_cov, by = c("Basin", "month_num")) %>%
  mutate(
    # scale predictors to avoid very different scales
    d18O_w_z   = as.numeric(scale(d18O_w)),
    Precip_z   = as.numeric(scale(Precip_mm)),
    d_excess_z = as.numeric(scale(d_excess))
  )

cat("\n========== RESIDUAL DATASET ==========\n")
print(df_residuals, n = 24)
str(df_residuals)

# -----------------------------------------------------------------------------
# 2) ERROR ANALYSIS
# NOTE: ICC is not appropriate here because only 2 projects are available.
# We therefore use linear models and compare AIC, MAE, RMSE, and R².
# -----------------------------------------------------------------------------

m0 <- lm(residual ~ project, data = df_residuals)

m1 <- lm(
  residual ~ project + d18O_w_z + Precip_z + d_excess_z,
  data = df_residuals
)

m2 <- lm(
  residual ~ project * (d18O_w_z + Precip_z + d_excess_z),
  data = df_residuals
)

cat("\n========== MODEL SUMMARIES ==========\n")
cat("\n--- m0 ---\n")
print(summary(m0))

cat("\n--- m1 ---\n")
print(summary(m1))

cat("\n--- m2 ---\n")
print(summary(m2))

cat("\n========== AIC COMPARISON ==========\n")
print(AIC(m0, m1, m2))

# error metrics before and after best model correction
mae_before  <- mean(abs(df_residuals$residual), na.rm = TRUE)
rmse_before <- sqrt(mean(df_residuals$residual^2, na.rm = TRUE))

df_residuals <- df_residuals %>%
  mutate(
    pred_m2 = predict(m2, newdata = .),
    resid_corrected = resid(m2)
  )

mae_after  <- mean(abs(df_residuals$resid_corrected), na.rm = TRUE)
# error metrics before and after best model correction
mae_before  <- mean(abs(df_residuals$residual), na.rm = TRUE)
rmse_before <- sqrt(mean(df_residuals$residual^2, na.rm = TRUE))

df_residuals <- df_residuals %>%
  mutate(
    pred_m2 = predict(m2, newdata = .),
    resid_corrected = residual - pred_m2
  )

mae_after  <- mean(abs(df_residuals$resid_corrected), na.rm = TRUE)
rmse_after <- sqrt(mean(df_residuals$resid_corrected^2, na.rm = TRUE))

cat("\n========== ERROR METRICS ==========\n")
cat("MAE before:", round(mae_before, 3), "\n")
cat("RMSE before:", round(rmse_before, 3), "\n")
cat("MAE after :", round(mae_after, 3), "\n")
cat("RMSE after:", round(rmse_after, 3), "\n")

error_summary <- tibble(
  Model = c("Raw residuals", "Corrected residuals (m2)"),
  MAE = c(mae_before, mae_after),
  RMSE = c(rmse_before, rmse_after)
)

write_csv(
  error_summary,
  file.path(TAB_CSV, "Table_error_model_summary.csv")
)

cat("\n========== ERROR METRICS ==========\n")
cat("MAE before:", round(mae_before, 3), "\n")
cat("RMSE before:", round(rmse_before, 3), "\n")
cat("MAE after :", round(mae_after, 3), "\n")
cat("RMSE after:", round(rmse_after, 3), "\n")

# Optional summary table
error_summary <- tibble(
  Model = c("Raw residuals", "Corrected residuals (m2)"),
  MAE = c(mae_before, mae_after),
  RMSE = c(rmse_before, rmse_after)
)

write_csv(error_summary, file.path(TAB_CSV, "Table_error_model_summary.csv"))

# -----------------------------------------------------------------------------
# 3) TABLE 5 — PRIMARY TABLE (MODELLED END-MEMBERS)
# This is the recommended main table for the manuscript
# -----------------------------------------------------------------------------

gen_tbl <- tribble(
  ~Project, ~date, ~Generation_GWh, ~Share_pct,
  "Manduriacu (a)", as.Date("2017-09-01"), 0.5, 2,
  "Manduriacu (a)", as.Date("2017-10-01"), 0.5, 3,
  "Manduriacu (a)", as.Date("2017-11-01"), 0.5, 3,
  "Manduriacu (a)", as.Date("2017-12-01"), 0.7, 3,
  "Manduriacu (a)", as.Date("2018-01-01"), 0.9, 5,
  "Coca Codo Sinclair (b)", as.Date("2017-09-01"), 19.3, 98,
  "Coca Codo Sinclair (b)", as.Date("2017-10-01"), 18.3, 97,
  "Coca Codo Sinclair (b)", as.Date("2017-11-01"), 17.9, 97,
  "Coca Codo Sinclair (b)", as.Date("2017-12-01"), 19.8, 97,
  "Coca Codo Sinclair (b)", as.Date("2018-01-01"), 19.2, 95
)

site_tbl <- tribble(
  ~Project, ~Intake_site, ~Headwater_site, ~Local_site,
  "Manduriacu (a)", "P2",  "P5",  "P1",
  "Coca Codo Sinclair (b)", "P12", "P9", "P14"
)

table5_model <- table5 %>%
  mutate(
    Month_en = case_when(
      Month %in% c("sept-17", "Sep-17") ~ "Sep-17",
      Month %in% c("oct-17",  "Oct-17") ~ "Oct-17",
      Month %in% c("nov-17",  "Nov-17") ~ "Nov-17",
      Month %in% c("dic-17",  "Dec-17") ~ "Dec-17",
      Month %in% c("ene-18",  "Jan-18") ~ "Jan-18",
      TRUE ~ Month
    ),
    date = case_when(
      Month_en == "Sep-17" ~ as.Date("2017-09-01"),
      Month_en == "Oct-17" ~ as.Date("2017-10-01"),
      Month_en == "Nov-17" ~ as.Date("2017-11-01"),
      Month_en == "Dec-17" ~ as.Date("2017-12-01"),
      Month_en == "Jan-18" ~ as.Date("2018-01-01"),
      TRUE ~ as.Date(NA)
    ),
    Headwater_pct = as.numeric(gsub("%", "", pct_head)),
    Local_pct     = as.numeric(gsub("%", "", pct_local))
  ) %>%
  left_join(gen_tbl, by = c("Project", "date")) %>%
  left_join(site_tbl, by = "Project") %>%
  transmute(
    Project,
    Month = Month_en,
    Intake_site,
    Headwater_site,
    Local_site,
    d18O_intake    = round(d18O_intake, 2),
    d18O_headwater = round(d18O_head, 2),
    d18O_local     = round(d18O_local, 2),
    Q_intake       = round(Q_turb, 1),
    Q_headwater    = round(Q_head, 1),
    Q_local        = round(Q_local, 1),
    Headwater_pct  = round(Headwater_pct, 0),
    Local_pct      = round(Local_pct, 0),
    Generation_GWh = round(Generation_GWh, 1),
    Share_pct      = round(Share_pct, 0)
  )

cat("\n========== TABLE 5 — MODELLED END-MEMBERS ==========\n")
print(table5_model, n = 20, width = Inf)

write_csv(
  table5_model,
  file.path(TAB_CSV, "Table5_modelled_endmembers.csv")
)
# -----------------------------------------------------------------------------
# 4) TABLE 5b — SENSITIVITY / VALIDATION TABLE (MEASURED END-MEMBERS)
# -----------------------------------------------------------------------------

q_tbl_clean <- q_tbl %>%
  rename(Q_intake = Q_turb)

table5_measured <- mix_measured %>%
  left_join(q_tbl_clean, by = c("Project", "date")) %>%
  mutate(
    Q_headwater = round(f_head * Q_intake, 1),
    Q_local     = round((1 - f_head) * Q_intake, 1),
    Headwater_pct = round(100 * f_head, 0),
    Local_pct     = round(100 * (1 - f_head), 0)
  ) %>%
  left_join(gen_tbl, by = c("Project", "date")) %>%
  left_join(site_tbl, by = "Project") %>%
  transmute(
    Project,
    Month,
    Intake_site,
    Headwater_site,
    Local_site,
    d18O_intake    = round(d18O_intake, 2),
    d18O_headwater = round(d18O_head, 2),
    d18O_local     = round(d18O_local, 2),
    Q_intake       = round(Q_intake, 1),
    Q_headwater    = Q_headwater,
    Q_local        = Q_local,
    Headwater_pct  = Headwater_pct,
    Local_pct      = Local_pct,
    Generation_GWh = round(Generation_GWh, 1),
    Share_pct      = round(Share_pct, 0)
  )

cat("\n========== TABLE 5b — MEASURED END-MEMBERS ==========\n")
print(table5_measured, n = 20)
write_csv(table5_measured, file.path(TAB_CSV, "Table5_measured_endmembers.csv"))

# -----------------------------------------------------------------------------
# 5) COMPARISON TABLE — MODEL vs MEASURED
# -----------------------------------------------------------------------------

table5_compare <- df_residuals %>%
  transmute(
    Project,
    Month,
    f_head_model    = round(100 * f_head_model, 1),
    f_head_measured = round(100 * f_head_measured, 1),
    residual_pctpt  = round(100 * residual, 1)
  )

cat("\n========== TABLE 5c — MODEL vs MEASURED ==========\n")
print(table5_compare, n = 24)
write_csv(table5_compare, file.path(TAB_CSV, "Table5_model_vs_measured_comparison.csv"))

# -----------------------------------------------------------------------------
# 6) FIGURE 6 — d2H vs H50 with LM, GAM and GAMM-like smooth
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# 6) FIGURE 6 — d2H vs H50 with LM, GAM and GAMM-like smooth
# -----------------------------------------------------------------------------

# use the dataset already built in your script:
# runoff_clean should exist and include dD, H50_m, Project, Id_site_u

df_plot <- runoff_clean %>%
  filter(is.finite(dD), is.finite(H50_m)) %>%
  mutate(
    label_mean = case_when(
      Project == "Manduriacu (a)" ~ "Mashpi annual mean",
      Project == "Coca Codo Sinclair (b)" ~ "Malo annual mean",
      TRUE ~ as.character(Project)
    )
  )

make_panel <- function(dat, panel_title = NULL) {
  
  dat <- dat %>%
    mutate(Site = factor(Id_site_u))
  
  # -----------------------------
  # MODELS
  # -----------------------------
  m_lm <- lm(dD ~ H50_m, data = dat)
  
  k_use <- min(5, max(3, length(unique(dat$H50_m)) - 1))
  m_gam <- gam(dD ~ s(H50_m, k = k_use), data = dat, method = "REML")
  
  m_gamm <- gam(
    dD ~ s(H50_m, k = k_use) + s(Site, bs = "re"),
    data = dat,
    method = "REML"
  )
  
  # -----------------------------
  # PREDICTION GRID
  # -----------------------------
  newdat <- data.frame(
    H50_m = seq(min(dat$H50_m, na.rm = TRUE),
                max(dat$H50_m, na.rm = TRUE),
                length.out = 200)
  )
  
  pred_lm <- predict(m_lm, newdata = newdat, se.fit = TRUE)
  pred_gam <- predict(m_gam, newdata = newdat, se.fit = TRUE)
  pred_gamm <- predict(
    m_gamm,
    newdata = cbind(newdat, Site = dat$Site[1]),
    se.fit = TRUE
  )
  
  # -----------------------------
  # DATA FRAMES FOR PLOTTING
  # -----------------------------
  lm_df <- newdat %>%
    mutate(
      fit   = pred_lm$fit,
      lo    = fit - 1.96 * pred_lm$se.fit,
      hi    = fit + 1.96 * pred_lm$se.fit,
      Model = "LM"
    )
  
  gam_df <- newdat %>%
    mutate(
      fit   = pred_gam$fit,
      lo    = fit - 1.96 * pred_gam$se.fit,
      hi    = fit + 1.96 * pred_gam$se.fit,
      Model = "GAM"
    )
  
  gamm_df <- newdat %>%
    mutate(
      fit   = pred_gamm$fit,
      lo    = fit - 1.96 * pred_gamm$se.fit,
      hi    = fit + 1.96 * pred_gamm$se.fit,
      Model = "GAMM"
    )
  
  pred_all <- bind_rows(lm_df, gam_df, gamm_df)
  
  # -----------------------------
  # ANNUAL MEAN POINT
  # -----------------------------
  annual_mean <- dat %>%
    summarise(
      H50_m = mean(H50_m, na.rm = TRUE),
      dD    = mean(dD, na.rm = TRUE)
    ) %>%
    mutate(label = unique(dat$label_mean)[1])
  
  # -----------------------------
  # COLORS
  # -----------------------------
  model_cols <- c(
    "LM"   = "#1f77b4",
    "GAM"  = "#d62728",
    "GAMM" = "#2ca02c"
  )
  
  # -----------------------------
  # PLOT
  # -----------------------------
  ggplot() +
    # observed points
    geom_point(
      data = dat,
      aes(x = H50_m, y = dD),
      shape = 21, size = 2.2,
      fill = "grey70", colour = "black",
      alpha = 0.85
    ) +
    
    # confidence bands
    geom_ribbon(
      data = pred_all,
      aes(x = H50_m, ymin = lo, ymax = hi, fill = Model),
      alpha = 0.14
    ) +
    
    # fitted lines
    geom_line(
      data = pred_all,
      aes(x = H50_m, y = fit, colour = Model),
      linewidth = 1.1
    ) +
    
    # annual mean point
    geom_point(
      data = annual_mean,
      aes(x = H50_m, y = dD),
      shape = 24, size = 4,
      fill = "black", colour = "black"
    ) +
    
    # annual mean label
    geom_text(
      data = annual_mean,
      aes(x = H50_m, y = dD, label = label),
      nudge_y = 2,
      size = 3.8,
      fontface = "plain"
    ) +
    
    scale_colour_manual(values = model_cols, name = "Model") +
    scale_fill_manual(values = model_cols, name = "Model") +
    
    labs(
      title = panel_title,
      x = "Mean basin elevation H50 (m)",
      y = expression("Monthly runoff "*delta^2*H~"(‰)")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5),
      axis.title = element_text(face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
}

p1 <- make_panel(
  filter(df_plot, Project == "Manduriacu (a)"),
  "Manduriacu (a)"
)

p2 <- make_panel(
  filter(df_plot, Project == "Coca Codo Sinclair (b)"),
  "Coca Codo Sinclair (b)"
)

fig6 <- p1 + p2 + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

print(fig6)

ggsave(
  file.path(FIG_PNG, "Figure6_d2H_vs_H50_LM_GAM_GAMM.png"),
  fig6,
  width = 12, height = 7, dpi = 300, bg = "white"
)

ggsave(
  file.path(FIG_TIF, "Figure6_d2H_vs_H50_LM_GAM_GAMM.tiff"),
  fig6,
  width = 12, height = 7, dpi = 600, bg = "white",
  compression = "lzw"
)

# -----------------------------------------------------------------------------
# 7) RECOMMENDED TEXTUAL INTERPRETATION
# -----------------------------------------------------------------------------

cat("\n========== RECOMMENDED INTERPRETATION ==========\n")
cat("Primary analysis: modelled end-members (integrated elevation-band signal)\n")
cat("Validation/sensitivity: measured monthly end-members\n")
cat("Best residual model: m2 if it shows lowest AIC and lower MAE/RMSE\n")
# -----------------------------------------------------------------------------
# 7) RECOMMENDED TEXTUAL INTERPRETATION
# -----------------------------------------------------------------------------

cat("\n========== RECOMMENDED INTERPRETATION ==========\n")
cat("Primary analysis: modelled end-members (integrated elevation-band signal)\n")
cat("Validation/sensitivity: measured monthly end-members\n")
cat("Best residual model: m2 if it shows lowest AIC and lower MAE/RMSE\n")
# ============================
# EXUTORIOS REALES Y RÍOS PRINCIPALES
# ============================

if (!requireNamespace("ggnewscale", quietly = TRUE)) {
  install.packages("ggnewscale")
}
library(ggnewscale)

shp_exutorios <- "C:/ECUADOR_DATOS/RASTER COCA MAND/CUENCAS_ACUMULADAS_RUNOFF_V2/exutorios_usados_v2_20260529_080104.shp"
shp_rios <- "C:/ECUADOR_DATOS/RASTER COCA MAND/MANUSCRIPT_FIGURE_1/rios_principales_manuscript.shp"

exut <- st_read(shp_exutorios)
rios <- st_read(shp_rios)
# NOTA: re-proyeccion de exut/rios se hace mas abajo, despues de crear cuencas2

# DEM: se carga mas abajo con la ruta correcta (PATH_DEM) despues de crear cuencas2

# exut2 se construye correctamente mas abajo, despues de definir cuencas2
# ============================
# MAPA FINAL — EXUTORIOS, END-MEMBERS Y SIMILITUD ISOTÓPICA
# ============================

if (!requireNamespace("terra", quietly = TRUE)) install.packages("terra")
if (!requireNamespace("tidyterra", quietly = TRUE)) install.packages("tidyterra")
if (!requireNamespace("ggnewscale", quietly = TRUE)) install.packages("ggnewscale")
if (!requireNamespace("ggrepel", quietly = TRUE)) install.packages("ggrepel")

library(terra)
library(tidyterra)
library(ggnewscale)
library(ggrepel)
# tidyterra enmascara dplyr::filter — se restaura explicitamente
filter <- dplyr::filter

# ----------------------------
# 1) Preparar cuencas con Id_site_u
# ----------------------------

cuencas2 <- cuencas_acc %>%
  mutate(
    Punto = toupper(trimws(Nombre_P)),
    Id_site_u = case_when(
      Punto == "P2A"  ~ "1MAND",
      Punto == "P2"   ~ "2MAND",
      Punto == "P1"   ~ "3MASH",
      Punto == "P4"   ~ "4GUAY",
      Punto == "P5"   ~ "5PITA",
      Punto == "P15"  ~ "6MICA",
      Punto == "P9"   ~ "7ANTI",
      Punto == "P8"   ~ "8GUAC",
      Punto == "P7"   ~ "9GRAN",
      Punto == "P16"  ~ "10TURB",
      Punto == "P17"  ~ "11TURB",
      Punto == "P18"  ~ "12RESR",
      Punto == "P14"  ~ "13MALO",
      Punto == "P12A" ~ "14COCA",
      Punto == "P12"  ~ "15COCA",
      Punto == "P13"  ~ "16SALA",
      Punto == "P11"  ~ "17QUIJ",
      Punto == "P10"  ~ "18TRIB",
      Punto == "P19"  ~ "19RESR",
      Punto == "P3"   ~ "20SAN",
      Punto == "P20"  ~ "21DSCH",
      Punto == "P6"   ~ "22CARI",
      TRUE ~ Punto
    )
  )

# Re-proyectar exut y rios al mismo CRS que cuencas2 (si es necesario)
if (!is.null(exut) && st_crs(exut) != st_crs(cuencas2))
  exut <- st_transform(exut, st_crs(cuencas2))
if (exists('rios') && !is.null(rios) && st_crs(rios) != st_crs(cuencas2))
  rios <- st_transform(rios, st_crs(cuencas2))

main_project_catchments <- cuencas2 %>%
  filter(Id_site_u %in% c("2MAND", "15COCA"))

reference_catchments <- cuencas2 %>%
  filter(Id_site_u %in% c("3MASH", "5PITA", "13MALO", "7ANTI")) %>%
  mutate(
    Ref_type = case_when(
      Id_site_u %in% c("5PITA", "7ANTI") ~ "Headwater reference catchment",
      Id_site_u %in% c("3MASH", "13MALO") ~ "Local reference catchment",
      TRUE ~ "Reference catchment"
    )
  )

# ----------------------------
# 2) Preparar DEM
# ----------------------------

PATH_DEM <- "C:/ECUADOR_DATOS/RASTER COCA MAND/MOSAICO_FINAL/COPERNICUS_DEM_30M_COCA_MAND.tif"

if (!file.exists(PATH_DEM)) {
  stop("No existe el DEM en PATH_DEM. Revisa la ruta.")
}

dem <- terra::rast(PATH_DEM)

if (!terra::same.crs(dem, terra::vect(cuencas2))) {
  dem <- terra::project(dem, terra::crs(terra::vect(cuencas2)))
}

map_ext <- terra::ext(
  718000,
  905000,
  9910000,
  10072000
)

dem_plot <- terra::crop(dem, map_ext)
dem_plot <- terra::aggregate(dem_plot, fact = 3, fun = mean, na.rm = TRUE)

cat("✅ dem_plot creado correctamente con extensión ampliada.\n")

# ----------------------------
# 3) Preparar exutorios
# ----------------------------

PATH_EXUTORIOS <- "C:/ECUADOR_DATOS/RASTER COCA MAND/CUENCAS_ACUMULADAS_RUNOFF_V2/exutorios_usados_v2_20260529_080104.shp"

if (!file.exists(PATH_EXUTORIOS)) {
  stop("No existe el shapefile de exutorios. Revisa PATH_EXUTORIOS.")
}

exut <- st_read(PATH_EXUTORIOS, quiet = FALSE)

if (!"ID" %in% names(exut)) {
  stop("El shapefile de exutorios no tiene columna ID.")
}

exut2 <- exut %>%
  mutate(
    Punto = toupper(trimws(as.character(ID))),
    Id_site_u = case_when(
      Punto == "P2A"  ~ "1MAND",
      Punto == "P2"   ~ "2MAND",
      Punto == "P1"   ~ "3MASH",
      Punto == "P4"   ~ "4GUAY",
      Punto == "P5"   ~ "5PITA",
      Punto == "P15"  ~ "6MICA",
      Punto == "P9"   ~ "7ANTI",
      Punto == "P8"   ~ "8GUAC",
      Punto == "P7"   ~ "9GRAN",
      Punto == "P16"  ~ "10TURB",
      Punto == "P17"  ~ "11TURB",
      Punto == "P18"  ~ "12RESR",
      Punto == "P14"  ~ "13MALO",
      Punto == "P12A" ~ "14COCA",
      Punto == "P12"  ~ "15COCA",
      Punto == "P13"  ~ "16SALA",
      Punto == "P11"  ~ "17QUIJ",
      Punto == "P10"  ~ "18TRIB",
      Punto == "P19"  ~ "19RESR",
      Punto == "P3"   ~ "20SAN",
      Punto == "P20"  ~ "21DSCH",
      Punto == "P6"   ~ "22CARI",
      TRUE ~ Punto
    ),
    Project = case_when(
      Id_site_u %in% c(
        "1MAND", "2MAND", "3MASH", "4GUAY",
        "5PITA", "8GUAC", "9GRAN", "22CARI"
      ) ~ "Manduriacu (a)",
      Id_site_u %in% c(
        "7ANTI", "10TURB", "11TURB", "12RESR",
        "13MALO", "14COCA", "15COCA", "16SALA",
        "17QUIJ", "18TRIB", "19RESR", "20SAN", "21DSCH"
      ) ~ "Coca Codo Sinclair (b)",
      Id_site_u == "6MICA" ~ "Interbasin transfer / Quito water supply",
      TRUE ~ NA_character_
    ),
    d18O_map = suppressWarnings(as.numeric(d18OMean)),
    d18O_obs = d18O_map
  ) %>%
  left_join(
    geom_new %>%
      select(Id_site_u, H50_m, H50_km, Source_zone),
    by = "Id_site_u"
  )

# ----------------------------
# 4) End-members y mezcla
# ----------------------------

em_values <- exut2 %>%
  st_drop_geometry() %>%
  filter(Id_site_u %in% c("5PITA", "3MASH", "7ANTI", "13MALO")) %>%
  mutate(
    em_role = case_when(
      Id_site_u %in% c("5PITA", "7ANTI") ~ "head",
      Id_site_u %in% c("3MASH", "13MALO") ~ "local",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Project), !is.na(em_role), is.finite(d18O_map)) %>%
  select(Project, em_role, d18O_map) %>%
  tidyr::pivot_wider(
    names_from = em_role,
    values_from = d18O_map,
    names_prefix = "d18O_"
  )

headwater_reference <- c("5PITA", "7ANTI")
local_reference <- c("3MASH", "13MALO")
auxiliary_high_elevation <- c("6MICA", "18TRIB", "19RESR")

exut2 <- exut2 %>%
  left_join(em_values, by = "Project") %>%
  mutate(
    f_head_raw = (d18O_map - d18O_local) / (d18O_head - d18O_local),
    f_head = pmin(pmax(f_head_raw, 0), 1),
    f_local = 1 - f_head,
    
    f_head_pct = case_when(
      Id_site_u == "2MAND"  ~ 57,
      Id_site_u == "15COCA" ~ 16,
      Id_site_u %in% headwater_reference ~ 100,
      Id_site_u %in% local_reference ~ 0,
      TRUE ~ 100 * f_head
    ),
    
    f_local_pct = 100 - f_head_pct,
    
    Source_similarity = case_when(
      Id_site_u %in% c("2MAND", "15COCA") ~ "Hydropower intake",
      Id_site_u %in% headwater_reference ~ "Headwater reference source",
      Id_site_u %in% local_reference ~ "Local reference source",
      Id_site_u %in% auxiliary_high_elevation ~ "High-elevation auxiliary point",
      TRUE ~ "Runoff exutorio"
    ),
    
    pt_size = case_when(
      Source_similarity == "Hydropower intake" ~ 7.5,
      Source_similarity %in% c(
        "Headwater reference source",
        "Local reference source",
        "High-elevation auxiliary point"
      ) ~ 6.5,
      TRUE ~ 5.2
    ),
    
    label_id = case_when(
      Id_site_u == "2MAND"  ~ "P2",
      Id_site_u == "15COCA" ~ "P12",
      Id_site_u == "3MASH"  ~ "P1\nMashpi",
      Id_site_u == "5PITA"  ~ "P5\nPita",
      Id_site_u == "6MICA"  ~ "P15\nMica Lake",
      Id_site_u == "7ANTI"  ~ "P9\nAntisana Stream",
      Id_site_u == "8GUAC"  ~ "P8\nGuachalá",
      Id_site_u == "9GRAN"  ~ "P7\nGranoble",
      Id_site_u == "13MALO" ~ "P14\nMalo",
      Id_site_u == "18TRIB" ~ "P10\nPapallacta tributary",
      Id_site_u == "19RESR" ~ "P19\nPapallacta Lake",
      Id_site_u == "22CARI" ~ "P6",
      Id_site_u == "4GUAY"  ~ "P4",
      Id_site_u == "17QUIJ" ~ "P11",
      Id_site_u == "16SALA" ~ "P13",
      Id_site_u == "20SAN"  ~ "P3",
      TRUE ~ Punto
    ),
    
    label_map = case_when(
      Id_site_u == "6MICA" ~ label_id,
      is.finite(f_head_pct) ~ paste0(label_id, "\n", round(f_head_pct, 0), "% H"),
      TRUE ~ label_id
    )
  )

# ----------------------------
# 5) Tabla de salida y etiquetas
# ----------------------------

exut_summary <- exut2 %>%
  st_drop_geometry() %>%
  select(
    Id_site_u, Punto, Project, H50_m, Source_zone,
    d18O_obs, d18O_map,
    d18O_head, d18O_local,
    f_head_pct, f_local_pct,
    Source_similarity
  ) %>%
  arrange(Project, desc(f_head_pct))

write.csv(
  exut_summary,
  file.path(TAB_CSV, "Exutorios_source_similarity_d18O.csv"),
  row.names = FALSE
)

print(exut_summary)

key_labels <- exut2 %>%
  filter(Id_site_u %in% c(
    "2MAND", "15COCA",
    "3MASH", "5PITA",
    "6MICA", "7ANTI", "18TRIB", "19RESR",
    "8GUAC", "9GRAN",
    "22CARI", "13MALO",
    "4GUAY", "17QUIJ", "16SALA", "20SAN"
  ))

# ----------------------------
# 6) Mapa final
# ----------------------------

out_png <- file.path(FIG_PNG, "Figure_runoff_isotope_similarity_map.png")
out_tif <- file.path(FIG_TIF, "Figure_runoff_isotope_similarity_map.tiff")

has_rios <- exists("rios")

p <- ggplot() +
  tidyterra::geom_spatraster(data = dem_plot, alpha = 0.95) +
  scale_fill_gradientn(
    name = "Elevation\n(m a.s.l.)",
    colours = c("grey98", "grey88", "grey72", "grey55", "grey35"),
    na.value = NA
  ) +
  ggnewscale::new_scale_fill() +
  geom_sf(
    data = cuencas2,
    fill = NA,
    color = "grey82",
    linewidth = 0.14,
    show.legend = FALSE
  )

if (has_rios) {
  p <- p +
    geom_sf(
      data = rios,
      color = "grey35",
      linewidth = 0.35,
      alpha = 0.55,
      show.legend = FALSE
    )
}

p <- p +
  geom_sf(
    data = main_project_catchments,
    fill = NA,
    color = "black",
    linewidth = 1.1,
    show.legend = FALSE
  ) +
  geom_sf(
    data = reference_catchments,
    aes(linetype = Ref_type),
    fill = NA,
    color = "black",
    linewidth = 0.75
  ) +
  geom_sf(
    data = exut2,
    aes(
      fill = f_head_pct,
      shape = Source_similarity,
      size = pt_size
    ),
    color = "black",
    stroke = 0.65
  ) +
  scale_size_identity() +
  scale_fill_gradientn(
    name = "Headwater-like\nsimilarity (%)",
    colours = c("#d95f02", "#66a61e", "#1f78b4"),
    limits = c(0, 100),
    breaks = c(0, 25, 50, 75, 100)
  ) +
  ggrepel::geom_label_repel(
    data = cbind(
      st_drop_geometry(key_labels),
      st_coordinates(key_labels)
    ),
    aes(x = X, y = Y, label = label_map),
    size = 2.6,
    fontface = "bold",
    label.size = 0.12,   # reemplaza linewidth (obsoleto en ggrepel)
    fill = "white",
    alpha = 0.93,
    min.segment.length = 0,
    max.overlaps = Inf,
    box.padding = 0.45,
    point.padding = 0.30,
    show.legend = FALSE
  ) +
  annotate(
    "label",
    x = 744000,
    y = 10059000,
    label = "P2 Manduriacu\nHeadwater-like ≈ 57%\nLocal-like ≈ 43%",
    size = 2.8,
    fontface = "bold",
    hjust = 0,
    linewidth = 0.18,
    fill = "white",
    alpha = 0.95
  ) +
  annotate(
    "label",
    x = 820000,
    y = 10059000,
    label = "P12 Coca\nHeadwater-like ≈ 16%\nLocal/intermediate-like ≈ 84%",
    size = 2.8,
    fontface = "bold",
    hjust = 0,
    linewidth = 0.18,
    fill = "white",
    alpha = 0.95
  ) +
  scale_linetype_manual(
    name = "Reference catchments",
    values = c(
      "Headwater reference catchment" = "solid",
      "Local reference catchment" = "dashed",
      "Reference catchment" = "dotdash"
    )
  ) +
  scale_shape_manual(
    name = "Exutorio type",
    values = c(
      "Hydropower intake" = 23,
      "Headwater reference source" = 24,
      "Local reference source" = 21,
      "High-elevation auxiliary point" = 21,
      "Runoff exutorio" = 21
    )
  ) +
  labs(
    title = "Runoff isotope similarity along hydropower contributing catchments",
    subtitle = expression(
      "Exutorios classified by mean annual runoff " *
        delta^{18} * "O using the same end-members as the mixing model"
    ),
    x = NULL,
    y = NULL,
    caption =
      "Filled symbols show percentage similarity to the headwater end-member. Only P5/Pita and P9/Antisana Stream are used as headwater reference end-members. P10, P19 and P15 are labelled as high-elevation auxiliary points, not as mixing end-members."
  ) +
  coord_sf(
    xlim = c(718000, 905000),
    ylim = c(9910000, 10072000),
    expand = FALSE,
    clip = "off"
  ) +
  theme_bw(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5, size = 9),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.caption = element_text(size = 8),
    plot.margin = margin(10, 25, 10, 25)
  )

print(p)

ggsave(
  out_png,
  p,
  width = 12,
  height = 7,
  dpi = 300,
  bg = "white"
)

ggsave(
  out_tif,
  p,
  width = 12,
  height = 7,
  dpi = 600,
  bg = "white",
  compression = "lzw"
)

cat("✅ Figura guardada en:\n", out_png, "\n", out_tif, "\n")
########################
# =============================================================================
# SUPPLEMENTARY TABLE S2 — Runoff isotope sampling sites
# =============================================================================

cat("\n========== SUPPLEMENTARY TABLE S2: Runoff sampling sites ==========\n")

supp_table_s2 <- geom_new %>%
  mutate(
    Site_code = Nombre_P,
    
    Site_name = case_when(
      Id_site_u == "1MAND"  ~ "Manduriacu intake downstream",
      Id_site_u == "2MAND"  ~ "Manduriacu intake upstream",
      Id_site_u == "3MASH"  ~ "Mashpi River",
      Id_site_u == "4GUAY"  ~ "Guayllabamba River",
      Id_site_u == "5PITA"  ~ "Pita River",
      Id_site_u == "6MICA"  ~ "Mica Lake",
      Id_site_u == "7ANTI"  ~ "Antisana Stream",
      Id_site_u == "8GUAC"  ~ "Guachalá River",
      Id_site_u == "9GRAN"  ~ "Granoble River",
      Id_site_u == "10TURB" ~ "Turbine discharge 1",
      Id_site_u == "11TURB" ~ "Turbine discharge 2",
      Id_site_u == "12RESR" ~ "Coca reservoir",
      Id_site_u == "13MALO" ~ "Malo River",
      Id_site_u == "14COCA" ~ "Coca intake upstream",
      Id_site_u == "15COCA" ~ "Coca intake downstream",
      Id_site_u == "16SALA" ~ "Salado River",
      Id_site_u == "17QUIJ" ~ "Quijos River",
      Id_site_u == "18TRIB" ~ "Tributary",
      Id_site_u == "19RESR" ~ "Papallacta Lake",
      Id_site_u == "20SAN"  ~ "San Pedro River",
      Id_site_u == "21DSCH" ~ "Downstream discharge",
      Id_site_u == "22CARI" ~ "Cariacu River",
      TRUE ~ Id_site_u
    ),
    
    Site_category = case_when(
      Id_site_u %in% c("5PITA", "22CARI", "7ANTI", "19RESR") ~ "High-elevation headwater",
      Id_site_u %in% c("3MASH", "13MALO") ~ "Low-elevation local analog",
      Id_site_u %in% c("2MAND", "14COCA") ~ "Hydropower intake upstream",
      Id_site_u %in% c("1MAND", "15COCA") ~ "Hydropower intake downstream",
      Id_site_u %in% c("10TURB", "11TURB", "21DSCH") ~ "Turbine/discharge point",
      Id_site_u %in% c("12RESR", "19RESR", "6MICA") ~ "Lake/reservoir",
      TRUE ~ "Tributary"
    ),
    
    Basin_domain = case_when(
      Basin == "MHP" ~ "Pacific-draining Manduriacu domain",
      Basin == "CCSHP" ~ "Amazon-draining Coca Codo Sinclair domain",
      Basin == "Interbasin transfer" ~ "Interbasin transfer / Quito water supply",
      TRUE ~ Basin
    )
  ) %>%
  transmute(
    Site_code,
    Internal_ID = Id_site_u,
    Site_name,
    Basin = Basin,
    Basin_domain,
    Site_category,
    Area_km2 = round(Area_km2, 3),
    H50_m = round(H50_m, 1),
    Altit_min_m = round(Altit_min, 1),
    Altit_max_m = round(Altit_max, 1),
    HI = round(HI, 3),
    Source_zone
  ) %>%
  arrange(Basin, Site_code)

print(supp_table_s2)

write_csv(
  supp_table_s2,
  file.path(TAB_SUP, "Supplementary_Table_S2_Runoff_sampling_sites.csv")
)

cat("✅ Supplementary Table S2 guardada en:", 
    file.path(TAB_SUP, "Supplementary_Table_S2_Runoff_sampling_sites.csv"), "\n")
########
# =============================================================================
# SUPPLEMENTARY TABLE S3 — DEM-derived catchment geomorphology
# =============================================================================

cat("\n========== SUPPLEMENTARY TABLE S3: DEM-derived catchment geomorphology ==========\n")

supp_table_s3 <- geom_new %>%
  transmute(
    Site_code = Nombre_P,
    Internal_ID = Id_site_u,
    Basin,
    Project,
    Area_km2 = round(Area_km2, 3),
    Perimeter_km = round(Perim_km, 3),
    Kc = round(Kc, 3),
    Hmin_m = round(Altit_min, 1),
    H50_m = round(H50_m, 1),
    Hmax_m = round(Altit_max, 1),
    Relief_m = round(Relieve_m, 1),
    HI = round(HI, 3),
    Mean_slope = round(Slope_mean, 3),
    Max_slope = round(Slope_max, 3),
    Main_channel_length_km = round(LenRio_km, 3),
    Drainage_density_km_km2 = round(Dren_dens, 3),
    Source_zone
  ) %>%
  arrange(Basin, Site_code)

print(supp_table_s3)

write_csv(
  supp_table_s3,
  file.path(TAB_SUP, "Supplementary_Table_S3_DEM_derived_catchment_geomorphology.csv")
)

cat("✅ Supplementary Table S3 guardada en:",
    file.path(TAB_SUP, "Supplementary_Table_S3_DEM_derived_catchment_geomorphology.csv"), "\n")

