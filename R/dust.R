## Generated by dust2 (version 0.3.24) - do not edit
SEIRV_Model <- structure(
  function() get("SEIRV_Model"),
  class = "dust_system_generator",
  name = "SEIRV_Model",
  package = "YEP",
  path = NULL,
  parameters = data.frame(
    name = c("time_inc", "t_incubation", "t_latent", "t_infectious", "FOI_spillover", "R0", "N_age", "vacc_rate_daily", "vaccine_efficacy", "year0", "S_0", "E_0", "I_0", "R_0", "V_0", "dP1_all", "dP2_all", "n_years", "n_t_pts"),
    type = c("real_type", "real_type", "real_type", "real_type", "real_type", "real_type", "int", "real_type", "real_type", "real_type", "real_type", "real_type", "real_type", "real_type", "real_type", "real_type", "real_type", "int", "int"),
    constant = c(FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE),
    required = c(TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE),
    rank = c(0L, 0L, 0L, 0L, 1L, 1L, 0L, 2L, 0L, 0L, 1L, 1L, 1L, 1L, 1L, 2L, 2L, 0L, 0L)),
  properties = list(
    time_type = "discrete",
    has_compare = FALSE,
    has_adjoint = FALSE),
  default_dt = 1)
