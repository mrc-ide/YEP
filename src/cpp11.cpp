// Generated by cpp11: do not edit by hand
// clang-format off


#include "cpp11/declarations.hpp"
#include <R_ext/Visibility.h>

// SEIRV_Model.cpp
cpp11::sexp dust_SEIRV_Model_gpu_info();
extern "C" SEXP _YEP_dust_SEIRV_Model_gpu_info() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_SEIRV_Model_gpu_info());
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_alloc(cpp11::list r_pars, bool pars_multi, cpp11::sexp r_time, cpp11::sexp r_n_particles, int n_threads, cpp11::sexp r_seed, bool deterministic, cpp11::sexp gpu_config, cpp11::sexp ode_control);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_alloc(SEXP r_pars, SEXP pars_multi, SEXP r_time, SEXP r_n_particles, SEXP n_threads, SEXP r_seed, SEXP deterministic, SEXP gpu_config, SEXP ode_control) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_alloc(cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(r_pars), cpp11::as_cpp<cpp11::decay_t<bool>>(pars_multi), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_time), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_n_particles), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_seed), cpp11::as_cpp<cpp11::decay_t<bool>>(deterministic), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(gpu_config), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(ode_control)));
  END_CPP11
}
// SEIRV_Model.cpp
cpp11::sexp dust_cpu_SEIRV_Model_capabilities();
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_capabilities() {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_capabilities());
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_run(SEXP ptr, cpp11::sexp r_time_end);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_run(SEXP ptr, SEXP r_time_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_run(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_time_end)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_simulate(SEXP ptr, cpp11::sexp time_end);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_simulate(SEXP ptr, SEXP time_end) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_simulate(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(time_end)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_set_index(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_set_index(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_set_index(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_update_state(SEXP ptr, SEXP r_pars, SEXP r_state, SEXP r_time, SEXP r_set_initial_state, SEXP index, SEXP reset_step_size);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_update_state(SEXP ptr, SEXP r_pars, SEXP r_state, SEXP r_time, SEXP r_set_initial_state, SEXP index, SEXP reset_step_size) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_update_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_pars), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_time), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_set_initial_state), cpp11::as_cpp<cpp11::decay_t<SEXP>>(index), cpp11::as_cpp<cpp11::decay_t<SEXP>>(reset_step_size)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_state(SEXP ptr, SEXP r_index);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_state(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(r_index)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_time(SEXP ptr);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_time(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_time(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// SEIRV_Model.cpp
void dust_cpu_SEIRV_Model_reorder(SEXP ptr, cpp11::sexp r_index);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_reorder(SEXP ptr, SEXP r_index) {
  BEGIN_CPP11
    dust_cpu_SEIRV_Model_reorder(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(r_index));
    return R_NilValue;
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_resample(SEXP ptr, cpp11::doubles r_weights);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_resample(SEXP ptr, SEXP r_weights) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_resample(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::doubles>>(r_weights)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_rng_state(SEXP ptr, bool first_only, bool last_only);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_rng_state(SEXP ptr, SEXP first_only, SEXP last_only) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<bool>>(first_only), cpp11::as_cpp<cpp11::decay_t<bool>>(last_only)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_set_rng_state(SEXP ptr, cpp11::raws rng_state);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_set_rng_state(SEXP ptr, SEXP rng_state) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_set_rng_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::raws>>(rng_state)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_set_data(SEXP ptr, cpp11::list data, bool shared);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_set_data(SEXP ptr, SEXP data, SEXP shared) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_set_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<cpp11::list>>(data), cpp11::as_cpp<cpp11::decay_t<bool>>(shared)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_compare_data(SEXP ptr);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_compare_data(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_compare_data(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_filter(SEXP ptr, SEXP time_end, bool save_trajectories, cpp11::sexp time_snapshot, cpp11::sexp min_log_likelihood);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_filter(SEXP ptr, SEXP time_end, SEXP save_trajectories, SEXP time_snapshot, SEXP min_log_likelihood) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_filter(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(time_end), cpp11::as_cpp<cpp11::decay_t<bool>>(save_trajectories), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(time_snapshot), cpp11::as_cpp<cpp11::decay_t<cpp11::sexp>>(min_log_likelihood)));
  END_CPP11
}
// SEIRV_Model.cpp
void dust_cpu_SEIRV_Model_set_n_threads(SEXP ptr, int n_threads);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_set_n_threads(SEXP ptr, SEXP n_threads) {
  BEGIN_CPP11
    dust_cpu_SEIRV_Model_set_n_threads(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<int>>(n_threads));
    return R_NilValue;
  END_CPP11
}
// SEIRV_Model.cpp
int dust_cpu_SEIRV_Model_n_state(SEXP ptr);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_n_state(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_n_state(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}
// SEIRV_Model.cpp
void dust_cpu_SEIRV_Model_set_stochastic_schedule(SEXP ptr, SEXP time);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_set_stochastic_schedule(SEXP ptr, SEXP time) {
  BEGIN_CPP11
    dust_cpu_SEIRV_Model_set_stochastic_schedule(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr), cpp11::as_cpp<cpp11::decay_t<SEXP>>(time));
    return R_NilValue;
  END_CPP11
}
// SEIRV_Model.cpp
SEXP dust_cpu_SEIRV_Model_ode_statistics(SEXP ptr);
extern "C" SEXP _YEP_dust_cpu_SEIRV_Model_ode_statistics(SEXP ptr) {
  BEGIN_CPP11
    return cpp11::as_sexp(dust_cpu_SEIRV_Model_ode_statistics(cpp11::as_cpp<cpp11::decay_t<SEXP>>(ptr)));
  END_CPP11
}

extern "C" {
static const R_CallMethodDef CallEntries[] = {
    {"_YEP_dust_SEIRV_Model_gpu_info",                    (DL_FUNC) &_YEP_dust_SEIRV_Model_gpu_info,                    0},
    {"_YEP_dust_cpu_SEIRV_Model_alloc",                   (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_alloc,                   9},
    {"_YEP_dust_cpu_SEIRV_Model_capabilities",            (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_capabilities,            0},
    {"_YEP_dust_cpu_SEIRV_Model_compare_data",            (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_compare_data,            1},
    {"_YEP_dust_cpu_SEIRV_Model_filter",                  (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_filter,                  5},
    {"_YEP_dust_cpu_SEIRV_Model_n_state",                 (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_n_state,                 1},
    {"_YEP_dust_cpu_SEIRV_Model_ode_statistics",          (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_ode_statistics,          1},
    {"_YEP_dust_cpu_SEIRV_Model_reorder",                 (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_reorder,                 2},
    {"_YEP_dust_cpu_SEIRV_Model_resample",                (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_resample,                2},
    {"_YEP_dust_cpu_SEIRV_Model_rng_state",               (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_rng_state,               3},
    {"_YEP_dust_cpu_SEIRV_Model_run",                     (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_run,                     2},
    {"_YEP_dust_cpu_SEIRV_Model_set_data",                (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_set_data,                3},
    {"_YEP_dust_cpu_SEIRV_Model_set_index",               (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_set_index,               2},
    {"_YEP_dust_cpu_SEIRV_Model_set_n_threads",           (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_set_n_threads,           2},
    {"_YEP_dust_cpu_SEIRV_Model_set_rng_state",           (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_set_rng_state,           2},
    {"_YEP_dust_cpu_SEIRV_Model_set_stochastic_schedule", (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_set_stochastic_schedule, 2},
    {"_YEP_dust_cpu_SEIRV_Model_simulate",                (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_simulate,                2},
    {"_YEP_dust_cpu_SEIRV_Model_state",                   (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_state,                   2},
    {"_YEP_dust_cpu_SEIRV_Model_time",                    (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_time,                    1},
    {"_YEP_dust_cpu_SEIRV_Model_update_state",            (DL_FUNC) &_YEP_dust_cpu_SEIRV_Model_update_state,            7},
    {NULL, NULL, 0}
};
}

extern "C" attribute_visible void R_init_YEP(DllInfo* dll){
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
}
