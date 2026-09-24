# Input contract (v0.1)

Each case is a pair of files: `case_<integer>_<name>.mat` and `case_<integer>_<name>_meta.json`. The MAT file must be MATLAB v5/v7 compatible and contain `data` (a numeric `N × C` matrix, `N ≥ 2`) and `header` (a `C × 2` cell array whose first column contains the signal names). Column order follows `header`; time is finite and strictly increasing. Convert channels to physical units before exporting.

| Signal | Unit | Requirement |
|---|---|---|
| `t_s` | seconds | All modes; source time axis |
| `Va_V`, `Vb_V`, `Vc_V` | phase voltage, instantaneous V | All modes |
| `Ia_A`, `Ib_A`, `Ic_A` | phase current, instantaneous A | All modes |
| `P_W` | W | All modes |
| `Q_var` | var | Required for PWR; useful elsewhere |
| `f_hz` | Hz | Required for frequency regulation, inertia, and frequency protection |
| `Vpos_pu` | pu | Optional; otherwise derived from three-phase voltage |
| `breaker` | 0/1 | Optional; used for protection trip timing |

The JSON requires `mode`, integer `case_index`, nonempty `case_name`, and `time_basis` (`recorded` or `absolute`). `recorded` estimates offset from a measured event edge and uses `record_start_nominal_s` from the project configuration as fallback. `absolute` preserves the original `t_s` values as the event metadata clock, even when the first sample is not zero; it applies no event-edge correction. `record_pre_s`, `record_post_s`, and `sim_time_s` are optional. Each mode also needs these numeric metadata fields:

| Mode | Event start / duration | Additional fields |
|---|---|---|
| `frt` | `fault_start_s`, `fault_duration_s` | `V_fault_pu`, `P_ref_pu` |
| `pwr` | `step_start_s`, `step_dur_s` | `P_base_pu`, `P_step_pu`, `Q_base_pu`, `Q_step_pu` |
| `freq_reg` | `event_start_s`, `event_hold_s` | `P_ref_pu`, `freq_deviation_hz` |
| `inertia` | `event_start_s`, `event_hold_s` | `P_ref_pu`, `freq_deviation_hz`, `freq_ramp_rate_hz_per_s` |
| `volt_prot` | `event_start_s`, `event_hold_s` | `P_ref_pu`, `V_target_pu` |
| `freq_prot` | `event_start_s`, `event_hold_s` | `P_ref_pu`, `f_target_hz` |
| `scr` | `disturbance_start_s`, `disturbance_duration_s` | `P_ref_pu`, `SCR_value` |

For a second PWR step, add `step2_start_s`, `step2_dur_s`, `P_step2_pu`, and `Q_step2_pu`. `test_type_code=4` identifies a Q step; otherwise PWR measures the P step. The project configuration JSON requires positive `power_base_w`, `voltage_base_peak_v`, `current_base_peak_a`, and `nominal_frequency_hz`; `record_start_nominal_s` is optional. Consult `examples/` for complete synthetic cases.
