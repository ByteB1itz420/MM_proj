# PROJECT_SPEC — Dynamic Thermal Modelling and Energy-Efficient Cooling Optimization of AI Data Centres Using Optimal Control

**Version:** v0.2 (28 Oct 2026). Supersedes v0.1; see the change log in §17.
**Team:** 5 students.
**Submission:** report (`.tex`, PDF, `.eps` figures) + notebook, confirmed due **before 2 Nov 2026**.
**Presentation:** 2 Nov 2026, 13:00–13:30 IST (10 min talk + 5 min Q&A).
**Course brief:** `Mathematical_Modeling_project_instructions (1).pdf`.

**Implementation status**

| Stage | What | Status |
|---|---|---|
| 1 | Thermal model + numerical verification | **Implemented, awaiting team review** (`src/thermal_model.py`, `notebook/01_thermal_model_verification.ipynb`) |
| FS-PI | Fixed-setpoint PI baseline (§7.1; the team calls this "Stage 2") | **Implemented, awaiting team review** (`src/controllers.py`, `notebook/02_fixed_setpoint_pi.ipynb`; 36/36 checks pass) |
| LA-PI | Load-adaptive PI (B2, §7.2) + cooling-electricity accounting + fair B1/B2 comparison (the team's "Stage 3") | **Implemented, awaiting team review** (`notebook/03_load_adaptive_pi_energy.ipynb`; 52/52 checks pass; adds the closed-cycle comparison variant S2c) |
| OC | Best-found perfect-foresight benchmark + certified LP bound (§7.4; the team's "Stage 4") | **Implemented, awaiting team review** (`src/optimal_control.py`, `notebook/04_perfect_foresight_benchmark.ipynb`; 197/197 checks pass; adds the separately labelled intersample-constrained variant OC-ISC5) |
| Sensitivity | B2 margin and B1 setpoint sweeps (§8.5; the team's "Stage 5") | **Implemented, awaiting team review** (`notebook/05_margin_setpoint_sensitivity.ipynb`; 152/152 checks pass; `figures/fig14_sensitivity_energy_headroom.{png,eps}`) |
| rest | Analysis write-up, packaging into one main notebook, report, slides; MPC bonus not planned | Not started |

### Labels

| Tag | Meaning |
|---|---|
| **[STD]** | Standard physical constant or textbook result. |
| **[ILL]** | **Illustrative.** Plausible order of magnitude only; *not* measured, *not* calibrated, *not* a certified limit. Every result that depends on it must say so. |
| **[SRC-n]** | From source *n* in §15, with its verification status. |
| **[TBD]** | Team decision needed (§16). |

**Ground rules**

- No invented data, references, results or energy-saving percentages. Savings are *measured in our model*, never assumed, and "no meaningful saving" is a valid result.
- **Numerical verification is not real-world validation.** The checks in §9 show that our code solves *our equations* correctly. They say nothing about whether the equations describe a real data centre.

---

## 1. The project in one paragraph

An AI data hall converts almost all of its IT electricity into heat, which a cooling system must remove to keep temperatures within limits. We build a **two-node lumped thermal model** (IT equipment + room air) from the first law of thermodynamics. We keep **heat removed** (W_th) separate from the **electricity spent removing it** (W_e); the two are linked by an illustrative temperature-dependent COP.

On the same simulated plant we compare:

1. a **fixed-setpoint PI** controller;
2. a **load-adaptive PI** controller whose setpoint rule comes from our analysis;
3. **one best-found perfect-foresight optimization**: a benchmark we could not run in reality, solved by a local solver.

A small sensitivity study tests how robust the conclusions are. If time remains, a bonus attempt adds **one MPC horizon with a forecast-error test**.

---

## 2. Model predictive control in plain English (bonus stage, later)

**Analogy: driving.** You look ahead and plan your steering for the next stretch of road, carry out only the first moment of the plan, then look again and re-plan.

**MPC for cooling.** Every few minutes the controller:

1. measures the temperatures;
2. predicts the next few hours from the model plus forecasts of IT load and outdoor temperature;
3. finds the cooling schedule that uses the least predicted electricity while keeping predicted temperatures in range;
4. applies **only the first step** of that schedule;
5. repeats.

Re-measuring each time makes MPC a feedback controller: it corrects for wrong forecasts and model error.

**How the others differ**

- **PI controller:** only reacts to the current error and does not look ahead.
- **Perfect-foresight optimization:** plans the whole period at once with exact knowledge of the future. That is impossible in practice, so it serves as a yardstick.

**Honest caveats**

- MPC can only help if looking ahead has value: COP must change over time *and* thermal mass must be able to store "coolth", or constraints must change predictably.
- A longer horizon is **not guaranteed** to perform better.
- The optimizer can fail to solve in time. MPC therefore needs a fallback law and logging (§7.5).

---

## 3. Scope and stage plan

### 3.1 Core

| Stage | Content | Key checks |
|---|---|---|
| **1** | Two-node model, exact ZOH propagation, verification | V0–V9 (§9.1), implemented |
| **2** | Analysis: equilibria, stability and time constants, PI closed-loop stability, Props. 1–5 (§6) | Written in the report; numerical cross-checks |
| **3** | Fixed-setpoint PI (**FS-PI**) and load-adaptive PI (**LA-PI**) | L2 sampling-interval check, L6 PI stability |
| **4** | **One** best-found perfect-foresight optimization (**OC**) + an LP lower bound on the *same* discretized problem | L1, L3, L4, L5, L7 |
| **5** | Small sensitivity study: one parameter at a time, about 5 parameters × 3 values | — |

### 3.2 Bonus (only after Stages 1–5 are complete)

| Stage | Content |
|---|---|
| **6** | MPC with **one** horizon [TBD]: perfect forecast vs noisy forecast. Report energy, violations, solver failures and fallback use. |

### 3.3 Deferred (not in the current plan; revisit only if the core is finished early)

- GPU-level identification on MIT Supercloud data (§11.3).
- Weather-data retrieval (ERA5 / Open-Meteo). Stages 1–5 use a synthetic outdoor temperature.
- H1 (18–22 °C) temperature band.
- Chilled-water storage (third state).
- κ ≠ 0, i.e. IT power that depends on temperature. **κ = 0 in all simulations**; Prop. 4 remains as an analytical remark only.
- Model mismatch in MPC.
- Large parameter sweeps.

### 3.4 Out of scope

- CFD and hot spots.
- Humidity.
- Liquid cooling.
- Chiller part-load curves.
- Economizers / free cooling.
- Electricity prices.
- Workload scheduling.
- Real-facility validation of energy savings: **every savings number will be a model result.**

---

## 4. Notation, variables, parameters, units

### 4.1 Time

| Symbol | Meaning | Unit |
|---|---|---|
| $t$ | time | s |
| $k$, $t_k$ | input-change index and breakpoint | –, s |
| $\Delta t_c$ | controller update interval | s. Candidate 300 s; **not yet shown to be safe**, see L2 |
| $N_h$ | MPC horizon (bonus) | steps |

### 4.2 States, inputs, disturbances

| Symbol | Role | Meaning | Unit |
|---|---|---|---|
| $T_s$ | state | **Lumped IT-equipment temperature** [ILL meaning]. One number standing for the average temperature of all servers, GPUs and heat sinks. It is **not** a measured GPU-die temperature or a component temperature. | °C |
| $T_a$ | state | Well-mixed room/inlet air temperature | °C |
| $Q_c$ | control | **Heat removed** from the room air by the cooling unit | W_th |
| $P_{IT}$ | disturbance | IT electrical power, all dissipated as heat | W_e |
| $T_{out}$ | disturbance | Outdoor air temperature | °C |

### 4.3 Energy quantities and units

**Unit conversion:**

- Internally, all energies are in joules (J).
- For reporting, $1\ \text{kWh} = 1000\ \text{W} \times 3600\ \text{s} = 3.6\times10^6\ \text{J}$.
- Subscripts mark the meaning, not the unit: **kWh_th** is heat, **kWh_e** is electricity. The two are linked only through COP and must never be added together.

| Symbol | Meaning | Unit |
|---|---|---|
| $E_Q = \int Q_c\,dt$ | heat removed by cooling | J_th (kWh_th) |
| $E_{IT} = \int P_{IT}\,dt$ | IT electricity = IT heat released | J (kWh) |
| $E_{env} = \int UA_e(T_{out} - T_a)\,dt$ | net heat entering through the envelope | J_th |
| $U = C_sT_s + C_aT_a$ | stored thermal energy relative to 0 °C | J |
| $\mathrm{COP}$, $g = 1/\mathrm{COP}$ | heat removed per unit electricity, and its inverse | – |
| $P_{cool} = Q_c/\mathrm{COP}$ | cooling **electrical** power (constant fan power omitted) | W_e |
| $E_{cool} = \int P_{cool}\,dt$ | cooling electricity, the objective (evaluated on *actual* trajectories, §7.4) | J_e (kWh_e) |
| $DH_a$, $DH_s$ | degree-hours outside the constraints | K·h |

### 4.4 Parameters

| Symbol | Meaning | Value | Unit | Tag |
|---|---|---|---|---|
| $\rho$ | air density | 1.2 | kg/m³ | [STD] |
| $c_p$ | air specific heat | 1005 | J/(kg·K) | [STD] |
| $V$ | hall air volume | 3000 | m³ | [ILL] |
| $C_a = \rho c_p V$ | air heat capacity | 3.618 × 10⁶ | J/K | derived |
| $C_s$ | lumped IT heat capacity | 2.5 × 10⁷ | J/K | [ILL] |
| $UA_s$ | IT→air conductance, $= 10^6/15$ (15 K rise at 1 MW) | 6.667 × 10⁴ | W/K | [ILL] |
| $UA_e$ | envelope conductance | 2 × 10³ | W/K | [ILL] |
| $Q_{\max}$ | cooling capacity | 1.5 × 10⁶ | W_th | [ILL] |
| $\eta_{II}$ | fraction of Carnot COP | 0.45 | – | [ILL] |
| $\Delta T_{ev}$, $\Delta T_{cd}$ | evaporator / condenser approach | 12, 8 | K | [ILL] |
| $\mathrm{COP}_{cap}$ | upper edge of COP-model domain (§5.5) | 10 | – | [ILL] |
| $T_{a,\min}$, $T_{a,\max}$ | inlet-air band | 18, 27 | °C | [SRC-5], secondary sources only |
| $T_{s,\max}$ | **illustrative modelling limit** on $T_s$ | 40 | °C | [ILL]. Not a manufacturer specification and not a certified safety limit |
| $\kappa$ | IT power sensitivity to temperature | **0** | 1/K | fixed (§3.3) |

---

## 5. Model derivation

### 5.1 Picture

```
     P_IT (electricity → heat)
           │
   ┌───────▼───────┐  UA_s (T_s − T_a)  ┌───────────────┐  UA_e (T_out − T_a)
   │ IT equipment  │ ─────────────────▶ │   Room air    │ ◀──────────────── outdoors
   │   T_s, C_s    │                    │   T_a, C_a    │
   └───────────────┘                    └───────┬───────┘
                                                │ Q_c (heat removed, W_th)
                                        ┌───────▼───────┐
                                        │ Cooling unit  │ electricity P_cool = Q_c / COP (W_e)
                                        └───────────────┘
```

### 5.2 Assumptions

| # | Assumption | What it hides |
|---|---|---|
| A1 | Each node is lumped (one uniform temperature) | Hot spots, recirculation, rack-to-rack spread |
| A2 | Constant $\rho, c_p, C_s, C_a, UA_s, UA_e$ | Variable airflow and fan speeds |
| A3 | All IT electricity becomes heat in the IT node | Negligible energy carried away in signals |
| A4 | IT→air transfer is $UA_s(T_s - T_a)$ | Forced-flow detail ($\dot m c_p \Delta T$) at fixed airflow |
| A5 | Cooling removes $Q_c \in [0, Q_{\max}]$ instantly from the air node | Compressor dynamics, ramp limits |
| A6 | COP is a Carnot fraction with fixed approach temperatures, valid only in its domain (§5.5) | Real chiller curves, free cooling |
| A7 | Sensible heat only | Humidity control |
| A8 | Envelope exchange $UA_e(T_{out} - T_a)$ | Solar gains, ground coupling |
| A9 | No other internal gains | UPS, lighting |
| A10 | Constant cooling-fan electricity, omitted from comparisons | Fan affinity laws |
| A11 | Controllers measure $T_s$, $T_a$ exactly | Sensor noise, state estimation |
| A12 | $P_{IT}$ is exogenous; κ = 0 | Workload shifting, temperature-dependent server power |

### 5.3 Derivation from the first law

For each lumped node, the first law of thermodynamics with constant heat capacity $C$ reads $C\,\dot T = \sum(\text{heat in}) - \sum(\text{heat out})$ [STD].

**Equation (1): IT node.** Electricity in at $P_{IT}$ is fully dissipated (A3); heat leaves to the air by convection (A4):

$$C_s\,\frac{dT_s}{dt} = P_{IT}(t) - UA_s\,(T_s - T_a). \tag{1}$$

*In plain English:* the equipment warms up when it makes heat faster than the air carries it away. The warmer it is relative to the air, the faster the air carries heat off.

**Equation (2): air node.** The air gains heat from the equipment and from outdoors, and loses heat to the cooling unit. The cooling unit returns air at $T_a$ and supplies it colder at $T_{sup}$, with mass flow $\dot m$. This defines $Q_c \equiv \dot m\,c_p\,(T_a - T_{sup}) \ge 0$, which we treat directly as the control (A5).

$$C_a\,\frac{dT_a}{dt} = UA_s\,(T_s - T_a) + UA_e\,(T_{out} - T_a) - Q_c(t). \tag{2}$$

*In plain English:* the room air warms when it receives more heat (from the servers and through the walls) than the cooling unit removes. Air has little heat capacity, so it reacts quickly.

**Matrix form** ($x = [T_s, T_a]^\top$, $u = Q_c$, $d = [P_{IT}, T_{out}]^\top$):

$$\dot x = Ax + Bu + Ed,\quad
A = \begin{bmatrix} -\frac{UA_s}{C_s} & \frac{UA_s}{C_s} \\ \frac{UA_s}{C_a} & -\frac{UA_s + UA_e}{C_a} \end{bmatrix},\;
B = \begin{bmatrix} 0 \\ -\frac{1}{C_a} \end{bmatrix},\;
E = \begin{bmatrix} \frac{1}{C_s} & 0 \\ 0 & \frac{UA_e}{C_a} \end{bmatrix}. \tag{3}$$

The units of every row are K/s.

**Energy identity.** Add (1) and (2); every term is in watts:

$$\frac{dU}{dt} = P_{IT} + UA_e\,(T_{out} - T_a) - Q_c. \tag{4}$$

The internal transfer between the nodes cancels.

### 5.4 Heat removed vs electricity used

$Q_c$ is heat moved out of the room. Eq. (4) fixes how much heat must be removed over any period. $P_{cool}$ is the electricity used to pump that heat to the warmer outdoors:

$$P_{cool}(t) = \frac{Q_c(t)}{\mathrm{COP}(t)}. \tag{5}$$

### 5.5 COP model and its valid domain

$$\mathrm{COP} = \eta_{II}\,\frac{T_{ev}}{T_{cd} - T_{ev}},\qquad
T_{ev} = T_a - \Delta T_{ev} + 273.15\ \text{K},\quad
T_{cd} = T_{out} + \Delta T_{cd} + 273.15\ \text{K}. \tag{6}$$

**Domain.**

- The formula is meaningful only for a positive "lift", $T_{cd} - T_{ev} > 0$.
- We further restrict it to $\mathrm{COP} \le \mathrm{COP}_{cap}$, i.e. lift $\ge \eta_{II}T_{ev}/\mathrm{COP}_{cap}$, which is about 12.6–13.0 K with these parameters.
- With the defaults, this requires $T_{out} \ge T_{out,\min}(T_a)$: about **10.6 °C at $T_a$ = 18 °C** and **20.0 °C at $T_a$ = 27 °C**.
- Below that, a real plant would use free cooling, which we do not model. The code therefore **raises an error instead of silently capping**, and synthetic $T_{out}$ profiles must stay above 20 °C.

Within the domain, COP increases with $T_a$ and decreases with $T_{out}$. The values are illustrative and not calibrated to any machine. (Moore et al. 2005 [SRC-7] report an empirical COP that rises with supply temperature; we could not verify its coefficients and **do not use them**.)

---

## 6. Where can energy savings genuinely come from?

Integrating (4) over $[0, t_f]$ gives

$$E_Q = E_{IT} + E_{env} - \Delta U, \qquad E_{cool} = \int_0^{t_f} \frac{Q_c(t)}{\mathrm{COP}(T_a(t), T_{out}(t))}\,dt. \tag{7}$$

1. **Prop. 1 (no free lunch).** With constant COP, $UA_e = 0$ and the same initial and terminal state, $E_{cool} = E_{IT}/\mathrm{COP}$ for **every** feasible controller. Clever timing saves nothing, so any saving the code reports in this case is a bug.
2. **Prop. 2 (envelope).** With constant COP, $E_Q$ decreases as $\int T_a\,dt$ increases, so a warmer room is cheaper.
3. **Prop. 3 (approximate timing bound).** Shifting heat removal to high-COP hours is limited by stored coolth, about $C_{tot}\,\Delta T_{band}\,(1/\mathrm{COP}_{low} - 1/\mathrm{COP}_{high})$ per cycle. This is an order-of-magnitude bound, not a theorem.
   *[ILL] back-of-envelope:* about 71 kWh_th storable vs about 19 MWh_th of daily heat. **Hypothesis:** timing leverage is small.
4. **Prop. 4 (break-even κ; analytical remark only, κ = 0 in simulations).** Raising $\bar T$ saves energy only if $\kappa < -g'(\bar T)/(1 + g(\bar T))$, which is about 0.7 %/K at the [ILL] point $T_a$ = 22 °C, $T_{out}$ = 32 °C. Savings are not guaranteed in reality.
5. **Prop. 5 (load-adaptive setpoint).** At quasi-steady state, $T_s - T_a = P_{IT}/UA_s$, so the warmest admissible air temperature is $\bar T_a^*(t) = \min(T_{a,\max},\ T_{s,\max} - P_{IT}(t)/UA_s)$. A fixed setpoint must be safe at peak load. *How much* this matters is computed, not assumed.

### Reporting decomposition

$$\underbrace{E_{FS} - E_{LA}}_{\text{adaptive-setpoint rule}} \;+\; \underbrace{E_{LA} - E_{OC}}_{\text{remaining value of optimization}} \;=\; E_{FS} - E_{OC}.$$

The effect of the setpoint *value* is studied separately in the sensitivity stage (FS-PI at several fixed setpoints). That way it is not credited to optimization.

---

## 7. Controllers (Stages 3–4; bonus Stage 6)

Within a comparison, all runs share the same plant, disturbance trace, horizon, actuator limits, $\Delta t_c$ and constraints.

**As implemented (Stages 3–5), initial states are *not* common.**

- Each controller starts at its **own** equilibrium (its own setpoint at the starting load) and must return to that state by the end of the common 17 h horizon (a closed cycle with terminal checks). The initial states differ between controllers and are reported with every result.
- OC is solved once per boundary protocol (B1's or B2's start/end state) and compared only with that controller.

### 7.1 FS-PI: fixed-setpoint PI

- Law: $Q_{c,k} = \mathrm{sat}_{[0,Q_{\max}]}\big(Q_{ff} + K_p e_k + K_i\Delta t_c\sum_j e_j\big)$, with $e_k = T_{a,k} - T_{set}$ and conditional-integration anti-windup.
- $T_{set} = \min(T_{a,\max},\ T_{s,\max} - P_{IT,\max}/UA_s) - m$.
- **As implemented, the gains are stability-tuned, not energy-tuned.** $K_p$ and $K_i$ were chosen once per update interval by the Stage 2 rule: minimum unsaturated closed-loop spectral radius on a grid, subject to stability at 2× the gains. They are stored as `STAGE2_GAINS` and reused unchanged by B1 and B2. No separate tuning day was used.
- The setpoint is 24 °C, i.e. the rule with $P_{IT,\max}=1$ MW and $m=1$ K (provisional). Stage 5 shows the effect of 23 and 25 °C.

### 7.2 LA-PI: load-adaptive PI

- The same PI law and the same stored gains track $\mathrm{clip}(\bar T_a^*(t) - m,\ 18,\ 27)$ (Prop. 5), with the setpoint recomputed from the measured load at every update.
- $m=1$ K is provisional; Stage 5 sweeps $m\in\{0.5,1,1.5,2\}$ K.

### 7.3 Sampling caution

Exact ZOH propagation says nothing about what happens *between* controller updates, or about whether a 5-minute update is safe. Stage 1 already shows that $T_a$ can move by several kelvin within about one fast time constant (≈ 46 s) after a load step. **L2** must check inter-sample violations on the continuous plant for $\Delta t_c \in \{60, 300\}$ s before we adopt 300 s.

### 7.4 OC: one best-found perfect-foresight optimization

**Discretized problem $\mathcal P_d$.** Exact ZOH dynamics at $\Delta t_c$; bounds $0 \le u_k \le Q_{\max}$; band constraints at sample points; the same initial state $x_0$; terminal condition $x_N = x_0$; and objective $J_d = \sum_k \Delta t_c\,u_k\,g(T_{a,k}, T_{out,k})$.

**Solver.** IPOPT via CasADi 3.8.2 (exact sparse Hessian, MUMPS), warm-started from the matched PI sequence (or, if that is infeasible, from the nearest feasible point) plus three seeded perturbed starts. The LP uses HiGHS through `scipy.optimize.linprog`; every configuration is logged, and the bound is certified by weak duality from the solver's multipliers (Stage 4 found that SciPy's default HiGHS path can fail on these LPs).

**Stage 4 additions.**

- **Boundary protocols.** OC is solved twice per scenario and $\Delta t_c$: once with B1's start/end state (compared only with B1) and once with B2's (compared only with B2). §10.1's "same initial state" is met within each pair. The Stage 3 own-equilibrium starts replace the recovery protocol, because every run here ends where it started.
- **Variant OC-ISC5.** It also imposes the limits every 5 s inside each interval, because sample-point feasibility alone allowed exceedances between samples in S1.

**Claim discipline:**

- The result is a **best-found** solution of $\mathcal P_d$, **not a proven global optimum**.
- Multi-start agreement is evidence, not proof.

**LP lower bound (valid for $\mathcal P_d$ only).** Every feasible point of $\mathcal P_d$ satisfies $u_k \ge 0$ and $T_{a,k} \le T_{a,\max}$. Because $g$ is non-increasing in $T_a$ on its domain, $u_k\,g(T_{a,k}, T_{out,k}) \ge u_k\,g(T_{a,\max}, T_{out,k})$.

- Replacing $g$ by this constant gives an LP with **identical** dynamics, constraints, $\Delta t_c$, initial and terminal conditions and quadrature rule. Its optimum satisfies $J_{LB} \le J_d^* \le J_d(\text{best})$.
- The relative gap $(J_d(\text{best}) - J_{LB})/J_{LB}$ is therefore an *upper bound on suboptimality within $\mathcal P_d$*. It does **not** bound continuous-time or inter-sample effects.
- The bound requires $(T_{a,\max}, T_{out,k})$ to lie inside the COP domain.
- **L5** asserts that the LP and NLP are built from the same matrices and constraint lists.

**Reported energy.** For every controller, including OC, the reported $E_{cool}$ is computed by **simulating the plant** under the applied input sequence and integrating $Q_c(t)/\mathrm{COP}(T_a(t), T_{out}(t))$ along the **actual** state trajectory. This is not the optimizer's objective value. **L1** checks this quadrature by refinement. (Exact state propagation does **not** make this integral exact, because COP is nonlinear in $T_a$.)

### 7.5 MPC (bonus Stage 6)

**Setup.**

- One horizon $N_h$ [TBD].
- Two forecast cases: perfect, and noisy (AR(1) error on $T_{out}$ and on $P_{IT}$, with [ILL] magnitudes).
- Soft state constraints with slack penalties.

**Soft constraints do not guarantee success.** They make $\mathcal P_d$ feasible for the solver's model, but the solver can still fail to converge, hit an iteration or time limit, or return a poor point. The plant can still violate limits, because of forecast error or actuator saturation.

**Fallback.** If the solver status is not "solved" or the time budget is exceeded:

1. apply the next step of the last successful plan, if one exists;
2. otherwise apply the FS-PI law.

**Logging, per step.** Solver status, solve time, slack values, whether the fallback was used, and plant constraint violations. Totals are reported.

**No monotonicity claim.** We do not require or claim that a longer horizon performs better. Only one horizon is planned.

---

## 8. Analysis (Stage 2)

### 8.1 Equilibria

$T_a^* = T_{out} + (P - Q_c)/UA_e$ and $T_s^* = T_a^* + P/UA_s$.

The cooling that holds $T_a$ at $\bar T$ is $Q_c = P + UA_e(T_{out} - \bar T)$. If $UA_e = 0$, there is no equilibrium unless $Q_c = P$ exactly.

### 8.2 Stability and time constants

- $\mathrm{tr}A < 0$ and $\det A = UA_s UA_e/(C_sC_a) > 0$.
- $A$ is similar to a symmetric matrix via $D = \mathrm{diag}(\sqrt{C_s}, \sqrt{C_a})$, so its eigenvalues are real and negative.
- With [ILL] parameters, verified in V3:

| Mode | Eigenvalue | Time constant |
|---|---|---|
| Fast (air–equipment exchange) | $\lambda_{fast} \approx -2.158\times10^{-2}$ s⁻¹ | **τ ≈ 46.3 s** |
| Slow (envelope / total heat storage) | $\lambda_{slow} \approx -6.83\times10^{-5}$ s⁻¹ | **τ ≈ 4.07 h** |

- The time scales differ by a factor of about 316. Forward Euler would be stable only for $h < 2/|\lambda_{fast}| \approx 93$ s. We therefore use exact ZOH.
- If $UA_e = 0$, one eigenvalue is zero: the hall is only marginally stable without control.

### 8.3 PI closed-loop stability

This analysis is continuous-time, linear and unsaturated. With integrator $\dot z = T_a - T_{set}$, the characteristic polynomial is $s^3 + c_2 s^2 + c_1 s + c_0$, where

$$c_2 = \tfrac{UA_s}{C_s} + \tfrac{UA_s + UA_e + K_p}{C_a},\qquad
c_1 = \tfrac{K_i}{C_a} + \tfrac{UA_s}{C_s}\,\tfrac{UA_e + K_p}{C_a},\qquad
c_0 = \tfrac{UA_s}{C_s}\,\tfrac{K_i}{C_a}.$$

Routh–Hurwitz ($c_i > 0$ and $c_2c_1 > c_0$) holds for all $K_i > 0$, $K_p > -UA_e$.

This does **not** cover sampling, saturation or delay. Those are checked with discrete-time eigenvalues and simulation (L6).

### 8.4 Propositions

Props. 1–5 come with proofs; Pontryagin's principle is an optional appendix [SRC-9].

### 8.5 Small sensitivity study (Stage 5)

- **As implemented (notebook 05):** one-at-a-time sweeps of B2's margin $m\in\{0.5,1,1.5,2\}$ K and of B1's fixed setpoint $\in\{23,24,25\}$ °C, on S1 and S2c at both update intervals.
  - Gains are unchanged.
  - Each run closes on its own equilibrium.
  - Safety is required on both the 5 s and the 1 s grid.
- **Not done (scope cut):** sensitivity to the physical parameters $\eta_{II}$, $C_s$, $UA_e$ and $UA_s$. Every conclusion remains conditional on the [ILL] parameter values.

---

## 9. Numerical checks

### 9.1 Stage 1 checks (implemented; outputs are in the notebook)

| ID | Check | Criterion |
|---|---|---|
| V0 | Matrix form (3) reproduces the term-by-term RHS (1)–(2); J↔kWh conversion | Agreement to round-off |
| V1 | Zero load, $Q_c = 0$, $T_{out} = T_s = T_a$: nothing changes | Drift < 1e-12 K over 24 h. Accumulators ≈ 0 relative to the size of the cancelling terms $UA_eT_0t_f$ (< 1e-9) |
| V2 | Constant load: analytic equilibrium = linear solve; it is preserved; decay from a perturbed state | Residuals tiny; converged after 30 τ_slow |
| V3 | Eigenvalues real and negative; τ_fast ≈ 46 s, τ_slow ≈ 4.07 h; agree with trace/determinant formula and with `eig(A)` | Within stated tolerances |
| V4 | Exact ZOH vs `solve_ivp` (Radau, rtol 1e-10), **restarted at every input change**, on an irregular piecewise-constant input sequence | Max \|ΔT\| < 1e-6 K; energies agree |
| V5 | Energy balance (4) using **consistently integrated** accumulators ($E_Q$, $E_{IT}$, $E_{env}$ propagated with the state) | Residual / $E_{IT}$ < 1e-9 |
| V6 | Contrast: $E_{env}$ by trapezoid on *sampled* states is **not** exact | Error reported (not a pass/fail) |
| V7 | COP domain: monotone in its domain; raises outside it; $T_{out,\min}$ values | As in §5.5 |
| V8 | Illustrative transient (§4 of notebook): ZOH vs `solve_ivp` at every 10 s output point | Max \|ΔT\| < 1e-6 K |
| V9 | Illustrative transient: energy balance (4) | Residual / $E_{IT}$ < 1e-9 |
| V10 | Invalid-input regression tests. `cop()` raises `COPDomainError` for non-finite temperatures, $T_{ev} \le 0$ K, and zero or negative lift. `simulate_zoh()` and `simulate_ivp()` raise `ValueError` for NaN/inf values, backwards, non-monotonic or repeated time points, wrong shapes, and fewer than 2 time points | Exactly the stated exception is raised (27 cases) |

**Result (28 Oct 2026, fresh kernel): 50 of 50 checks pass.** The measured values are printed in the notebook.

**What these checks do *not* establish:**

- that the electrical-energy integral is exact (L1);
- that a 5-minute controller update is safe (L2);
- anything about a real facility.

### 9.2 Later checks (planned)

| ID | Stage | Check |
|---|---|---|
| L1 | 3–4 | Quadrature convergence of $E_{cool}$ on actual trajectories (sub-sampling refinement). **Done** (notebooks 03, 04: 5 s vs 1 s, plus an independent 0.5 s trapezoid) |
| L2 | 3 | Inter-sample constraint violations vs $\Delta t_c$ (60 vs 300 s) on the continuous plant |
| L3 | 4 | OC objective as $\Delta t_c$ is refined; multi-start spread. **Done** (notebook 04: 300 vs 60 s; spread reported per problem) |
| L4 | 3–4 | Recovery protocol (§10) achieves matched terminal states within tolerance. **Replaced** by closed-cycle runs from each controller's own equilibrium, with terminal checks (notebooks 03, 04) |
| L5 | 4 | LP bound and NLP use identical discretized problem data; $J_{LB} \le J_d(\text{best})$. **Done** (notebook 04: same array objects, NLP Jacobian = A, certified dual bound) |
| L6 | 3 | PI stability: Routh–Hurwitz (continuous) and discrete eigenvalues vs simulation |
| L7 | 4 | Prop. 1 invariance: constant COP, $UA_e = 0$, matched states → equal $E_{cool}$ across controllers. **Done** (notebook 04, PI and OC at both intervals) |
| L8 | 6 | MPC: solver-failure and fallback counts logged; violations counted on the plant |

---

## 10. Fair-comparison protocol

1. **Initial state (as implemented).** Each run starts at its **own** equilibrium at the initial disturbance values, so initial states differ between controllers and are reported. (The original plan of a common $x_0$ was not used.)
2. **Terminal states (as implemented):** every run is a closed cycle. The disturbances return to their starting values, the horizon extends 8 h past the last change, and each run must return to its own start state within $|\Delta T|\le 10^{-3}$ K and $|\Delta U|\max g\le10^{-4}E_{cool}$. No recovery controller and no $\Delta U$ correction are used.
   *Original plan, not implemented:* **matched terminal states via an explicit recovery protocol.** Runs can end in different states, and dividing $\Delta U$ by a mean COP is **not** an exact correction, because COP varies along the path. Instead:
   - after the evaluation window, every run hands over to the **same recovery controller** (the FS-PI law, integrator reset, common setpoint $T_{rec}$) under the same disturbance continuation, for $t_{rec}$ [TBD];
   - recovery electricity is **counted**, using the actual COP;
   - L4 asserts that terminal states match within $\varepsilon_T$ [TBD, e.g. 0.01 K]. If they do not, extend $t_{rec}$ or report the remaining mismatch as an explicit uncertainty, not as a correction;
   - OC additionally enforces $x_N = x_0$ inside $\mathcal P_d$, and is still evaluated through the same protocol.
3. **Same plant for all.** The continuous plant (exact ZOH plus fine evaluation), never a controller's internal model. Inter-sample violations count.
4. **Equal tuning effort (as implemented):** B1 and B2 use the identical stability-designed Stage 2 gains. Neither is tuned for energy, and there was no separate tuning day.
5. **Energy is always reported together with violations:** $E_{cool}$, $E_Q$, $DH_a$, $DH_s$, maximum violation, peak $P_{cool}$, solve time and failure counts.
6. **Several test days** (planned, **not done**): only the two synthetic closed-cycle scenarios S1 and S2c were run, so there is no mean ± s.d.
7. **Metrics are fixed before results are seen.** Null and negative results are reported.
8. **OC is labelled "best-found perfect-foresight benchmark (not realizable)".**

---

## 11. Experiments

### 11.1 Synthetic only (Stages 1–6)

| ID | Purpose | Inputs |
|---|---|---|
| S0 | Verification (Stage 1; L7) | Piecewise-constant test inputs; constant COP for L7 |
| S1 | FS-PI vs LA-PI vs OC (as implemented) | Step trace: load 0.6 → 1.0 MW at 1 h → 0.6 MW at 9 h, $T_{out}$ 30 °C, 17 h [ILL]. No diurnal profile was used |
| S2c | The same comparisons | 0.8 MW; $T_{out}$ 30 → 35 °C at 1 h → 30 °C at 9 h, 17 h [ILL] |
| S2 | Sensitivity (§8.5) | S1 and S2c traces |
| S3 | MPC bonus | S1 inputs + forecast noise |

**All energy figures are model results with [ILL] parameters.** They are not evidence about real data centres.

### 11.2 Verification vs validation

| | What it is | Where |
|---|---|---|
| **Verification** | The code solves the stated equations correctly | §9 |
| **Validation** | The equations describe reality | **Not attempted** in the current scope |

The report will state this distinction explicitly.

### 11.3 Deferred real-data items (kept for reference, not planned)

**GPU-level RC identification** on MIT Supercloud V1 telemetry [SRC-1, SRC-2]:

- **Already checked:** columns include `temperature_gpu` (whole °C) and `power_draw_W`, at about 100 ms; files are 0.5–50 MB; licence CC BY-NC-ND 4.0, so no redistribution.
- **Planned protocol:** split by job, free-run tests on held-out jobs, MAE/RMSE vs persistence/static baselines.
- **What it would validate:** only the lumped-capacitance *structure* at device level. It would not validate $T_s$ as defined here, the room node, or any saving.

**ERA5 outdoor temperature** via Open-Meteo [SRC-3, SRC-4]. This would make the inputs realistic; it would not validate the model.

---

## 12. Rubric mapping

| Rubric item (marks) | Delivery |
|---|---|
| **LaTeX report (3)** | Intro/motivation → objectives → assumptions (§5.2) → derivation (§5) + schematic → analysis (§6, §8) → numerical methods and verification (§9) → controllers (§7) → results (§10–11) → limitations (incl. verification ≠ validation) → AI acknowledgement → references. Figures saved as `.eps` without transparency. |
| **Notebook (2)** | Executable; coloured Markdown; equations and algorithms before code; labelled plots; asserts. Notebook 01 (Stage 1) exists; later stages add notebooks or sections. |
| **Presentation (5)** | One block per member (§13.2); rehearsed Q&A |
| **Bonus (≤5)** | Optimality gap within $\mathcal P_d$ (LP bound); MPC forecast-error test with failure logging; verification metrics. Real-data validation is **deferred**, so the bonus case is weaker without it; revisit only if the core is done (D8). |
| **Integrity** | AI use acknowledged truthfully (`AI_USE_LOG.md`, notebook header, report section); every external source cited |

---

## 13. Deliverables, roles, timeline

### 13.1 Files

```
PROJECT_SPEC.md
AI_USE_LOG.md
requirements.txt
src/thermal_model.py                               Stage 1 (done)
notebook/01_thermal_model_verification.ipynb       Stage 1 (done)
figures/fig01_transient_illustrative.{eps,png}     Stage 1 (done)
src/controllers.py                                 FS-PI, LA-PI (done)
notebook/02_fixed_setpoint_pi.ipynb                FS-PI (done)
notebook/03_load_adaptive_pi_energy.ipynb          LA-PI + comparison (done)
src/optimal_control.py                             OC + LP bound (done)
notebook/04_perfect_foresight_benchmark.ipynb      OC (done)
figures/fig02-fig13 (.png)                         Stages 2-4 (done; .eps for the report later)
report/  slides/                                   later
```

### 13.2 Roles [D9]

| Member | Owns |
|---|---|
| 1 | Physics, derivation, COP domain |
| 2 | Analysis (stability, propositions) |
| 3 | Numerics & verification (Stage 1, L-checks) |
| 4 | Controllers + OC + LP bound |
| 5 | Sensitivity study; MPC bonus; report assembly lead |

### 13.3 Timeline (submission before 2 Nov)

| Date | Milestone |
|---|---|
| Wed 28 Oct | Stage 1 implemented; **team reviews it** |
| Thu 29 Oct | Stage 2 write-up; Stage 3 FS-PI / LA-PI + L2, L6 |
| Fri 30 Oct | Stage 4 OC + LP bound + recovery protocol; L1, L3–L5, L7 |
| Sat 31 Oct | Stage 5 sensitivity; figures; full report draft; Stage 6 MPC only if Stages 1–5 are done |
| Sun 1 Nov | Final report + notebook **submitted**; slides; rehearsals |
| Mon 2 Nov | Present 13:00 IST |

**Cut order if behind:**

1. MPC (Stage 6);
2. multi-start depth;
3. sensitivity reduced to 3 parameters;
4. OC replaced by the LP (stated honestly).

---

## 14. Risks

| Risk | Mitigation |
|---|---|
| NLP local minima / solver failure | Best-found wording, multi-start, LP bound on the same $\mathcal P_d$ |
| 300 s update unsafe | L2; fall back to 60 s |
| Terminal-state mismatch biases comparisons | Recovery protocol with an explicit tolerance (L4) |
| Over-claiming | "Model result, illustrative parameters" on every number; verification ≠ validation |

---

## 15. Sources and verification status (checked 28 Oct 2026)

| # | Source | Used now? | Status |
|---|---|---|---|
| SRC-1 | S. Samsi *et al.*, "The MIT Supercloud Dataset," *IEEE HPEC 2021*, pp. 1–8, doi:10.1109/HPEC49654.2021.9622850; arXiv:2108.02037 | Deferred | Verified (paper text; README citation) |
| SRC-2 | MIT Supercloud V1, `s3://mit-supercloud-dataset/datacenter-challenge/202201/` | Deferred | Verified (listing, CSV header, CC BY-NC-ND 4.0) |
| SRC-3 | P. Zippenfenig, Open-Meteo.com Weather API (2023), doi:10.5281/ZENODO.7970649 | Deferred | Verified (API and licence pages; CC BY 4.0) |
| SRC-4 | H. Hersbach *et al.*, ERA5 hourly data on single levels (2023), doi:10.24381/cds.adbb2d47 | Deferred | Citation as listed by Open-Meteo; DOI page not opened |
| SRC-5 | ASHRAE TC 9.9, *Thermal Guidelines for Data Processing Environments*, 5th ed. (2021) | Yes (18–27 °C) | Secondary sources only; primary not accessed. Write "as reported in …" |
| SRC-6 | The Green Grid WP#49, PUE | Context | Exists; full text not read |
| SRC-7 | Moore, Chase, Ranganathan, Sharma (2005), "Making Scheduling 'Cool' …" | Context only | Existence via a citing article only; coefficients **not used** |
| SRC-8 | Rawlings, Mayne, Diehl, *Model Predictive Control*, 2nd ed., Nob Hill, 2017 | MPC background | Bibliographic details verified |
| SRC-9 | D. E. Kirk, *Optimal Control Theory: An Introduction*, Dover, 2004, ISBN 9780486434841 | Optional appendix | Bibliographic details verified |
| — | A thermodynamics textbook the team has actually opened | First law, Carnot COP | Team to choose and verify |

No source enters `refs.bib` unless a team member has opened it and checked the specific claim.

---

## 16. Decisions for the team

| ID | Decision | Recommendation |
|---|---|---|
| D1 | OC solver: CasADi/IPOPT vs SciPy only | **Used in Stage 4:** CasADi 3.8.2/IPOPT; LP via HiGHS in `scipy.optimize.linprog`, with a dual certificate. Team to confirm |
| D2 | $\Delta t_c$ | Partly answered by notebook 02 (FS-PI only). The 60 s gains are unstable at 300 s; the retuned 300 s loop exceeds 27 °C briefly in the S1 load step, while 60 s does not. Leaning towards 60 s |
| D3 | FS-PI setpoint rule and margin $m$ | As §7.1; $m=1$ K provisional, not tuned; the effect is shown by the Stage 5 sweep |
| D4 | Synthetic $P_{IT}$ and $T_{out}$ profiles | Diurnal sinusoids + noise; $T_{out}$ kept ≥ 20 °C (COP domain) |
| D5 | Recovery protocol: $T_{rec}$, $t_{rec}$, $\varepsilon_T$ | **Superseded:** closed cycles from each controller's own equilibrium, with terminal checks (§10) |
| D6 | MPC horizon (bonus) | One value, e.g. 6 h |
| D7 | Sensitivity parameters | $\eta_{II}$, $C_s$, $UA_e$, $UA_s$, FS setpoint |
| D8 | Revisit deferred real-data items after the core? | Only if Stages 1–5 are done by Sat noon |
| D9 | Roles (§13.2) | Assign now |
| D10 | Review the [ILL] parameter values | Especially $C_s$, $UA_s$, $T_{s,\max}$, $\eta_{II}$ |

---

## 17. Change log

**v0.2.2 (30 Oct 2026) — Stage 5 and documentation corrections**

- §8.5: records the sensitivity actually run (B2 margin and B1 setpoint), and that physical-parameter sensitivity was not done.
- §7, §7.1, §7.2, §10, §11.1, D3, D5 now describe what was actually implemented:
  - own-equilibrium closed cycles (not a common initial state and not a recovery protocol);
  - stability-tuned Stage 2 gains (not energy-tuned, no tuning day);
  - step traces S1 and S2c (not diurnal profiles);
  - two scenarios only.
- Notebook 04: two text cells corrected without re-running; its code and outputs are unchanged.
  - OC-ISC5's extra cost refers to J_d. Its replayed electricity is lower than base OC's at 60 s.
  - The 1 s-grid exceedances now sit next to the S1 OC-ISC5 results, which are stated to be not continuously safe.

**v0.2.1 (30 Oct 2026) — Stage 4 implementation notes**

- §7.4: solver choice recorded; boundary protocols; the OC-ISC5 variant; LP bound certified by weak duality.
- §9.2: L1, L3, L5 and L7 done; L4 replaced (closed cycles from own equilibria, with terminal checks).
- `controllers.simulate_closed_loop` gained an optional `inv_cop` argument, used only by the invariance test. Its default is unchanged, and notebooks 01–03 re-run with identical output.

**v0.2 (28 Oct 2026)**

- Removed the deadline alarm; submission is confirmed before 2 Nov.
- Reduced scope: FS-PI, LA-PI, one best-found OC and a small sensitivity study. MPC is now a later bonus with one horizon and a forecast-error test. Deferred: GPU identification, weather retrieval, H1, storage, κ ≠ 0, large sweeps, model mismatch. Dropped B0 as a separate controller; the effect of the setpoint value moves to the sensitivity study.
- **Correction:** the NLP result is "best-found", not a proven optimum. The LP bound must use the identical discretized problem, and energy is evaluated on actual trajectories with the actual COP (§7.4).
- **Correction:** removed the "ΔU / mean COP" terminal correction. Comparisons now use the same initial state plus an explicit recovery protocol with counted recovery energy and matched terminal states (§10).
- **Correction:** removed the claim that MPC improves monotonically with horizon, and the claim that soft constraints guarantee feasibility. Added fallback control and failure/violation logging (§7.5).
- **Correction:** $T_s$ and the 40 °C limit are labelled illustrative (not GPU-die temperature, not a certified limit). Added the COP domain and replaced the silent cap with an error (§5.5). Added "verification ≠ validation" (ground rules, §11.2).
- Updated τ_slow to 4.07 h (v0.1 said "≈ 4 h").
- Added Stage 1 implementation and verification checks V0–V9, and the later checks L1–L8.

**v0.1 (28 Oct 2026)** — initial draft.
