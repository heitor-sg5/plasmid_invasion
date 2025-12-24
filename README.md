# Plasmid Transfer Dynamics Models

## Overview

Plasmid invasion dynamics were studied using two modelling methods: a deterministic system of ordinary differential equations (ODEs) and a stochastic simulation algorithm (SSA). Both approaches track the population dynamics of plasmid-free bacteria ($F$) and plasmid-bearing bacteria ($P$), incorporating population growth, density dependence, plasmid transfer via conjugation, segregational loss, and plasmid-associated fitness costs. The deterministic model captures average population behaviour, while the stochastic model resolves demographic noise and rare-event dynamics that arise in near-zero populations.

---

# Table of Contents

- [Project Structure](#project-structure)
- [Introduction](#introduction)
  - [Horizontal Gene Transfer](#horizontal-gene-transfer)
  - [Cost and Transfer Fitness Trade-off](#cost-and-transfer-fitness-trade-off)
  - [Selection and Invasion Thresholds](#selection-and-invasion-thresholds)
  - [Study Aims and Objectives](#study-aims-and-objectives)
- [Methods](#methods)
  - [Deterministic ODE Model](#deterministic-ode-model)
  - [Equilibria and Invasion Analysis](#equilibria-and-invasion-analysis)
  - [Stochastic SSA Model](#stochastic-ssa-model)
  - [Outputs and Visualizations](#outputs-and-visualizations)
- [Results](#results)
  - [Deterministic Invasion Threshold and Equilibria](#deterministic-invasion-threshold-and-equilibria)
  - [Deterministic Dynamics Across Selective Pressure](#deterministic-dynamics-across-selective-pressure)
  - [Stochastic Effects and SSA Outcomes](#stochastic-effects-and-ssa-outcomes)
- [Conclusions](#conclusions)
- [References](#references)

---

## Project Structure

```
├── main.py
├── ode_model/
│   ├── model.py
│   ├── charts.py
│   └── figures/
├── ssa_model/
│   ├── model.py
│   ├── charts.py
│   └── figures/
├── appendices/
├── requirements.txt
└── README.md
```

---

## Introduction

### Horizontal Gene Transfer

Plasmids are self-replicating, extra-chromosomal DNA molecules that are widely distributed among bacterial populations. In addition to vertical transmission through cell division, many plasmids are capable of horizontal transfer between neighbouring cells, allowing them to spread independently of host reproduction. Through this process, plasmids can carry and accessory genes that enable bacteria to adapt to new environments or stressful conditions, including genes conferring antibiotic resistance, virulence, or specialised metabolic functions. 

Among the mechanisms of horizontal gene transfer, conjugation-direct plasmid transfer from a donor to a recipient cell via cell–cell contact is the best studied and is believed to play a central role in plasmid persistence at the population level. From both ecological and clinical perspectives, conjugative plasmids are of particular importance, as they are major drivers of the spread of antibiotic resistance genes and have contributed substantially to the global antimicrobial resistance crisis. 

The plasmid life cycle is governed by three key processes: intracellular replication and copy-number control, partitioning between daughter cells at division, and, for conjugative plasmids, horizontal transmission between cells. Together, these processes determine whether plasmids are stably maintained or lost from a bacterial population. Because plasmid inheritance is imperfect, plasmid-free cells can arise through segregational loss, and these cells may have a fitness advantage if plasmid carriage is costly. As a result, plasmid persistence is not guaranteed even when plasmids encode beneficial traits [1, 2].

### Cost and Transfer Fitness Trade-off

Although plasmids can provide adaptive benefits to their hosts, plasmid carriage is frequently associated with a fitness cost. These costs are commonly attributed to the metabolic burden of plasmid replication, maintenance, and expression, and can vary widely across plasmid-host combinations. In the absence of positive selection, plasmid-free cells often grow faster than plasmid-bearing cells, allowing plasmids to be out-competed and eliminated from the population. As a result, plasmid persistence is not guaranteed even when plasmids are capable of horizontal transfer. 

Two broad mechanisms have been proposed to explain how plasmids persist over evolutionary timescales. First, plasmids may be maintained through positive selection, where the benefits conferred by plasmid-encoded traits outweigh their associated costs. Second, plasmids may persist as infectious elements, spreading horizontally at a sufficiently high rate to compensate for reduced host growth. In this latter scenario, plasmids behave as genetic parasites whose persistence depends critically on the balance between host fitness costs and the rate of horizontal gene transfer. 

Conjugation dynamics are governed by two key processes: the efficiency of plasmid transfer and the growth dynamics of plasmid-bearing cells. Both processes depend on intrinsic plasmid properties and environmental conditions. Experimental studies have shown that conjugation rates can vary by orders of magnitude depending on host physiology, while plasmid-associated fitness effects can range from strongly deleterious to beneficial. Together, these factors determine whether horizontal transfer can offset plasmid cost and segregational loss. 

In addition to long-term fitness costs measured in stable plasmid-bearing lineages, plasmid acquisition itself can impose an immediate metabolic burden on newly formed transconjugants. Following conjugation, transient disruptions to gene regulation and resource allocation can reduce growth rates or prolong lag phases before exponential growth resumes. Although these short-term acquisition costs are difficult to quantify and are often excluded from traditional fitness measurements, they further reinforce the trade-off between plasmid transfer and host growth [3, 4].

### Selection and Invasion Thresholds

A central question in plasmid ecology and antibiotic resistance management is whether resistance can be eliminated by removing selective pressure, for example through reduced antibiotic use. In principle, if resistance-conferring plasmids impose a fitness cost, plasmid-free cells should outcompete plasmid-bearing cells once selection is relaxed, leading to resistance reversal. However, both experimental and theoretical studies have shown that this outcome is far from guaranteed. 

Several mechanisms can allow plasmids to persist in the absence of positive selection, including genetic co-selection, compensatory evolution that reduces plasmid cost, and horizontal gene transfer. In particular, conjugation has been proposed as a mechanism by which plasmids can be maintained even when costly, provided that transfer occurs sufficiently rapidly. This idea leads naturally to the concept of an invasion threshold: a critical conjugation rate above which plasmids increase when rare and below which they are eliminated. 

Early mathematical models of conjugative plasmids formalised this idea by identifying conditions under which horizontal transfer compensates for fitness cost and segregational loss. These analyses showed that plasmid fate is governed by the relative magnitudes of three key processes: plasmid-associated fitness, plasmid loss during cell division, and conjugation efficiency. Fast transfer and low cost favour plasmid persistence, while slow transfer and high cost promote elimination. Importantly, theoretical results predict that even purely parasitic plasmids can persist if the transfer rate exceeds a critical value. 

Despite these predictions, the ecological relevance of conjugation-mediated persistence has been debated. Some studies have argued that the transfer rates required to maintain costly plasmids are unrealistically high, while others have demonstrated experimentally that conjugation can indeed be fast enough to support plasmid invasion. Differences in experimental systems, host strains, plasmids, and environmental conditions complicate direct comparisons further.

Mathematical models provide a natural framework for deriving such criteria. By analysing plasmid invasion at low frequency, it can be determined whether plasmids can establish in an otherwise plasmid-free population. This approach yields explicit invasion thresholds that separate regimes of plasmid extinction from persistence and dominance. These thresholds depend not only on plasmid cost and loss, but also on host demographic parameters and selective pressures acting on plasmid-free cells. As a result, selection can lower invasion barriers and qualitatively change long-term population outcomes, even without directly benefiting plasmid-bearing cells [5, 6].

### Study Aims and Objectives

The aim of this study is to investigate the conditions under which conjugative plasmids can invade and persist in bacterial populations, despite imposing fitness costs on their hosts. Using a combination of analytical and computational approaches, we seek to clarify how plasmid cost, segregational loss, horizontal transfer, and environmental selection jointly determine plasmid fate.

Specifically, this study aims to: 
- Derive analytical invasion thresholds that define the minimum plasmid transfer rate required for persistence. 
- Quantify how plasmid cost and loss contribute to these thresholds and how selective pressure modifies invasion conditions. 
- Compare deterministic predictions from ordinary differential equation models with outcomes from stochastic simulations. 
- Identify parameter regimes in which stochastic extinction or rare survival deviates from deterministic expectations. 

These objectives provide a unified workflow for interpreting plasmid invasion dynamics and for assessing when conjugation-mediated persistence is likely to occur.

---

## Methods
### Deterministic ODE Model

The deterministic model describes the temporal evolution of plasmid-free and plasmid-bearing populations using coupled nonlinear ODEs (Appendix A.1). Population growth follows logistic dynamics with a shared carrying capacity, while plasmid carriage imposes a growth cost on plasmid-bearing cells. Plasmid transfer is modelled using mass-action kinetics proportional to the product of plasmid-free and plasmid-bearing cells, and plasmid loss during cell division is represented as a constant segregational loss rate. 

All model parameters and their biological interpretations are summarised in Appendix A.2, and the core modelling assumptions are outlined in Appendix A.3. These assumptions include homogeneous mixing, time-invariant parameters, and full immunity of plasmid-bearing cells to selective pressures acting on plasmid-free cells. 

The ODE system was numerically integrated using standard time-stepping methods implemented in Python. Time-series simulations were performed to characterise transient dynamics and long-term outcomes across parameter regimes.

### Equilibria and Invasion Analysis

Analytical expressions for the plasmid-free and plasmid-bearing equilibria were derived by solving the steady-state conditions of the ODE system (Appendix B.1). Local stability of these equilibria was assessed using the Jacobian matrix of the system (Appendix B.2). 

Plasmid invasion conditions were determined by evaluating the invasion eigenvalue at the plasmid-free equilibrium (Appendix B.3). This analysis yields a threshold condition on the conjugation rate, defining the minimum plasmid transfer rate required for successful invasion. Expressions for critical values of the transfer rate and selective pressure were derived (Appendix B.4), including linearization of the critical transfer rate in respect to plasmid loss and cost (Appendix B.5). 

To visualise these invasion thresholds, ODE simulations were used to generate bifurcation plots showing plasmid persistence as a function of the conjugation rate $β$, as well as two-parameter heatmap sweeps illustrating the combined effects of plasmid transfer rate and plasmid cost on long-term plasmid prevalence.

### Stochastic SSA Model

To capture stochastic effects not represented in the deterministic system, a Gillespie-style SSA was implemented (Appendix C). The SSA explicitly simulates individual-level birth, death, plasmid loss, and conjugation events, with event propensities derived from the same biological processes as the ODE model. 

The full set of reactions and propensity functions is given in Appendix C.2. Logistic growth was incorporated by scaling birth rates with the total population size. The simulation proceeds by sampling reaction times from an exponential distribution and selecting reaction events probabilistically according to their relative propensities (Appendix C.3). 

Multiple stochastic realisations were generated to produce time series of plasmid-free and plasmid-bearing populations. These simulations were used to examine variability between runs, extinction events, and stochastic delays in plasmid invasion.

### Outputs and Visualizations

For both models, time-series outputs of $F$ and $P$ were recorded and visualised. ODE simulations were further used to construct $β$ - $c$ heatmaps and $β$-dependent bifurcation plots to summarise invasion and persistence regimes. All figures were generated using the plotting utilities contained within the respective `ode_model/` and `ssa_model/` directories.

---

## Results
### Deterministic Invasion Threshold and Equilibria

Analysis of the deterministic ODE model yields a clear invasion criterion for plasmid-bearing cells. Linear stability analysis at the plasmid-free equilibrium shows that plasmids can invade only if the conjugation rate $β$ exceeds a critical threshold:

$$
β_c = \frac {δ+μ(c/1-s)} {(1 - μ/r(1-s))}
$$

where $δ$ is the segregational loss rate, $c$ is the plasmid cost, $μ$ is the mortality rate, $r$ is the baseline growth rate, and $s$ represents selective pressure acting against plasmid-free cells.

In the absence of selective pressure ($s=0$), this simplifies to:

$$
β_c = \frac {δ+μc} {(1 - μ/r)}
$$

Under the standard parameter values, this critical transfer rate evaluates to $β_c=0.0167$.

For $β<β_c$, plasmids cannot invade and the plasmid-free equilibrium is locally stable; for $β>β_c$, plasmids grow when rare and persist in the population. 

The threshold expression reveals that the minimum transfer rate increases additively with both plasmid loss and plasmid cost, amplified by the demographic factor $1/(1−μ/r)$. Even in the absence of segregational loss ($δ=0$), plasmid invasion requires a transfer rate proportional to plasmid cost, while loss imposes an additional independent constraint. 

Existence conditions further constrain long-term outcomes. The plasmid-free equilibrium exists only when $μ<r(1−s)$, while the plasmid-bearing-only equilibrium requires $μ+δ<r(1−c)$. As mortality or plasmid cost increases, these equilibria disappear via feasibility loss, independently of invasion stability.

### Deterministic Dynamics Across Selective Pressure

ODE simulations are consistent with the analytical invasion condition. When $s=0$ and $β$ is below the critical value, plasmids fail to invade and decay deterministically to extinction. Increasing selective pressure lowers the effective invasion threshold by penalising plasmid-free cells, allowing plasmids to persist at progressively lower transfer rates. 

At intermediate selective pressure ($s=0$.2), plasmids invade successfully when $β>β_c$, eventually becoming the dominant population while maintaining coexistence with plasmid-free cells. At higher selective pressure ($s=0.6$), plasmids rapidly dominate the system, driving plasmid-free cells to very low equilibrium densities. 

As $s→r$, the existence condition for the plasmid-free equilibrium $μ<r(1−s)$ fails, making extinction of plasmid-free cells increasingly likely even in the absence of transfer dynamics.

### Stochastic Effects and SSA Outcomes

Stochastic simulations broadly reproduce the qualitative behaviour of the deterministic model but reveal important deviations near extinction boundaries. At $s=0$, most SSA realisations lead to plasmid extinction (98%), consistent with the ODE prediction. However, a small fraction of runs show transient plasmid persistence at very low population sizes ($0.0007±0.0049$), reflecting demographic noise and rare-event survival that cannot occur in the deterministic framework. 

At $s=0.2$, the majority of stochastic realisations closely match the ODE dynamics: plasmids successfully invade, increase in frequency, and become the dominant population ($0.5513±0.0368$ versus $0.3376±0.0340$). No plasmid extinctions were observed under these conditions. 

At high selective pressure ($s=0.6$), plasmids effectively take over the population in all stochastic runs ($0.7624±0.0202$), with plasmid-free cells persisting only at very low abundances ($0.1251±0.0141$). No plasmid extinctions occur in this regime, and stochastic variability primarily affects the timing rather than the outcome of invasion. 

Overall, stochasticity is most influential near the deterministic invasion threshold and when plasmid populations are rare, while strong selective pressure stabilises plasmid persistence and suppresses extinction risk.

## Conclusions

This study demonstrates that plasmid invasion and persistence are determined by a balance between plasmid cost, segregational loss, conjugation rate, and selective pressures on plasmid-free cells. The deterministic analysis yields a clear invasion threshold, showing that plasmids require a minimum transfer rate proportional to both their fitness cost and loss rate to establish in a population. ODE simulations confirm that selective pressure lowers this threshold, facilitating plasmid persistence even when costs are substantial. 

Stochastic simulations reveal that demographic noise introduces variability near the invasion boundary, allowing rare plasmid survival or delayed extinction in low-density populations. However, under moderate to strong selective pressures, stochastic outcomes closely match deterministic predictions, with plasmids reliably invading and dominating the population. 

Overall, these results highlight that conjugation-mediated plasmid persistence is not solely determined by fitness cost, but by the interplay of horizontal transfer, loss, and selection. The derived critical transfer rate provides a general criterion to predict plasmid fate, consistent with other theoretical studies on plasmid dynamics and antibiotic resistance spread.

## References

1. Hernández-Beltrán, J.C.R., San Millán, A., Fuentes-Hernández, A. and Peña-Miller, R. (2021). Mathematical Models of Plasmid Population Dynamics. Frontiers in Microbiology, 12. Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC8600371/ 
2. Stewart, F.M. and Levin, B.R. (1977). THE POPULATION BIOLOGY OF BACTERIAL PLASMIDS: A PRIORI CONDITIONS FOR THE EXISTENCE OF CONJUGATIONALLY TRANSMITTED FACTORS. Genetics, 87(2), pp.209–228. Available at: https://pmc.ncbi.nlm.nih.gov/articles/instance/1213735/ 
3. Lopez, J.G., Donia, M.S. and Wingreen, N.S. (2021). Modeling the ecology of parasitic plasmids. The ISME Journal, 15(10), pp.2843–2852. Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC8443676/ 
4. Prensky, H., Gomez‐Simmonds, A., Uhlemann, A. and Lopatkin, A.J. (2021). Conjugation dynamics depend on both the plasmid acquisition cost and the fitness cost. Molecular Systems Biology, 17(3). Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC7919528/
5. Lopatkin, A.J., Meredith, H.R., Srimani, J.K., Pfeiffer, C., Durrett, R. and You, L. (2017). Persistence and reversal of plasmid-mediated antibiotic resistance. Nature Communications, 8(1), pp.1–10. Available at: https://pmc.ncbi.nlm.nih.gov/articles/PMC5698434/ 
6. Ibarguen-Mondragón, E., Esteva, L., M. Victoria Otero-Espinar, Vega, E. and Miller Cerón-Gómez (2025). On qualitative properties of replication and transfer of conjugative plasmids encoding antibiotic resistance genes. Computational and Applied Mathematics, 44(5). Available at: https://link.springer.com/article/10.1007/s40314-025-03231-w