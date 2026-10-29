# Post-hoc pseudobulk corrections that lost points on T1:val

Base file: **R1** — LS-WOT within-type continuation of E9.5 donors, λ 1, κ 0.63 (binomial thinning), seed 7. Board headline **57.6485** (submission `69c77e4c`, 4 Oct 2026; per-metric skill DES 47.5 · DCS 59.4 · MMD 57.4 · CSS 68.6). R1's pseudobulk showed exactly the pattern in the README: top-200 expressed genes displaced by −0.199 on average, 82 "technical" genes (ribosomal, mitochondrial, heat-shock, haemoglobin, sex) among the 300 largest |Δ|. Every entry below tried to remove that pattern after the cells were emitted.

| variant | what it did | board headline | Δ vs R1 | per-metric skill (DES · DCS · MMD · CSS) | submission |
|---|---|---|---|---|---|
| TA | gene-wise rescaling of non-zero values to hit pseudobulk targets | 56.3239 | −1.32 | 46.1–47.0 · 57.8–58.9 · 53.5–55.5 · 64.6–67.1 (range over the four) | `719bad76`, 2 Oct |
| TB | same, alternative target | 56.1369 | −1.51 | " | `bb900917`, 2 Oct |
| TBc | TB with compartment constraint | 56.3623 | −1.29 | " | `074c3cfc`, 2 Oct |
| TB15 | TB at 1.5× strength | 55.1301 | −2.52 | " | `17672b86`, 2 Oct |
| LTB6000 | un-thin the 6,000 most abundant genes | 50.8274 | −6.82 | 45.9 · 54.1 · 52.3 · 50.5 | `9a328f98`, 2 Oct |
| C3 | λ 0.25 with gene means held at E9.5 | 54.6288 | −3.02 | 43.1 · 58.1 · 53.0 · 67.2 | `c8f5c795`, 25 Sep |
| **S1** | frozen-gene donor replacement + two-compartment renormalisation (top-200 + frozen genes) | 56.1952 | −1.45 | 44.66 · 59.45 · 54.80 · 68.63 | `de194a61`, 26 Oct |
| S2 (control arm) | maturation donor tilt β 0.5, frozen genes excluded, *no* pseudobulk fix | 57.3612 | −0.29 | 47.68 · 59.43 · 56.63 · 67.98 | `4bf1ac01`, 26 Oct |

S1 is the instructive one. Pre-registered with hypothesis and falsifier, it met every local objective — top-200 mean Δ −0.199 → −0.057, technical genes in the top-300 |Δ| 82 → 4, pseudobulk correlation with the measured step unchanged — and lost **2.81 DES and 2.55 MMD** on the board. The falsifier ("headline < R1") fired and was honoured.

Variants built but dropped before submission because local diagnostics already showed the renormalisation undoing them: a log-space level de-trend of R1 (4 passes, top-200 mean Δ only −0.16 to −0.18); count-space continuation with binomial-down / Poisson-up moves, median-centred per level bin (ρ(Δ, level) returned to −0.77 after thinning).

## What the board was telling us

1. The level-dependent downward displacement in R1 is mostly the *continued E9.5 → E10.5 depth drop*, inherited by any model of the step — plus κ 0.63 thinning of the output, which on its own produces ρ(Δ, level) ≈ −0.94. Treating it as a technical artifact to be scrubbed removed real signal together with the depth signature.
2. Per-cell renormalisation with cell-dependent factors moves every cell in the PC basis in which MMD is fitted. The MMD loss is mechanical and independent of the pseudobulk intent (S1: −2.55; the TA/TB family: −2 to −4).
3. Local pseudobulk diagnostics did not discriminate: all of them called S1 an improvement. The only reliable readout was the board, with one pre-registered falsifier per file.

Fix that worked locally: thin the earlier stage to the donor's depth *before* fitting the step (`depth_neutral_fit.md`).

Sources: our board records (`board_R1_detail`, trajectory logs, 25-Oct board GET replies) as compiled in our calibration set; scorer `veckit` 0.1.1 commit 46d41e6.
