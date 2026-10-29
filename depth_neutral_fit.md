# Depth-neutral fit of a stage-to-stage change

**Problem.** The released stages differ in per-cell depth (E9.5 median 11,080 UMI; E10.5 6,723). Any model that estimates the change between a *previous* stage and a *donor* stage and then applies it to the donor cells (shift, velocity, transport, continuation) estimates depth-drop + biology and applies both.

**Recipe.** Before fitting, binomially thin the deeper stage's integer counts to the shallower stage's median depth, renormalise, and fit the change on the thinned pair. The emitted cells are the donor cells at their native depth; nothing is thinned on output.

```
ratio = median_UMI(donor) / median_UMI(previous)
if ratio < 1:  previous_counts ~ Binomial(previous_counts, ratio)   # per entry, one RNG stream
fit  Δ_state = mean(donor | state) − mean(thinned previous | state)   # in whatever space the model uses
emit donor cells moved by λ·Δ_state
```

**Evidence (local replica, 5 seeds, veckit 0.1.1 commit 46d41e6; local anchors):** reference E8.5, donors E9.5 *thinned to E10.5 depth* (ratio 0.607, so that the replica has the test board's depth structure), target E10.5.

| candidate | headline | DES | MMD | CSS |
|---|---|---|---|---|
| thinned donors verbatim | 64.2 ± 0.8 | | | |
| naive transport λ 1 (step fitted on raw E8.5 → E9.5) | 63.3 ± 0.6 | | +2 worse | 72 → 64 |
| **depth-neutral transport λ 1** | **66.8 ± 0.7** | +4.8 | +3.6 | +1.6 (vs naive) |
| depth-neutral, λ 0.5 | 66.4 | | | |

The naive row reproduces the public-board signature of our P3 shot 1 (CSS loss, MMD loss, no DES gain). Two further replica facts worth knowing: on a 2-stage replica, matching donor depth to the target by *selection* of donors was worth +6.6 and by thinning +4.4; thinning when the target is not shallower cost −2.4 (κ 0.8) to −7.0 (κ 0.63). Depth matching the *output* is therefore a bet on the hidden stage's depth; depth-neutralising the *fit* is not.

**Caveats.** Thinning in count space is only exact if counts are recoverable (they are for this release: max rounding error 0.0012). The depth-neutral step still contains the real E9.5 → E10.5 change, of which about two thirds survives (314 of 462 pseudobulk calls). Whether E12.5 is shallower or deeper than E10.5 is unknown; the board contract's cell count suggests deeper.

Reference implementation: [`depth_neutral.py`](depth_neutral.py) (recover counts → thin → renormalise → per-type Δ; ~30 lines, numpy/scipy only).
