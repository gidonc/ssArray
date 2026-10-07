# ssArray: basic expectations for the six parameterisations (single area)

Working summary, updated 2026-10-07 after experiments 04 to 24 and the recheck of 16 and 17 with the current model. Evidence is from gidonc/ssArray, branch real-table-experiment (PR #7). Short runs (500+500 draws, two chains). Treat numbers as indicative.

**Current model** (stan/single_area.stan; density on tables unchanged by any of these):
- margin coordinates: column log-ratios against the largest column (K_margin; 'whitened' and 'last' also available);
- centred interior for adjusted row and adjusted table (centred_interior = 1): contrasts with no reference cell;
- split basis evaluated without its matrix (split_basis = 1); adjusted table's Jacobian by the determinant lemma.
Sections marked (old model) have not been rerun and used last-column margins and reference-cell interiors.

## Principle

HMC cost per independent draw grows with (largest posterior scale / smallest posterior scale) left after the metric (toy Gaussian, experiment 12).

- A tight direction that is a single parameter costs nothing: a diagonal metric rescales it. Stan's metric regularisation puts a floor of about 0.005 on the adapted sd, so tighter parameters must be scaled in the model (experiment 11).
- A fixed linear mixture of parameters needs a dense metric or a linear change of coordinates in the model. **Any reference category is a case: its noise enters every parameter measured against it.** This was found three times: the reference column of the margins (20-22), the reference column of each row in adjusted row, and the reference row and column in adjusted table (24).
- Scales that change across the posterior (nonlinearity, kinks) cannot be fixed by any metric. What matters is the nonlinearity of the map from parameters to each quantity the density is built from, over a posterior-sized step, weighted by how tightly the density constrains that quantity.

## Components by parameterisation

Nonlinearity r = |actual change - linear prediction| / |linear prediction| over one posterior sd; 0 = linear.

| component | ILR + log volume | log cells | position logit | odds-ratio logit | adjusted row | adjusted table |
|---|---|---|---|---|---|---|
| margins (tight, from data) | nonlinear, r 0.09-0.12 | nonlinear, r 0.12-0.17 | exact coordinate | exact coordinate | exact coordinate | exact coordinate |
| row and column effects (prior) | linear | linear | r 0.09-0.16 | r 0.08-0.14 | r 0.07-0.11 | r 0.09-0.12 |
| interactions (prior) | linear | linear | r 0.05-0.28 | r 0.04-0.13 | r 0.02-0.04 | 0 exactly |
| kinks | none | none | yes | none | none | none |
| tight directions isolated in parameters | no | no | yes | yes | yes | yes |
| time per gradient at 10x10, current model (ms) | 0.03 | 0.02 | 0.06 | 0.08 | 0.11 | 0.18 |

Nonlinearity and isolation rows: experiments 07, 10, 13 on the Scottish grand total (old model; the nonlinearity of a map does not depend on a linear change of its parameters, so these should stand).

## Results with the current model

Experiment 16: Scotland 2007 and New Zealand 2017 constituencies (5x5, 20,000-40,000 votes), senc (3x3, about 1,600 per precinct, some margins of 3-8), redistrict (3x3, 500-650, margins only). Five areas x two seeds, eps 1, sigma_b 0.5 and 2, four prior set-ups; 1,296 sequential fits, direct schemes from the earlier run (unchanged by the model changes). Experiment 17: simulated 5x5 baseline (total 30,000, sigma_b 0.5), one feature changed per arm, three jittered replicates.

1. **The raking schemes are now the most efficient per gradient at both prior widths** (this replaces "prior width separates the schemes"). Elections, median ESS per 1,000 gradients with the raked national centre: sigma_b 0.5: position 61-84, odds-ratio 57-103, adjusted row 121-137, adjusted table 122-136; sigma_b 2: position 2-3, odds-ratio 5-7, adjusted row 25-36, adjusted table 38-41. Best scheme per cell: adjusted row or table in 35 of 40 (Scotland) and 40 of 40 (NZ) at sigma_b 0.5, and 40 of 40 in both at sigma_b 2. Against the old model adjusted row gained a median 2.5-3.5 times and adjusted table 4.5-6.4 times in the elections (1.3-2.4 times in the 3x3 data); the cell-by-cell schemes gained 1.0-1.4 times.
2. **On time, prior width still decides.** sigma_b 0.5: position or odds-ratio logit fastest in 34 of 40 (Scotland) and 22 of 40 (NZ) cells, because their gradient is cheaper and all schemes take 7-12 leapfrogs. sigma_b 2: adjusted row or table fastest in 39 of 40 and 40 of 40; adjusted row usually (62 cells against 17). Median ESS per second at sigma_b 2 (raked centre, Scotland and NZ): position 141 and 94, odds-ratio 299 and 206, adjusted row 488 and 772, adjusted table 406 and 496.
3. **In the 3x3 data the four are close.** Per gradient the raking schemes win 18 of 38 (senc) and 17 of 29 (redistrict) cells at sigma_b 0.5 and 24 of 39 and 26 of 29 at sigma_b 2; on time the cell-by-cell schemes win most cells at sigma_b 0.5.
4. **Wide prior:** leapfrogs per draw at sigma_b 2 in the elections: position 110-140, odds-ratio 31-39, adjusted row 14, adjusted table 12. Simulated arm (relative to baseline): position x0.06, odds-ratio x0.05, adjusted row x0.29, adjusted table x0.41.
5. **Diagonal centre with a wide prior** (simulated): position x0.02, odds-ratio x0.04, adjusted row x0.13, adjusted table x0.30; adjusted table ahead of adjusted row in all three replicates (41 against 23 ESS per 1,000 gradients).
6. **Prior on interactions only:** adjusted table's parameters are then the prior's own coordinates. Elections: adjusted table gains 1.4-2.3 times and is the best scheme per gradient (Scotland 188 at sigma_b 0.5 and 82 at sigma_b 2, against 91 and 32 for adjusted row) and level or ahead on time; adjusted row about level; position and odds-ratio lose (x0.2-0.5 at sigma_b 2). **Still fails where a margin is tiny:** senc 2-8 of 10 runs per scheme not converged, redistrict 1-3 of 10; simulated tiny-margin arm adjusted table 3 of 3 failed. With equal sds the same margins are harmless.
7. **Margin conflict from a common centre:** sequential schemes 1 non-converged run of 312 with conflict below 14 prior sds, 5 of 8 at 18.9 sds. Cost where it converges: about 20% in Scotland, none in NZ and senc. Simulated: 5 sds x0.5-1.1, 20 sds x0.1-0.3.
8. **Margin tightness decides direct against sequential, and the direct schemes' range has shrunk.** Per gradient the best sequential scheme is ahead of the best direct scheme down to about 300 votes per table in the elections (59-73 against 18-29) and level at 16 in senc; direct is ahead only in redistrict at 5 units per table. On time direct is ahead only at the lowest volume in each dataset (5-350 per table). At about 30,000 per table direct is 40-90 times worse per gradient, and in the millions it is several hundred times worse or fails. Simulated: total 30 gives direct x20-29 and all sequential schemes x0.3-0.55.
9. **Table size.** Simulated 3x3 and 8x8 relative to 5x5: position x2.0 and x0.73, odds-ratio x1.2 and x0.48, adjusted row x1.0 and x0.46, adjusted table x1.3 and x0.66. Large tables (24): every scheme converges up to 48x48 (2,304 cells) in 30-110 s for 2 chains x 1,000 iterations, with 15-63 leapfrogs per draw (about 6,000 per margin, sigma_b 0.5, independence centre, one replicate).

## Model changes and what they fixed

10. **Margin reference column (20-22).** Column margins were log-ratios against the last column; with it at 0.4% of the total the column parameters correlate at 0.99 and ESS per gradient falls about ninefold. Largest-column reference or whitening removes it; the two perform within 2-5% of each other. A residual of about two times remains for which column sits last in the fill order.
11. **Reference cells in the raking schemes (24).** Correlation of about 0.5 between a row's parameters, so conditioning worsens with the number of columns: adjusted row's min ESS fell from 202 at 10x10 to 14 (not converged) at 32x32. Centred interior: 650-1,450 at every size.
12. **Cost (24).** Dense basis product: no measurable cost up to 10x10; 2 times at 16x16, 3.5 times at 24x24, 13 times at 32x32. Adjusted table's Jacobian: (R-1)(C-1) square determinant reduced to (R+C-1) square; time per gradient at 10x10 from 1.15 ms to 0.18 ms; adjusted table is now about 2.5 times odds-ratio logit per gradient at every size.
13. **Margin scaling (11).** Margin parameters need scaling by their penalty width whenever their posterior sd is below about 0.005-0.01, because Stan's adapted metric cannot go lower.

## What carries over to the EI model (ssEI not changed)

- Reference categories: the ssEI margin parameterisation uses the last column as reference (mp_s_obs); any raking-style interior there would have the reference row/column issue. Cost grows with the number of categories: mild at 3x3, clear by 5x5.
- The dense basis product is per area and does not matter at EI table sizes. What grows in EI is the number of areas: any matrix spanning areas is the analogue (in the ssEI snapshot ROT_red, about 1,483 square for senc, is multiplied in at every gradient when lflag_rot_llrep = 1).
- The determinant-lemma Jacobian matters for adjusted table from about 5x5.
- Margin scaling and the metric floor (item 13).

## Simulator and other findings (old model unless stated)

- experiments/simdesign.py: separate dials for size (shape; total, per-margin or per-cell count held fixed; eps), margins (balance and family per side, alignment, jitter), interior (diagonal log odds ratio) and prior (sigma_b, sd on row/column effects, conflict); each design reports derived features.
- For ILR the count per margin is the volume that matters and the largest margin sets the cost; a tiny margin does not (19). With counts per margin fixed, ILR's margin ESS is flat across sizes and its cell ESS rises and saturates by about 8x8 (18).
- A 'tied margins' dial did not control bound flips; no design yet isolates kinks in position logit.
- Diagonal centre with a narrow prior: log cells 2.7-3.8 times ILR (17).

## Not yet tested

- Priors with unequal sds across interaction coordinates, and the choice of basis within the interaction part.
- A moderate (not near-flat) prior on row and column effects, as the fix for tiny margins under an interactions-led prior.
- The fill-order residual in item 10; large tables with a wide prior, unbalanced margins or large counts.
- The single-table diagnostics (07, 10, 13) and the first regime sweep (08) with the current model.
- An estimated centre: cross-area models (E_rc, sigma_jrc), centring, and the coupling between hyperparameters and areas.
- Areas with zero margins or structural zeros (excluded so far: 17 of 71 NZ electorates, 29 of 212 senc precincts, 83 of 150 redistrict precincts).
