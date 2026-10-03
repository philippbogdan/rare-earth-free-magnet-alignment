# Goal

Establish, for the rare-earth-free hard magnetic phases under serious study, which processing routes can align their easy axes, the most remanence each route could deliver, and why, so that candidates can be ranked by the magnet they could become rather than by intrinsic properties alone. Aspire to the best that could exist, not a good one. What is known today is a floor, not a target.

## Specific

1. **The account:** for each phase, the routes that can align it and those that cannot, each with its mechanism and the remanence ceiling it implies.
2. **The tests:** for each prediction, the single measurement that would confirm or refute it.
3. **`FINDINGS.md`** (what the work establishes, each claim with its evidence), **`REPORT.md`** (each criterion against the floor, what failed, what remains open), and one command that regenerates every number.

## Measurable

Every number is checked against references this work does not produce.

- **Fidelity:** the account reproduces the known cases (field-aged AlNiCo, the lost texture in printed MnAl(C), textured MnAl-C by hot extrusion, field-pressed NdFeB), with their published remanence where measured.
- **Accuracy:** each predicted ceiling against published measured remanence or texture wherever one exists, including cases held back while the account was built.
- **The frontier:** the non-dominated set of phases covered, routes covered and predictive accuracy. After every success, push it on whichever criterion is furthest from its limit.
- **Limits:** only physics, or an unpublished quantity shown to decide the answer with the measurement that would settle it. A scope, cut-off or tolerance chosen here and anything left untried are not limits.
- **Evidence:** every number comes from one named commit and a named source; no claim is stronger than its data (an abstract is an abstract, one sample is one sample); one consistent method across phases, not a choice made per case.

## Relevant

Whoever designs and screens rare-earth-free magnets needs to know which candidates can become aligned magnets, and by what route, before spending on processing trials.

## Done

Done means `COMPLETE.md` shows, through the measurements above, that all three deliverables hold, the known cases are reproduced, and every remaining limit belongs to the problem itself. There is no deadline and no other way to finish.

## Rules

- Every choice is yours.
- Never invent results or references; label estimates as estimates.
- Contact nobody and publish nothing: no email, posts, issue comments, pushes or pull requests.
- Work only inside this directory; do not read or search local files elsewhere. The internet is yours.
- Commit as you go, so the repository alone can resume the work. Do not edit PROBLEM.md or GOAL.md.
- The machine: Apple M4, 10 cores, 16 GB memory, no FP64 GPU, shared with other long runs. Past 75% memory use the largest process any run started here is stopped, and disk is short, so keep data and environments here small.
- Quick scripts run here directly. Heavier work goes through `gpurun`, run from inside this directory: to Imperial's clusters (GPUs, CPU servers and a pool of lab machines), or to this Mac with `--on mac`, where jobs wait their turn so the Mac stays responsive. `gpurun --help` says what each offers and how to ask for it. Results come back by themselves when a job ends; `gpurun pull` fetches a file sooner.
