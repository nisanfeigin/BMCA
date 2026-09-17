# Status against Nature's Code and Software Submission Checklist

| Checklist item | Status |
|---|---|
| Source code | **Done** — `R/`, five scripts with headers |
| Small dataset to demo the code | **Done** — `demo/`, simulated, ~0.2 MB |
| 1. System requirements: dependencies and OS, with version numbers | **Needs one command** — run `Rscript R/capture_session_info.R` and commit the output; the README points at it for exact versions |
| 1. Versions the code has been tested on | **Needs you** — same step; commit `sessionInfo.txt` and `package_versions.csv` |
| 1. Non-standard hardware | **Done** (none required) |
| 2. Installation guide: instructions | **Done** — README §2 |
| 2. Typical install time | **Done** — ~10–15 min binaries, ~30–45 min compiled |
| 3. Demo: instructions to run on data | **Done** — README §3, two commands |
| 3. Expected output | **Done** — described in README, reference copies in `demo/expected_output/` |
| 3. Expected run time for demo | **Done** — measured, ~4 seconds total |
| 4. How to run the software on your data | **Done** — README §4, with the object schemas |
| 4. Reproduction instructions (optional) | **Done** — README §4, run order given |

The only outstanding step is running `Rscript R/capture_session_info.R` on the
machine you use for the final figure run and committing `sessionInfo.txt` and
`package_versions.csv`. Those two files carry your R version, package versions
and OS release, which is what the first checklist section asks for, and the
README points to them rather than repeating the numbers.

Two things the checklist does ask for that are no longer stated anywhere: the
RAM needed for the full analyses, and the approximate total run time. Neither is
worth guessing at, but a reviewer may ask, so add a line each once you know
them. A licence and the data DOI/accessions also need adding before the
repository goes public.

The demo was written and executed here: both scripts run clean on R 4.3.3, are
seeded, and produce identical output on repeated runs. Run them once on your own
machine before submitting, since the checklist explicitly asks that a colleague
unfamiliar with the code install and run it.

---

# Code issues found but deliberately not auto-fixed

The cleanup was mechanical: no analysis logic, object name, parameter or
plotting call was altered, so the scripts still correspond exactly to the code
that produced the submitted figures. The items below need your judgement.

### 1. Filename case mismatch — will break on Linux

Three scripts read `BMCA/csv/all_MPS.xlsx`; `cell_state_tables.R` reads
`BMCA/csv/all_MPs.xlsx`. macOS is case-insensitive so this worked on your laptop,
but a reviewer on Linux will hit a file-not-found. Pick one spelling and make all
four scripts agree.

### 2. `figure_3.R` reads from outside the project tree

`readRDS("code/inferGEP/data/hotmap.rds")` appears 6 times and points at a
separate `inferGEP` directory. Nobody cloning the repo will have it. Since it is
only a colour palette, copy the object to `BMCA/data/hotmap.rds` and update the
six lines, or inline the palette as hex codes in `R/config.R`.

### 3. Missing random seeds

There are 18 calls to `sample()` (permutation tests, co-occurrence null models)
and only 5 `set.seed()` calls, in `figure_3.R` and `figure_5.R`. The permutation
analyses in the other scripts are therefore not exactly reproducible — which is
precisely what the checklist's reproduction-instructions item is about.

I did **not** add seeds, because doing so would change which random draws occur
and could shift p-values away from the submitted values. Add `set.seed()`
yourself, re-run those blocks, and confirm the numbers are unchanged within
tolerance before committing.

### 4. Patient-level clinical tables

`figure_5.R` reads `Clinical information for paired samples.csv` and
`Sample table with date-BrM-0730.csv`, which carry survival and date fields. The
`.gitignore` excludes `BMCA/csv/` entirely, so nothing is committed by default —
but check before the first push that no version of these files is staged, and
confirm what your data-sharing agreements allow before depositing them.

### 5. Typos baked into object filenames

`proccesed` (for *processed*) appears in several `.RDS` names, e.g.
`metadata_states_proccesed_filtered.RDS`, and `Mesenchimal` in
`score_tib_Mesenchimal.RDS`. Harmless, but if you rename the data files for
deposition, rename the references in the scripts at the same time.

---

# What the cleanup did

| | Removed |
|---|---|
| Duplicate/scattered `library()` calls | 373 → one declared block per script |
| `setwd()` and `getwd()` boilerplate with your home paths | 19 lines |
| Bare console-exploration calls (`unique(x)`, `table(x)`, `str(x)`, …) | 196 lines |
| Trailing whitespace, runs of 3+ blank lines, stray `#####` bars | throughout |

Every remaining line was verified present verbatim in your originals — nothing
was rewritten, reordered or reformatted, and bracket balance is unchanged in all
five scripts.

Each script opens with a header giving its purpose, the figure panels it
produces, every input file it reads (including paths built at run time, shown as
`BMCA/data/UMI_pb_<cell_type>.RDS`), and its original filename.
