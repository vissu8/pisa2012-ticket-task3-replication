# Replication materials

**Same answer, different processes: behavioral types among correct responders in PISA problem-solving log data**

This folder contains the derived data and the analysis code that reproduce the results reported in the manuscript (Ticket Task 3, item CP038Q03, PISA 2012 computer-based problem-solving assessment).

The original data are the publicly available PISA 2012 data (https://www.oecd.org/pisa/data/). All files below are derived from them.

## How to run

All scripts use the Python 3 standard library only. Place the scripts and the CSV files in the same folder and run each script from that folder:

```
python3 person_level_analysis_python.py
python3 substitution_audit_python.py
python3 raw_vs_substituted_python.py
python3 cluster_sensitivity_python.py   # takes several minutes
```

## Code

| File | Reproduces | Input files |
|---|---|---|
| `person_level_analysis_python.py` | Tables 7–10, cluster-level ANOVA, Cohen's d, person-level regressions with country-clustered standard errors, comparison before and after the first cancellation message (Section 4.4, Table 9), the N-CHECK definition check (Section 3.6) | `SeqN_CP2.csv`, `data2_cluster_kmedoids_k8.csv`, `score1_seqofseq.csv` |
| `substitution_audit_python.py` | Rank of the assigned key in the similarity ordering and agreement of the rules with the ordering (Section 4.1; Online Resource 1, Table S2), comparison with edit distance (Table S5) | `cancel_score1_notrank1.csv` |
| `raw_vs_substituted_python.py` | Dependence on the standardization method: rules alone and embedding alone versus the final assignment (Section 4.1; Online Resource 1, Table S6) | `score1_raw_sequences.csv`, `score1_seqofseq.csv`, `SeqN_CP2.csv`, `cancel_score1_notrank1.csv`, `data2_cluster_kmedoids_k8.csv` |
| `cluster_sensitivity_python.py` | Sensitivity of the cluster solution to K = 6 to 10 (Section 4.3) | `data2_cluster_kmedoids_k8.csv`, `SeqN_CP2.csv` |

## Values that should be reproduced

- Cluster-level ANOVA (cluster membership reconstructed from the five variables): F(7, 10060) = 81.36, η² = .054.
- Person-level regression: N-CHECK b = 0.57 (SE = 0.04); with log total time and other action sequences b = 0.54 (SE = 0.04); with country fixed effects b = 0.50 (SE = 0.03); before / after the first message b = 0.70 / 0.42.
- Variance in CP2: N-CHECK alone .066, five clustering variables (linear) .073, eight clusters .054. Cohen's d: 0.66 (I-1 vs. III-1), 0.50 (I-1, II-1, II-2 vs. III and IV).
- 11,143 correct responders, 5,742 distinct SOS, mean SOS length 5.29; 37.3% of comparison actions before the first message.
- Substitution: 2,475 sequences, 56 keys; for the 2,427 non-key sequences, mean rank of the first admissible key 5.50 (chance 12.64), first-ranked key admissible 41.7% (chance 11.1%), admissible key among the first three 68.9% (chance 28.7%).
- Rules alone: clusters reproduced for 99.5% of test takers, N-CHECK for 97.7%. Embedding alone: same code for 91.8% of all Cancel-concluding sequence instances and 53.1% of the 8,398 substituted instances; clusters reproduced for 82.8%.
- Cluster sensitivity: CP2 variance explained for K = 6 to 10: .033, .043, .054, .055, .059; adjusted Rand index with K = 8: .93 (K = 7) and .94 (K = 9); the K = 8 run reproduces the reported clusters.

## Data

`SeqN_CP2.csv` (10,068 rows): test takers who answered Ticket Task 3 correctly and have complete CP2 responses.
Columns: `ID` (PISA student ID), `country`, `seq` (sequence-of-sequences, SOS, as space-separated codes after substitution), `seqn` (SOS length), `totaltime` (s), `LastPause`, `FirstPause`, `Cancel` (number of Cancel clicks), `score` (Ticket Task 3, all 1), `CP2` (sum score of the CP2 item group excluding Ticket Task 3).
The file contains no cluster label; cluster membership is reconstructed in the code from the pattern of the five clustering variables, which maps to exactly one cluster.

`data2_cluster_kmedoids_k8.csv` (11,143 rows): clustering input and result for all correct responders.
Columns: `SCIBO2`, `SCDBO`, `BOBO`, `BOBOBO` (counts of the corresponding action-sequence codes), `DA` (= N-CHECK in the manuscript), `clustering` (cluster 1–8).
Cluster numbers: 1 = I-1, 2 = II-1, 3 = III-1, 4 = III-2, 5 = II-3, 6 = IV-1, 7 = IV-2, 8 = II-2. No ID column; the row order is identical to `score1_seqofseq.csv`.

`score1_seqofseq.csv` (11,143 rows): SOS data of all correct responders (one column, `sentence`: the substituted codes in temporal order). Row-aligned with the clustering file.

`score1_raw_sequences.csv` (11,143 rows): the action sequences of all correct responders before substitution, in the same row order as `score1_seqofseq.csv`. One column, `sequences`: the action sequences separated by ` | `; a sequence consisting only of a Cancel click without any selection is written as `cancel`.

`cancel_score1_notrank1.csv` (136,125 rows): cosine similarities (Max embedding) between each of the 2,475 distinct Cancel-concluding sequences of correct responders and the key sequences.
The row of the finally assigned key is not included. Columns: `seq`, `key`, `score` (cosine similarity), `rank` (1 = most similar), `c` / `c2` (number of Buy-Ok clicks of sequence / key), `what` / `whatk` (final-ticket code; `nb` = no purchase), `check1` (same number of Buy-Ok clicks), `check2` (same final-ticket code), `allocation` (rank of the finally assigned key; one value per sequence), `exact` (always FALSE).
The 56 keys are the 20 most frequent sequences plus 36 Buy-Ok-repeat variants (9 with two and 27 with three Buy-Ok clicks), which were consolidated into six codes after assignment.
