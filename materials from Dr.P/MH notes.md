| Iteration | $\alpha$ | $\theta$ | $\beta_0$ | $\beta_\text{infected}$ | $\beta_\text{female}$ | $\beta_\text{age14}$ | $\beta_\text{age28}$ |
| --------- | -------- | -------- | --------- | ----------------------- | --------------------- | -------------------- | -------------------- |
| 1         | 2        | 1.5      | 2         | 0                       | 0                     | 0                    | 0                    |
| 2         | Step 1   | Step 1   | Step 2    | Step 2                  | Step 2                | Step 2               | Step 2               |
| ...       |          |          |           |                         |                       |                      |                      |
| 20,000    |          |          |           |                         |                       |                      |                      |

#### Step 1 -- MH on $\alpha,\beta$

- proposed $log(\alpha^{*},\beta^{*})$ from old $log(\alpha^{t-1},\theta^{t-1})$
- compute $r = \text{post}(\alpha^{*},\beta^{*}, [\beta^{t-1}]) / \text{post}(\alpha^{t-1},\beta^{t-1}, [\beta^{t-1}])$ and the jacobian
- update $(\alpha^{t},\beta^{t})$

#### Step 2 -- MH on $[\beta]$

- proposed $[\beta^{*}]$ from old $[\beta^{t-1}]$ 
- compute $r = \text{post}(\alpha^{t},\beta^{t}, [\beta^{*}]) / \text{post}(\alpha^{t},\beta^{t}, [\beta^{t-1}])$
- update $[\beta^{t}]$ 

#### Step 3 -- thin draws 

#### Step 4 -- intervals

- construct a table for each subset (e.g. ACO female age 14)

| Iteration after thinning | Parameters                | Day 2                   | Day 3  | ...  | Day End |
| ------------------------ | ------------------------- | ----------------------- | ------ | ---- | ------- |
| 1                        | $(\alpha,\theta,[\beta])$ | cdf                     | cdf    |      | cdf     |
| 2                        | $(\alpha,\theta,[\beta])$ | cdf                     | cdf    |      | cdf     |
| ...                      | ...                       | ...                     | ...    | ...  | ...     |
| 20,000/50                | $(\alpha,\theta,[\beta])$ | cdf                     | cdf    |      | cdf     |
|                          |                           | Bounds(0.025,0.5,0.975) | Bounds | ...  | Bounds  |
