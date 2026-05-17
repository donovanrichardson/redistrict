# New York Census Tract Redistricting Experiments

**State:** New York (FIPS 36) | **Geography:** Census tracts (5,411) | **Districts:** 26 | **ncuts:** 10 | **niter:** 20

## Formula Comparison

| Run ID | Formula | Edge Weight | Worst Deviation | Notes |
|--------|---------|-------------|-----------------|-------|
| 26 | `original` | `SCALE / (dist/(2√pop_a) + dist/(2√pop_b))` | −4.4% (district 3) | Population-informed; one outlier |
| 27 | `uniform` | `SCALE` for all land edges | −8.6% (district 8) | **Preferred** — most compact, avoids urban fragmentation |
| 28 | `original_clamped` | Original + constant C so w_max/w_min = 4 | −2.5% (district 10) | Best population balance |
| 29 | `blend` | 50% norm(1/cost_orig) + 50% norm(1/dist) | −3.6% (district 25) | Middle ground |

## QGIS-Curated Result: run 43 (**preferred**)

| Run ID | Geography | Formula | water_penalty | Curated edges | Worst Deviation | Notes |
|--------|-----------|---------|--------------|---------------|-----------------|-------|
| 43 | block_groups | `uniform` | 40.0 (8×) | 62 water links removed (of 1,107) | −0.8% | Best result — uniform weights, QGIS-curated water links, zero-pop nodes excluded |

Parent run 41 (pending stub, block_groups, 15,739 active nodes, 331 zero-pop excluded).
Continued via: `redistrict --continue 41 --formula uniform --water-penalty 40.0`
GeoJSON: `run_43_curated_uniform_water8x.geojson`

## Selected Formula: `uniform` (run 27)

Run 27 produced the most geographically compact districts with the least unnecessary fragmentation of dense urban areas. The absence of a population signal in edge weights means METIS cuts purely on graph topology — preferring geographic boundaries rather than demographic ones.

## District Populations — All Runs

### Run 26: original

| District | Population | Deviation |
|----------|-----------|-----------|
| 0 | 782,191 | +0.7% |
| 1 | 779,458 | +0.3% |
| 2 | 773,560 | −0.4% |
| 3 | 742,731 | −4.4% |
| 4 | 772,014 | −0.6% |
| 5 | 782,152 | +0.7% |
| 6 | 782,579 | +0.7% |
| 7 | 779,817 | +0.4% |
| 8 | 773,947 | −0.4% |
| 9 | 774,942 | −0.3% |
| 10 | 779,135 | +0.3% |
| 11 | 781,165 | +0.5% |
| 12 | 780,940 | +0.5% |
| 13 | 771,938 | −0.6% |
| 14 | 782,093 | +0.7% |
| 15 | 783,154 | +0.8% |
| 16 | 778,807 | +0.2% |
| 17 | 780,613 | +0.5% |
| 18 | 782,689 | +0.7% |
| 19 | 781,434 | +0.6% |
| 20 | 774,386 | −0.3% |
| 21 | 781,074 | +0.5% |
| 22 | 779,072 | +0.3% |
| 23 | 771,458 | −0.7% |
| 24 | 778,197 | +0.2% |
| 25 | 771,703 | −0.7% |

### Run 27: uniform

| District | Population | Deviation |
|----------|-----------|-----------|
| 0 | 779,383 | +0.3% |
| 1 | 781,841 | +0.6% |
| 2 | 782,177 | +0.7% |
| 3 | 783,216 | +0.8% |
| 4 | 778,639 | +0.2% |
| 5 | 778,784 | +0.2% |
| 6 | 773,845 | −0.4% |
| 7 | 770,971 | −0.8% |
| 8 | 709,816 | −8.6% |
| 9 | 779,666 | +0.3% |
| 10 | 782,597 | +0.7% |
| 11 | 779,835 | +0.4% |
| 12 | 780,086 | +0.4% |
| 13 | 777,287 | 0.0% |
| 14 | 776,976 | 0.0% |
| 15 | 781,967 | +0.6% |
| 16 | 778,092 | +0.1% |
| 17 | 780,003 | +0.4% |
| 18 | 778,071 | +0.1% |
| 19 | 777,871 | +0.1% |
| 20 | 781,550 | +0.6% |
| 21 | 782,652 | +0.7% |
| 22 | 782,856 | +0.8% |
| 23 | 779,243 | +0.3% |
| 24 | 783,134 | +0.8% |
| 25 | 780,691 | +0.5% |

### Run 28: original_clamped

| District | Population | Deviation |
|----------|-----------|-----------|
| 0 | 775,174 | −0.2% |
| 1 | 775,056 | −0.2% |
| 2 | 778,074 | +0.1% |
| 3 | 776,162 | −0.1% |
| 4 | 782,713 | +0.7% |
| 5 | 777,137 | 0.0% |
| 6 | 772,373 | −0.6% |
| 7 | 773,824 | −0.4% |
| 8 | 771,634 | −0.7% |
| 9 | 773,249 | −0.5% |
| 10 | 757,855 | −2.5% |
| 11 | 779,541 | +0.3% |
| 12 | 782,401 | +0.7% |
| 13 | 776,417 | −0.1% |
| 14 | 778,138 | +0.1% |
| 15 | 773,418 | −0.5% |
| 16 | 777,245 | 0.0% |
| 17 | 779,612 | +0.3% |
| 18 | 781,501 | +0.6% |
| 19 | 782,809 | +0.8% |
| 20 | 782,855 | +0.8% |
| 21 | 781,890 | +0.6% |
| 22 | 782,484 | +0.7% |
| 23 | 782,715 | +0.7% |
| 24 | 773,171 | −0.5% |
| 25 | 773,801 | −0.4% |

### Run 29: blend

| District | Population | Deviation |
|----------|-----------|-----------|
| 0 | 774,719 | −0.3% |
| 1 | 782,887 | +0.8% |
| 2 | 779,299 | +0.3% |
| 3 | 778,312 | +0.2% |
| 4 | 781,898 | +0.6% |
| 5 | 781,746 | +0.6% |
| 6 | 776,073 | −0.1% |
| 7 | 772,530 | −0.6% |
| 8 | 777,535 | +0.1% |
| 9 | 782,749 | +0.7% |
| 10 | 782,523 | +0.7% |
| 11 | 780,766 | +0.5% |
| 12 | 781,753 | +0.6% |
| 13 | 774,636 | −0.3% |
| 14 | 782,158 | +0.7% |
| 15 | 782,861 | +0.8% |
| 16 | 771,566 | −0.7% |
| 17 | 782,668 | +0.7% |
| 18 | 776,549 | −0.1% |
| 19 | 776,785 | 0.0% |
| 20 | 772,915 | −0.5% |
| 21 | 780,579 | +0.5% |
| 22 | 773,530 | −0.4% |
| 23 | 771,774 | −0.7% |
| 24 | 773,276 | −0.5% |
| 25 | 749,162 | −3.6% |
