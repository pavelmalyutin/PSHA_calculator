Вот дополненная документация:

```markdown
# PSHA_calculator
This program, based on the Gusev-Shumilina model, calculates the intensity according to a synthetic catalog. The approach used in this program is fully consistent with the 2016 General Seismic Zoning approach (OSR-2016)


## Compilation

Basic:
g++ -O3 -fopenmp -std=c++17 -o psha_intensity psha_intensity.cpp -lm
Optimized for your CPU:
g++ -O3 -march=native -mtune=native -fopenmp -std=c++17 -ffast-math -funroll-loops -o psha_intensity psha_intensity.cpp -lm


Execution
./psha_intensity [your_synthetic_catalog] [parameters_table] [years_in_one_cycle] [number_of_cycles] [grid_file] [output_file]

Example:
./psha_intensity CTL.txt param_table.txt 10000 100 grid.geg output.txt
This means: catalog covers 10000 years per cycle × 100 cycles = 1,000,000 years of seismicity.

---

## File [your_synthetic_catalog]
Tab-separated file with 9 columns and header:

INDZ	MAG	L	W	AZ	DIP	PHI1	LMD1	H1
12	6.23	15.2	8.1	45.0	75.0	55.123	160.456	12.5
12	5.80	10.5	5.8	120.0	60.0	54.987	159.234	15.0

| Column | Description |
|--------|-------------|
| INDZ | Zone number according to "Explanatory note..." ([reference](https://scholar.google.ru/citations?view_op=view_citation&hl=ru&user=XpqWIJIAAAAJ&cstart=20&pagesize=80&sortby=pubdate&citation_for_view=XpqWIJIAAAAJ:vDijr-p_gm4C)). **Use 12 for non-Kamchatka regions** |
| MAG | Moment magnitude Mw |
| L | Rupture length (km) |
| W | Rupture width (km) |
| AZ | Azimuth of the upper edge (degrees, 0-360). Direction of strike measured clockwise from North |
| DIP | Dip angle (degrees, 0-180). Angle between horizontal plane and rupture plane, measured clockwise when viewed along strike direction |
| PHI1 | Latitude of the rupture corner X1 (degrees) |
| LMD1 | Longitude of the rupture corner X1 (degrees) |
| H1 | Depth of the rupture corner X1 (km) |
**Note:** PHI1, LMD1, H1 define the corner point X1 of the rupture plane (not the centroid).

## File [parameters_table]

Tab-separated file with zone parameters and header:
ind		sdevm	sdevi	parfln
12		0.5	    0.7	    45imr	
15		0.4	    0.6	    11imr	

| Column | Description |
|--------|-------------|
| ind | Zone number (must match INDZ in catalog) |
| sdevm | Standard deviation for magnitude-dependent intensity scatter |
| sdevi | Standard deviation for individual event intensity scatter |
| parfln | Attenuation model name (see below) |


## File [grid_file]

Space or tab-separated file with calculation points (latitude, longitude):

55.0 160.0
55.1 160.0
55.2 160.0
55.0 160.1


## Output file

Tab-separated file with 3 columns:
55.000000	160.000000	7.25
55.100000	160.000000	7.18
55.200000	160.000000	7.05

| Column | Description |
|--------|-------------|
| 1 | Latitude (degrees) |
| 2 | Longitude (degrees) |
| 3 | Intensity with 500-year return period (T500) |

## Algorithm overview

1. **Load synthetic catalog** — events with rupture geometry
2. **Project coordinates** — geographic to local Cartesian (km) using GEDECCON projection
3. **For each grid point:**
   - Loop over all events
   - Calculate 3D distance to rupture plane
   - Apply finite-fault correction (FINCOR)
   - Compute intensity using Gusev-Shumilina attenuation model
   - Add random scatter (σ_M and σ_I)
   - Accumulate intensity histogram
4. **Extract T500** — intensity exceeded with probability 1/500 per year

### Intensity formula
I = I_bas + α_M × (M - M_bas) + α_A × log₁₀(FCOR × ATT(R) / DENOM)
Where:
- `I_bas` — reference intensity
- `α_M` — magnitude scaling coefficient
- `α_A` — attenuation scaling coefficient  
- `FCOR` — finite-fault correction factor
- `ATT(R)` — distance attenuation function
- `DENOM` — normalization factor

Random scatter is added:
I_result = I + ε_i × σ_I + ε_M × σ_M
