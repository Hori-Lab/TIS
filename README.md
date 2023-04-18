
## TIS Models
|     | Name  |  M+/M-  |   M++   | Exv. | #sep-P | #sep-S | #sep-B |
|-----|-------|---------|---------|------|--------|--------|--------|
| RNA | DT13  | DH      | N/A     | WCA  |    3   |    3   |    2   |
| RNA | DT15  | Coulomb | Coulomb | mWCA |    3   |    3   |    2   |
| RNA | DHT18 |         | N/A     |      |        |        |        |
| RNA | NHT19 | DH      | PMF     | mWCA |    4   |    4   |    2   |
| DNA | CHT18 | DH      | N/A     | WCA  |    3   |    3   |    2   |
| DNA | TBA   | DH      | PMF     |      |        |        |        |

## Molecular topology
 + Linear S: linear nucleic acid whose 5'-end starts from sugar
 + Linear P: linear nucleic acid whose 5'-end starts from phosphate
 + Circular: circular nucleic acid

#### Number of particles
|             |Linear S   |Linear P   |Circular (always P)   |
|-------------|------|------|------|
| nucleotides | N    | N    | N    |
| P           | N-1  | N    | N    |
| S           | N    | N    | N    |
| B           | N    | N    | N    |
| nmp         |3N-1  | 3N   | 3N   |

#### Number of bonds
|             |Linear S   |Linear P   |Circular (always P)   |
|-------------|------|------|------|
| PS          | N-1  | N    | N    |
| SB          | N    | N    | N    |
| SP          | N-1  | N-1  | N    |
| Total       | 3N-2 | 3N-1 | 3N   |

#### Number of angles
|             |Linear S   |Linear P   |Circular (always P)   |
|-------------|------|------|------|
| PSB         | N-1  | N    | N    |
| PSP         | N-2  | N-1  | N    |
| BSP         | N-1  | N-1  | N    |
| SPS         | N-1  | N-1  | N    |
| Total       | 4N-5 | 4N-3 | 4N   |

## Tertiary stacks
`tst-dist`
 1. share with secondary stack (`exclusive = False`)
 2. exclusively tertiary stack alone (`exclusive = True`)
