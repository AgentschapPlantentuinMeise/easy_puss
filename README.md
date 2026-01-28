### easy_puss  
*A lightweight toolkit for processing, cleaning, and interpreting ICP-OES analytical data.*

easy_puss provides a set of convenience classes to streamline the workflow around ICP-OES data. It focuses on reading raw exports, cleaning and harmonizing element measurements, applying LOD/LOQ rules, and preparing results for downstream analysis or reporting.

The package is designed to be simple, predictable, and easy to integrate into existing Python data pipelines.

---

## 🌱 Features

- **Load PUSS result files** into tidy pandas DataFrames  
- **Automatic type conversion** (strings → floats, numeric cleanup, NA handling)  
- **LOD/LOQ masking**  
  - < LOD formatting  
  - > LOD formatting  
  - Optional flagging mode (*, **)  
- **Element-wise operations** with proper alignment of LOD/LOQ tables  
- **Helpers for MultiIndex structures**  
- **Utility functions** for common analytical workflows

---

## 📦 Installation

```bash
pip install easy_puss
```

Or install directly from GitHub:

```bash
pip install git+https://github.com/AgentschapPlantentuinMeise/easy_puss
```

---

## 🚀 Quick Start

```python
import easy_puss as ep
import pandas as pd

# Load results and LOD/LOQ tables
project = ep.ICPOES("results.csv")

# Apply masking
masked = project.apply_lod_loq_mask(project.results, project.lodq)

# Or apply flagging instead
flagged = project.apply_lod_loq_flags(project.results, project.lodq, flag = True)
```

---

## 🔬 LOD/LOQ Handling

easy_puss implements the standard interpretation rules:

| Condition | Output (mask mode) | Output (flag mode) |  
|----------|---------------------|---------------------|  
| value < LOD | < LOD | value* |  
| LOD ≤ value < LOQ | > LOD | value** |  
| value ≥ LOQ | numeric | numeric |

All operations are vectorized and align element names automatically.

---

## 📁 Typical Workflow

1. Export raw ICPOES data  
2. Load into Python using easy_puss  
3. Clean and harmonize element columns  
4. Apply LOD/LOQ rules  
5. Export masked or flagged results

## Disclaimer

Developed together with Copilot AI
