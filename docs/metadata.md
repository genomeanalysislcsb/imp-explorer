# Sample metadata

Clinical, diet and other sample data enter a release with `imp-explorer export --metadata` and the cohort analysis
with `imp-explorer cohort --metadata`. `imp-explorer metadata` prepares and checks such a table first, so that a new
table can go in without surprises: a template with a **data dictionary**, and a **check** of ids, types, units,
categories and missing values.

## Template

```bash
imp-explorer metadata template --out metadata/ --data /path/to/release
```

writes

- `metadata.dictionary.tsv`: one row per column with `type`, `unit`, `allowed` (categories `a|b|c`, or a range
  `min..max`), `missing` (extra missing codes, e.g. `-9|999`), `category` (sample, clinical, medication, diet,
  lifestyle, technical) and `description`;
- `metadata.template.tsv`: the columns of the dictionary, one row per sample of the release (ids, groups and the
  metadata already in it filled in).

The dictionary starts with columns typical of a human gut microbiome cohort (sex, age, BMI, disease duration,
constipation, Bristol stool scale, antibiotics, PPI, levodopa / LEDD, methotrexate, smoking, alcohol, diet pattern,
energy, fibre and macronutrient intake, coffee, sampling date). Edit it for the study: delete columns, add others,
set units and the categories as they are spelled in the table.

| type | values | in the cohort analysis |
|---|---|---|
| `id` | sample id | – |
| `categorical` | categories, e.g. `F|M`, `never|former|current`; also codes like `0|1|2` | categorical |
| `binary` | two categories, e.g. `yes|no` | categorical |
| `integer` | whole numbers, optional range, e.g. `1..7` | numeric |
| `numeric` | numbers (a decimal comma is read as a dot), optional range | numeric |
| `date` | `YYYY-MM-DD` | not a covariate |
| `text` | free text | not a covariate |

## Check

```bash
imp-explorer metadata check --metadata clinical.tsv --dictionary metadata/metadata.dictionary.tsv \
    --data /path/to/release --report issues.tsv --write-clean clinical.clean.tsv
```

It lists per column the type, number of values and missing values and a summary (median and range, or the
categories), then the issues:

| level | examples |
|---|---|
| ERROR | duplicated or empty sample ids; values that are not numbers, not whole numbers, not dates; categories not in the dictionary |
| WARN | decimal commas; values outside the range (wrong unit?); categories that differ only in case; metadata ids not in the release and release samples without metadata; more than 30 % missing; groups (per block) with fewer than 5 samples |
| INFO | columns not in the dictionary; columns without values; constant columns; numeric codes without a dictionary entry |

The exit code is 1 when there are errors. Ids that differ from the release (`P-001` vs `P001`, other case) are
translated with the same `--id-map`, `--id-regex` and `--id-case` as for `export`. Missing values are empty cells or
one of `NA`, `n/a`, `nan`, `none`, `null`, `-`, `?`, `nd`, `n.d.`, `unknown` (any case) plus the dictionary's codes.

`--write-clean` writes the table as it would be used: ids translated, missing codes empty, numbers with a dot,
categories spelled as in the dictionary. Values with errors become empty there; fix them in the source table.

## Use in the cohort analysis

```bash
imp-explorer cohort --data /path/to/release --metadata clinical.clean.tsv \
    --dictionary metadata/metadata.dictionary.tsv --model "group+site, group+site+sex+age, group+site+fibre_g_day"
```

or `dictionary:` in the settings file. Without a dictionary, a column of numbers with more than two values is
numeric: `smoking` coded `0/1/2` would enter the models as a number. With the dictionary its type decides
(categorical), dates and free text are no covariates, and the dictionary's missing codes become NA.
