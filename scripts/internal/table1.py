"""Aggregate patient characteristics for Table 1 using recorded values."""
import math


def number(value):
    if value is None or str(value).strip() in ('', 'NA', 'NaN'):
        return None
    return float(value)


def median_iqr(values):
    values = sorted(x for x in values if x is not None and math.isfinite(x))
    if not values:
        return 'Not available'
    def q(p):
        pos = (len(values) - 1) * p
        lo, hi = math.floor(pos), math.ceil(pos)
        return values[lo] + (values[hi] - values[lo]) * (pos - lo)
    def fmt(x):
        return f'{x:.2f}'.rstrip('0').rstrip('.')
    return f'{fmt(q(.5))} ({fmt(q(.25))}–{fmt(q(.75))}); n={len(values)}'


def build_table1(rows):
    groups = {g: [r for r in rows if r['group'] == g] for g in 'ABCD'}
    output = []
    def continuous(label, field, healthy='Not applicable'):
        output.append([label] + [median_iqr([number(r[field]) for r in groups[g]]) for g in 'ABC'] + [healthy])
    def categorical(label, field, target, include_healthy=False):
        vals = []
        for g in 'ABCD':
            if g == 'D' and not include_healthy:
                vals.append('Not applicable' if field != 'sex' else 'Not available')
                continue
            available = [r[field] for r in groups[g] if r[field] not in ('NA', '', None)]
            n = sum(str(v) == str(target) for v in available)
            vals.append(f'{n} ({100*n/len(available):.1f})' if available else 'Not available')
        output.append([label] + vals)
    output.append(['Participants, n'] + [str(len(groups[g])) for g in 'ABCD'])
    categorical('Experimental batch 1, n (%)', 'Batch', 'Sep3', True)
    categorical('Experimental batch 2, n (%)', 'Batch', 'Dec2', True)
    continuous('Age, years', 'age', 'Not available')
    categorical('Male, n (%)', 'sex', '男')
    categorical('Female, n (%)', 'sex', '女')
    continuous('APACHE II score, median (IQR)', 'severity')
    continuous('ICU stay days, median (IQR)', 'icu_days')
    for field, label in [('steroids', 'Corticosteroid treatment'), ('anticoagulants', 'Anticoagulant treatment'), ('transfusion', 'Transfusion'), ('ECMO', 'ECMO'), ('death', 'Deaths')]:
        categorical(label + ', n (%)', field, '1')
    for field, label in [('platelet', 'Platelet count, 10⁹/L'), ('pf_ratio', 'PaO₂/FiO₂, mmHg'), ('il6_value', 'IL-6'), ('wbc_value', 'White blood cell count, 10⁹/L'), ('bilirubin_value', 'Total bilirubin'), ('creatinine_value', 'Creatinine')]:
        continuous(label, field, 'Not available')
    return output
