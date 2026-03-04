import re

def test_parse(item):
    # Mimic JS split logic
    tf_name = item.split(' (')[0].strip()
    
    prob = 0
    effect = '?'
    l_dir = '?'
    b_dir = '?'
    l_cons = ''
    b_cons = ''

    # Regex patterns
    # Full: TF (0.90, + [+(0.98),+(0.95)])
    full_pattern = r'\(([\d\.]+),\s*([\+\-])\s*\[([\+\-])\(([\d\.]+)\),([\+\-])\(([\d\.]+)\)\]\)'
    # Simple: TF (0.90, + [+,+])
    simple_pattern = r'\(([\d\.]+),\s*([\+\-])\s*\[([\+\-]),([\+\-])\]\)'
    # ML: TF (0.90)
    ml_pattern = r'\(([\d\.]+)\)'

    m = re.search(full_pattern, item)
    if m:
        prob, effect, l_dir, l_cons, b_dir, b_cons = m.groups()
        print(f"MATCH FULL: name='{tf_name}' prob={prob} effect={effect} L={l_dir}({l_cons}) B={b_dir}({b_cons})")
        return
    
    m = re.search(simple_pattern, item)
    if m:
        prob, effect, l_dir, b_dir = m.groups()
        print(f"MATCH SIMPLE: name='{tf_name}' prob={prob} effect={effect} L={l_dir} B={b_dir}")
        return

    m = re.search(ml_pattern, item)
    if m:
        prob = m.group(1)
        print(f"MATCH ML: name='{tf_name}' prob={prob}")
        return
    
    print(f"FAIL: {item}")

# Test Cases
test_cases = [
    "STAT5B (0.85, + [+(0.98),+(0.95)])",
    "ZNF169 (0.72, - [-,+])",
    "KLF1 (0.99)",
    "ZFP64 (1.00), ",
    "FOXP1 (0.82, + [+(0.90),-(0.10)])"
]

for tc in test_cases:
    test_parse(tc)
