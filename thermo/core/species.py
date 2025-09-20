import re

def parse_CH(s: str):
    m = re.fullmatch(r'C(\d*)H(\d+)', s)
    if not m: return None, None
    C = int(m.group(1)) if m.group(1) else 1
    H = int(m.group(2))
    return C, H
