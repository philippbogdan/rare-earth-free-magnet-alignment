"""Markdown rendering helpers for generated documents."""


def fmt(v, nd=3):
    if v is None:
        return '–'
    if isinstance(v, float):
        return f'{v:.{nd}f}'
    return str(v)


def table(headers, rows):
    out = ['| ' + ' | '.join(headers) + ' |', '|' + '|'.join('---' for _ in headers) + '|']
    for r in rows:
        out.append('| ' + ' | '.join(str(x).replace('|', '/').replace('\n', ' ') for x in r) + ' |')
    return '\n'.join(out)


def fill(template, values):
    """Replace {{key}} placeholders; fail loudly on any unfilled placeholder."""
    out = template
    for k, v in values.items():
        out = out.replace('{{' + k + '}}', str(v))
    if '{{' in out:
        i = out.index('{{')
        raise KeyError('unfilled placeholder: ' + out[i:i + 40])
    return out
