# -*- coding: utf-8 -*-


"""
Add ANSI color to the levelname, for nicer printing of log records.

Taken from: https://stackoverflow.com/a/384125
"""


import sys
import string
import re
from copy import copy
import logging
from UItools.shellchrome import RESET_SEQ, COLOR_SEQ, BOLD_SEQ, COLOR, \
                                BG_COLOR, COLORS_I


LVL_I = {
    'WARNING':  COLORS_I.YELLOW,
    'INFO':     COLORS_I.GREEN,
    'DEBUG':    COLORS_I.WHITE,
    'CRITICAL': COLORS_I.YELLOW,
    'ERROR':    COLORS_I.RED}
LVL_COLNAME = {  # For Html.
    'WARNING':  'gold',  # Orange/GoldenRod
    'INFO':     'green',
    'DEBUG':    'gray',
    'CRITICAL': 'orange',
    'ERROR':    'red'}

BASIC_FORMAT = '$LVL%(levelname)s$RESET:$BOLD%(name)s$RESET:%(message)s'
#BASIC_HTML_FORMAT = '<span style="$LVL">%(levelname)s</span>'

class ColoredFormatter(logging.Formatter):
    """Colorize logs using terminal color codes.
    Works in terminal and Jupyter notebook.
    """

    RESET_SEQ = RESET_SEQ
    BOLD_SEQ = BOLD_SEQ
    BG_LVLCOLOR = {lvl: COLOR_SEQ % (40+i) for lvl,i in LVL_I.items()}
    LVLCOLOR    = {lvl: COLOR_SEQ % (30+i) for lvl,i in LVL_I.items()}
    LVLCOLOR.update(ERROR=COLOR['RED'], CRITICAL=BG_COLOR['bgred'])
    COLOR = COLOR
    BG_COLOR = BG_COLOR

    def __init__(self, fmt=None, datefmt=None, style='%'):
        self.fmt_color_template = fmt
        if fmt is not None:
            colored_template = string.Template(fmt)
            fmt = colored_template.substitute(
                            BGLVL=('%(bglvlcol)s' if style=='%' else '{bglvlcol}'),
                            LVL=('%(lvlcol)s' if style=='%' else '{lvlcol}'),
                            RESET=self.RESET_SEQ,
                            BOLD=self.BOLD_SEQ,
                            **self.COLOR, **self.BG_COLOR)
        super(ColoredFormatter, self).__init__(fmt, datefmt, style)


    def format(self, record):
        levelname = record.levelname
        if levelname in self.LVLCOLOR:
            record = copy(record)
            record.__dict__.update(bglvlcol=self.BG_LVLCOLOR[levelname],
                                   lvlcol=self.LVLCOLOR[levelname])
        return super(ColoredFormatter, self).format(record)

    @classmethod
    def install(cls, *loggers, stream=sys.stderr, format=BASIC_FORMAT, **formatter_kwargs):
        sh = logging.StreamHandler(stream)
        sh.setFormatter(cls(format, **formatter_kwargs))
        if not loggers:
            logging.basicConfig(handlers=[sh])
        else:
            for logger in loggers:
                logger.addHandler(sh)

RESET_REGEX = re.compile(r'\$RESET(?![A-Za-z0-9])|\$\{RESET\}')


def template_resetall(fmt, opencodes):
    """Update RESET codes so that they close *all* the opened environments.
    For example, if RED and then BOLD were opened, the next RESET is supposed
    to close both: so it is replaced two RESETs.
    In html, this ensures that the right number of closing </span> are written.
    """
    opening = re.compile(r'|'.join(r'\$'+code+r'(?![A-Za-z0-9_])|\$\{'+code+r'\}' for code in opencodes))
    start = 0
    closed = 0
    it = 0
    while start < len(fmt):
        if it>1000:
            raise RuntimeError('While loop limit exceeded')
        it += 1
        match = RESET_REGEX.search(fmt[start:])
        if match:
            pos = start + match.end()
            matchstart = start + match.start()
            closed += 1
        else:
            pos = matchstart = len(fmt)
        opened = len(opening.findall(fmt, start, pos))
        #logger.debug('$RESET at index=%d, after %d openings', pos, opened)
        if matchstart>0 and fmt[matchstart-1] == '$':
            # the reset is escaped. Do not close now, but memorize the openings
            closed -= (opened + 1)
            start = pos
            #logger.debug('Escaped. opened=%d -> closed=%d; next start=%d', opened, closed, start)
            continue
        while RESET_REGEX.match(fmt[pos:]):
            if closed < opened:
                closed += 1
            pos += RESET_REGEX.match(fmt[pos:]).end()
        #TODO: ignore if already the right number of resets
        fmt = fmt[:pos] + '${RESET}'*(opened-closed) + fmt[pos:]
        start = pos + 8 * max(0, opened-closed)
        #logger.debug('New: "%s" len=%d next start: %d closed: %d', fmt, len(fmt), start, closed)
        closed = 0
    return fmt


convert_html_char = {'<': '&lt;', '>': '&gt;', '&': '&amp;'}
escape_html_char = {'<': r'\<', '>': r'\>', '&': r'\&'}
escape_html_trans = str.maketrans(escape_html_char)
#def escape_html(txt):
#    return txt.translate(escape_html_trans)

#trans_special_html_char = str.maketrans(special_html_char)
html_escaped = re.compile(r'\\(<|>|\\|&)')
html_converted = re.compile(r'(?<!(?<!\\)\\)(<|>|&)')
def to_html(txt):
    txt = html_converted.sub(lambda m: convert_html_char[m.group(1)], txt)
    return html_escaped.sub(lambda m: m.group(1), txt)
#TODO: how to escape quotes in raw strings: e.g. this backslashed quotes print as is: r"thing = \"value\""


def instyle(style: str):
    return ('<span style="%s">' % style).translate(escape_html_trans)


class HtmlColoredFormatter(ColoredFormatter):
    RESET_SEQ = '</span>'.translate(escape_html_trans)
    BOLD_SEQ = instyle("font-weight:bold")
    BG_LVLCOLOR = {lvl: instyle("background:%s" % col) for lvl,col in LVL_COLNAME.items()}
    LVLCOLOR    = {lvl: instyle("color:%s" % col) for lvl,col in LVL_COLNAME.items()}
    LVLCOLOR.update(ERROR=instyle("color:red;font-weight:bold"),
                    CRITICAL=instyle("background:red"))
    COLOR = {col.lower(): instyle("color:%s" % col.lower()) for col in COLOR}
    COLOR.update({col.upper(): instyle("color:%s;font-weight:bold" % col.lower())
                  for col in COLOR})
    BG_COLOR = {('BG' if col.isupper() else 'bg')+col: tag.replace('color:', 'background:')
                for col, tag in COLOR.items()}
    
    def __init__(self, fmt=None, datefmt=None, style='%'):
        if fmt is not None:
            fmt = template_resetall(fmt, ['LVL', 'BGLVL', 'BOLD'] + list(COLOR) + list(BG_COLOR))
        super(HtmlColoredFormatter, self).__init__(fmt, datefmt, style)
            
    
    def format(self, record):
        return to_html(super(HtmlColoredFormatter, self).format(record))

# Pb: if there were width specifications, they aren't visually respected anymore (because ANSI escape codes consume width internally, but not visibly)...
# Better: using newstyle string formatting, then:
#import string
#formatparse = string.Formatter().parse
#
#class newColoredFormatter(logging.Formatter):
#    def __init__(self, fmt=None, datefmt=None, bg=True):
#        self.base = 40 if bg else 30
#        self.colors = {level: code+self.base for level,code in COLORS.items()}
#        
#        for txt, fieldname, format_spec, converter in formatparse(fmt):
#            if fieldname == 'levelname':
#                break
#        else:
#            logging.Formatter.__init__(self, fmt, datefmt, style='{')
#            return
#
#        levelname_fmtstr = '{' + field_name
#        if converter: levelname_fmtstr += '!' + converter
#        if format_spec: levelname_fmtstr += ':' + format_spec
#        levelname_fmtstr += '}'
#
#        levelname_fmtstr = .. + levelname_fmtstr + ..
#        fmt = 
#        logging.Formatter.__init__(self, fmt, datefmt, style='{')
#
#    def format(self, record):
#        levelname = record.levelname
#        if levelname in COLORS:
#            record = copy(record)
#            record.levelname = (COLOR_SEQ % (self.base + COLORS[levelname])
#                               + levelname + RESET_SEQ)
#        return logging.Formatter.format(self, record)
