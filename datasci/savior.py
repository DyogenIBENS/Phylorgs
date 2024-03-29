#!/usr/bin/env python3
# -*- coding: utf-8 -*-


import sys
import os
import os.path as op
from datetime import datetime as dt
from contextlib import contextmanager
from io import StringIO, BytesIO
from pprint import pprint, pformat
import base64
import matplotlib.pyplot as plt
import pandas as pd
from pandas.io.formats.style import Styler
from markdown import Markdown
from IPython.display import display_html, display_markdown
import traceback
import warnings
import logging
logger = logging.getLogger(__name__)
try:
    from UItools.colorlog import ColoredFormatter, HtmlColoredFormatter, BASIC_FORMAT
except ImportError:
    ColoredFormatter = HtmlColoredFormatter = logging.Formatter
    BASIC_FORMAT = logging.BASIC_FORMAT


trans_special_html_char = str.maketrans({'<': '&lt;', '>': '&gt;', '&': '&amp;'})


def format_p(txt):
    return '<p>' + txt + '</p>\n'

def format_pre(txt):
    """Display unformatted text (fixed-width)"""
    return '<pre>\n' + txt + '\n</pre>\n'

def format_mkd(txt):
    pass


def format_fig(fig, format='svg', **kwargs):
    with StringIO() as out:
        fig.savefig(out, format=format, **kwargs)

def figure(*content):
    return '<figure>\n' + '\n'.join(content) + '\n</figure>\n'

def format_fig_embed(fig, format='png', bare=False, **kwargs):
    kwargs = {'bbox_inches': 'tight', **kwargs}
    open_out = StringIO if format=='svg' else BytesIO
    with open_out() as out:
        fig.savefig(out, format=format, **kwargs)
        imgdata = out.getvalue()

    # image/svg+xml
    if format == 'svg':
        #TODO: the <?xml> opening tag should be stripped out.
        #          and the '<!DOCTYPE svg PUBLIC' also.
        #          We should use an XML parser for that, but we will assume
        #          That the 3 first lines are to be removed.
        imglines = imgdata.split('\n', 10)  # Should be within the 10 first lines?
        for i, line in enumerate(imglines):
            if line.startswith('<!-- ') or line.startswith('<svg '):
                break
            assert line.startswith('<?xml ') or line.startswith('<!DOCTYPE ') or line.startswith(' ')

        img = imglines[i:]
        #<object type="image/svg+xml" data="image.svg"></object>

    else:
        imgdata = base64.encodebytes(imgdata).decode()
        img = ['<img src="data:image/%s;base64,%s" />' % (format, imgdata)]
    return img[0] if bare else figure(*img)
    # alt=
    # style="width:500px;height:600px;" || width="500" height="600"
    # <figcaption>

def format_fig_extern(fig, filename, bare=False, **kwargs):
    kwargs = {'bbox_inches': 'tight', **kwargs}
    fig.savefig(filename, **kwargs)
    format = op.splitext(filename)[1].lstrip('.')
    img = '<img src="%s" />' % filename
    return img if bare else figure(img)
# For interactive svg, I should use <iframe> or <object>.

# See help of IPython.display.display:
#  - `_repr_html_`: return raw HTML as a string, or a tuple (see below).
#  - `_repr_json_`: return a JSONable dict, or a tuple (see below).
#  - `_repr_jpeg_`: return raw JPEG data, or a tuple (see below).
#  - `_repr_png_`: return raw PNG data, or a tuple (see below).
#  - `_repr_svg_`: return raw SVG data as a string, or a tuple (see below).
#  - `_repr_latex_`: return LaTeX commands in a string surrounded by "$",
#                    or a tuple (see below).
#

def files_read_lines(filenames):
    #TODO: when reading a javascript or css file, remove comments.
    lines = []
    for filename in filenames:
        with open(filename) as f:
            lines.extend(line.rstrip() for line in f)
    return lines


from abc import ABC, abstractmethod

class Report(ABC):
    """Base class for HtmlReport and CellReport. Do not instantiate directly."""
    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, tb):
        logger.debug('Exited context: %r %r %r', exc_type, exc_value, tb)

    @abstractmethod
    def raw(self, txt):
        """Write raw string"""
    # Facing problem here: to be able to pass Report to the file arg of print(),
    # I'd need the `write` method to behave as `print()`

    @abstractmethod
    def write(self, txt):
        """Write string. Implement it so that it can be used by print calls."""

    @abstractmethod
    def pprint(self, *args):
        """pretty print"""

    @abstractmethod
    def show(self, fig=None, **kwargs):
        """Show figure"""
        #return plt.gcf() if fig is None else fig

    @abstractmethod
    def html(self, obj):
        """Write the html representation of an object"""
        #return obj.to_html()

    @abstractmethod
    def mkd(self, *strings):
        pass

    def print(self, *args, sep=' ', end='\n'):
        """Write unformatted text"""
        self.write(sep.join(str(a) for a in args) + end)

    def output(self, *args):
        """Automatically dispatch figures, dataframes and text (as printed)."""
        for arg in args:
            if isinstance(arg, (plt.Figure)):  # TODO: seaborn.axisgrid.FacetGrid
                self.show(arg)
            elif isinstance(arg, (pd.DataFrame, Styler)):
                self.html(arg)
            else:
                self.print(arg)

    def display(self, *args):
        """like output(), but fall back to markdown text."""
        for arg in args:
            if isinstance(arg, (plt.Figure)):  # TODO: seaborn.axisgrid.FacetGrid
                self.show(arg)
            elif isinstance(arg, (pd.DataFrame, Styler)):
                self.html(arg)
            else:
                self.mkd(arg)


css_dark_style = '''body {
background-color: black;
color: #E5E5E5;
}
'''

# For default Jupyter Notebook styles,
# see /static/style/style.min.css?v=29c09309dd70e7fe93378815e5f022ae
# or /Users/grant/Sites/jupyter/notebook/notebook/static/notebook/less/renderedhtml.less
class HtmlReport(Report):

    @classmethod
    def get_default_styles(cls):
        default_styles = ['''
            body {
              max-width: 940px;
              margin-right: auto;
              margin-left: auto;
            }
            pre {
              display: inline-block; /* ensures that nested pre lines do not have to much vertical spacing */
              white-space: pre-wrap;
              width: max(900px, 80vw);
            }
            pre.log {
              margin-bottom: 0;
              margin-top: 0;
            }
            table.dataframe {
              margin: 0 auto;
              border: none;
              border-collapse: collapse;
              border-spacing: 0;
              /*color: black;*/
              table-layout: fixed;
            }
            table > *, table.dataframe > * {
              font-size: x-small;
            }
            .dataframe thead {
              border-bottom: 1px solid black;
              vertical-align: bottom;
            }
            .dataframe tr, .dataframe th, .dataframe td {
              text-align: right;
              vertical-align: middle;
              padding: 0.5em 0.5em;
              line-height: normal;
              white-space: normal;
              max-width: none;
              border: none;
            }
            .dataframe tbody tr:nth-child(2n+1) {
              background: #7f7f7f33;
            }
            .dataframe tbody tr th:only-of-type {
              vertical-align: middle;
            }
            div>table[class="dataframe"]+p {
              font-family: monospace;
              font-size: small;
              text-align: center;
            }
''']
        default_styles += ['''
            .column{0} {{
              float: left;
              width: {0}%;
            }}'''.format(p) for p in range(10,100,10)]
        default_styles += ['''
            /* Clear floats after the columns */
            .columns:after {
              content: "";
              display: table;
              clear: both;
            }
            ''']

        return default_styles

    def __init__(self, filename=None, mode='w', figformat='png', title=None,
                 css=None, scripts=None, scripts_embedded=None,
                 mathjax=True, style=None, metas=None,
                 postscripts=None, postscripts_embedded=None, external_content=False):
        self.filename = filename
        self.mode = mode
        self.figformat = figformat
        self.printing = False  # Whether a <pre> is already opened
        if css is None: css = []
        if scripts is None: scripts = []
        if scripts_embedded is None: scripts_embedded = []
        if postscripts is None: postscripts = []
        if postscripts_embedded is None: postscripts_embedded = []
        for argname, arg in [('css', css), ('scripts', scripts),
                             ('scripts_embedded', scripts_embedded),
                             ('postscripts', postscripts),
                             ('postscripts_embedded', postscripts_embedded)]:
            if not isinstance(arg, (list, tuple)):
                raise ValueError(argname+' argument should be a list/tuple.')
        if metas is None: metas = ['charset="UTF-8"']
        self._mkd = Markdown(extensions=['fenced_code', 'attr_list', 'tables', 'sane_lists'])  # Use a single instance for the whole report

        self.external_content = external_content
        if external_content:
            # Dir content relative to the Html file.
            self.content_dir = op.basename(op.splitext(self.filename)[0])
            self.figcounter = 0
            def format_fig(fig, format='png', bare=False, **kwargs):
                self.figcounter += 1
                # if fig._suptitle is not None:
                figfile = op.join(self.content_dir, '%03d.%s' %(self.figcounter, format))
                return format_fig_extern(fig, figfile, bare=bare, **kwargs)
            self.format_fig = format_fig
            # TODO: external CSS as well.
        else:
            self.format_fig = format_fig_embed

        default_styles = False
        self.default_css = None  # External file (embed by default).
        if style is not False:
            # Default styling
            default_styles = self.get_default_styles()
            if external_content:
                self.default_css = op.join(self.content_dir, 'default.css')
                css.insert(0, self.default_css)

        styles = []
        if style and style is not True:
            styles.append(style)

        self.head = '\n'.join(
            ['<!DOCTYPE html>', '<html>', '<head>',
             '<title>%s</title>' if title else ''] +
            ['<meta %s />' % m for m in metas] +
            (['<style>\n%s\n</style>' % '\n'.join(default_styles)] if default_styles and not self.default_css else []) +
            ['<link rel="stylesheet" type="text/css" href="%s" />' % c for c in css] +
            (['<script type="text/javascript" async src="https://cdnjs.cloudflare.com/ajax/libs/mathjax/2.7.5/latest.js?config=TeX-MML-AM_CHTML"></script>'] if mathjax else []) +
            ['<script type="text/javascript" src="%s" async></script>' % s for s in scripts] +
            (['<script type="text/javascript">', *files_read_lines(scripts_embedded), '</script>'] if scripts_embedded else []) +
            ['<style>\n%s\n</style>' % '\n'.join(styles),
            '</head>',
            '<body>\n'])
        self.foot = '\n'.join(['<script type="text/javascript" src="%s"></script>' % s
                               for s in postscripts] +
                              (['<script>', *files_read_lines(postscripts_embedded), '</script>']
                               if postscripts_embedded else []) +
                              ['</body>', '</html>\n'])
        self.closed = True
        self.begun = []  # opened html tags inside body.

    def setlogfmt(self, fmt=None, datefmt=None, style='%'):
        self.loghandler.setFormatter(
                HtmlColoredFormatter(r'\<pre class="log"\>' + fmt + r'\</pre\>',
                                     datefmt, style))

    def open(self, mode):
        self.handle = StringIO() if self.filename is None else open(self.filename, mode)
        self.loghandler = logging.StreamHandler(self.handle)
        self.setlogfmt(BASIC_FORMAT)
        self.closed = False
        if self.external_content:
            # Filename should not be None
            os.makedirs(self.content_dir, exist_ok=True)
            if self.default_css:
                with open(self.default_css, 'w') as css_f:
                    css_f.write('\n'.join(self.get_default_styles()) + '\n')

        logger.debug('Opened %s(filename=%r, mode=%r) -> closed=%s, handle=%r',
                     self.__class__.__name__, self.filename, mode, self.closed,
                     self.handle)

    def close(self):
        logger.debug('Closing </body>, then file handle.')
        self.loghandler.close()
        logging.root.removeHandler(self.loghandler)
        for logr in logging.root.manager.loggerDict.values():
            try:
                logr.removeHandler(self.loghandler)
            except AttributeError:
                # Some stored loggers have class PlaceHolder.
                pass
        footer_date = '\n<br /><br /><br /><br />\n<footer><p>Generated: %Y/%m/%d %H:%M:%S</p></footer>' 
        self.handle.write(dt.now().strftime(footer_date) + self.foot)
        if self.filename is None:
            display_html(self.handle.getvalue(), raw=True)
        self.handle.close()
        self.closed = True

    def flush(self):
        self.handle.flush()

    def __enter__(self):
        self.open(self.mode)
        if self.mode == 'w':
            logger.debug('W %d bytes (head)', self.handle.write(self.head))
        logger.debug('Entered.')
        return self

    def __exit__(self, exc_type, exc_value, tb):
        if self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        if exc_type is not None:
            fmtted_exc = [line.translate(trans_special_html_char) for line in
                          traceback.format_exception(exc_type, exc_value, tb)]
            print('\n'.join(fmtted_exc) + '\n', file=sys.stderr)
            err_name = exc_type.__name__.strip("'")
            colored_etype = '<span style="color:red;">%s</span>' % err_name
            fmtted_exc[-1] = fmtted_exc[-1].replace(err_name, colored_etype)
            self.handle.write('<pre><hr style="margin:0;padding:0;border:1px dashed red;" />\n' + colored_etype +'\n')
            self.handle.write('\n'.join(fmtted_exc) + '\n</pre>\n')

        self.close()
        super().__exit__(exc_type, exc_value, tb)

    def raw(self, txt):
        """Write raw string to file (should be valid Html)."""
        #logger.debug('W %d bytes (body)', )
        self.handle.write(txt)

    def write(self, txt):
        """Write unformatted text (<pre> environment)"""
        if not self.printing:
            self.handle.write('<pre>\n')
            self.printing = True
        #logger.debug('W %d bytes (body)', )
        self.handle.write(txt.translate(trans_special_html_char))

    def echo(self, *args, **kwargs):
        """Print to the file and to stdout in the same time."""
        print(*args, **kwargs)
        self.print(*args, **kwargs)

    def pprint(self, *args):
        """pretty print (in <pre> environment)"""
        self.print(*(pformat(a) for a in args), sep='\n')

    def show(self, fig=None, **kwargs):
        """Show figure"""
        if self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        if fig is None:
            if not plt.get_fignums():
                raise RuntimeError('No current figure!')
            fig = plt.gcf()
        self.handle.write(self.format_fig(fig, self.figformat, **kwargs))

    def html(self, obj):
        if self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        obj_html = obj._repr_html_()
        if obj_html is None:
            raise ValueError('method _repr_html_() returned None.')
        self.handle.write(obj_html + '\n')

    def mkd(self, *strings):
        if self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        self.handle.write('\n'.join(self._mkd.convert(s) for s in strings) + '\n')

    def begin(self, tag, class_=None, id_=None, attrs=None, print_off=True):
        """Open an html tag"""
        if print_off and self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        tagstr = tag
        if class_:
            tagstr += (' class="%s"' % class_)
        if id_:
            tagstr += (' id="%s"' % id_)
        if attrs:
            tagstr += ' ' + ' '.join('%s="%s"' %(k,v) for k,v in attrs.items())
        self.handle.write('<%s>\n' % tagstr)
        self.begun.append((tag, class_, id_, attrs))

    def end(self, tag=None, class_=None, id_=None, print_off=True):
        if print_off and self.printing:
            self.handle.write('</pre>\n')
            self.printing = False
        if not self.begun:
            raise ValueError('Trying to close tag but none were opened.')
        for spec_name, spec, prev_spec in zip(('tag', 'class', 'id'),
                                              (tag, class_, id_),
                                              self.begun[-1][:3]):
            if spec is not None and spec != prev_spec:
                raise ValueError('Mismatching closing %s %r (last opened %r)'
                                 % (spec_name, spec, prev_spec))
        self.handle.write('\n</%s>\n' % tag)
        self.begun.pop()

    def ends(self, *tags):
        """Close html tags, if matching opened ones. Otherwise raise ValueError."""
        for tag in tags:
            self.end(tag)

    def end_all(self):
        while self.begun:
            self.handle.write('\n</%s>' % self.begun.pop()[0])

    def columns(self):
        self.begin('div', class_='columns')

    def column(self, p=50):
        self.begin('div', class_='column%d' % p)

    def end_column(self, p=None):
        if p is None:
            begun, begun_class = self.begun[-1][:2]
            if not begun_class.startswith('column') or begun_class == 'columns':
                raise ValueError('Mismatching closing column (last opened "%s .%s")'
                                 % (begun, begun_class))
        else:
            begun_class = 'column%2d' % p
        self.end('div', begun_class)
        self.handle.write('<!-- end %s -->\n' % begun_class)

    def end_columns(self):
        self.end('div', 'columns')
        self.handle.write('<!-- end columns -->\n<br />\n')

    # Haven't found elegant way to use `cols` with or without `with ...` using @contextmanager
    @contextmanager
    def col(self, p=50):
        self.column(p)
        yield
        self.end_column(p)

    @contextmanager
    def cols(self):
        self.columns()
        yield
        self.end_columns()


class CellReport(Report):
    #def __init__(self):
    #    pass

    def __enter__(self):
        self.loghandler = logging.StreamHandler(sys.stdout)
        self.loghandler.setFormatter(ColoredFormatter(fmt=BASIC_FORMAT))
        return self
    def __exit__(self):
        self.loghandler.close()
        logging.root.removeHandler(self.loghandler)
        for logr in logging.root.manager.loggerDict.values():
            try:
                logr.removeHandler(self.loghandler)
            except AttributeError:
                # Some stored loggers have class PlaceHolder.
                pass

    def raw(self, txt):
        """Write raw string"""
        print(txt, end='')

    def write(self, txt):
        #stdout.write(txt)
        print(txt, end='')

    def print(self, *args, sep=' ', end='\n'):
        """Write unformatted text"""
        #TODO: turn off pretty printing
        print(*args, sep=sep, end=end)

    def pprint(self, *args):
        """pretty print"""
        for arg in args:
            pprint(arg, sep=sep, end=end)

    def show(self, fig=None, **kwargs):
        """Show figure"""
        if fig is None:
            plt.show(**kwargs)
        else:
            fig.show()  # Not sure this works as expected.

    def html(self, obj):
        obj_html = obj._repr_html_()
        if obj_html is None:
            raise ValueError('method _repr_html_() returned None.')
        display_html(obj_html, raw=True)

    def mkd(self, *strings):
        display_markdown(*strings)


class DoubleReport(HtmlReport, CellReport):
    """Double displaying: in an HTML document, and in the cell output."""
    # __init__ inherited from HtmlReport
    def write(self, txt):
        """Write raw string"""
        HtmlReport.write(self, txt)  # Will call both parents only if HtmlReport uses a super() call...
        CellReport.write(self, txt)

    def print(self, *args, sep=' ', end='\n'):
        """Write unformatted text"""
        HtmlReport.print(self, *args, sep=sep, end=end)
        CellReport.print(self, *args, sep=sep, end=end)

    def pprint(self, *args):
        """pretty print"""
        HtmlReport.pprint(self, *args)
        CellReport.pprint(self, *args)

    def show(self, fig=None, **kwargs):
        """Show figure"""
        HtmlReport.show(self, fig, **kwargs)
        CellReport.show(self, fig, **kwargs)

    def html(self, obj):
        HtmlReport.html(self, obj)
        CellReport.html(self, obj)

    def mkd(self, *strings):
        HtmlReport.mkd(self, *strings)
        CellReport.html(self, *strings)


def generate_slideshow(hr: HtmlReport, nmax: int=100, **attrs):
    """Needs JS and CSS from `atavistic-doc-tools/myslideshow`"""
    hr.begin('div', class_="slideshow", id_=attrs.pop('id_', None), attrs=attrs)
    hr.begin('ul')
    try:
        for i in range(nmax):
            hr.begin('li')
            try:
                yield i
                # Example:
                #plt.plot(np.random.random(5), color=color)
                #hr.show(bare=True)
                #plt.close()
            # Capture a break:
            except StopIteration as err:
                # Unsure but this might be closing hr.
                break
            finally:
                hr.end('li')
    except GeneratorExit:
        hr.end('ul')
        hr.end('div', class_='slideshow')
        if i==nmax-1:
            warnings.warn('Reached max number of slides. This iterator should be explicitly stopped before.')

def slideshow_generator(hr: HtmlReport, **gkwargs):
    # Mimic the implementation of UItools.jupytertricks.make_folds()
    def slideshow_iter(iterable):
        for elem, slide in zip(iterable, generate_slideshow(hr, **gkwargs)):
            yield elem
    return slideshow_iter


_LOGGING_STYLECLASS_CODE = {constructor[0]: codestring for codestring, constructor in logging._STYLES.items()}


@contextmanager
def reroute_loggers(hr: HtmlReport, *loggers):
    """Update the handlers of loggers to write into the HtmlReport"""
    if len(loggers)==1 and isinstance(loggers[0], (tuple, list)):
        loggers, = loggers  # Unpack.
    logger_states = {}
    with hr as hr_entered:
        for logr in loggers:
            #if logger.hasHandlers():
            logger_states[logr.name] = logr.handlers
            logr.handlers = [hr_entered.loghandler]
            # Preserve the log format:
            if logger_states[logr.name]:
                orig_formatter = logger_states[logr.name][0].formatter
                try:
                    orig_fmt = orig_formatter.fmt_color_template  # class UItools.colorlog.ColoredFormatter
                except AttributeError:
                    orig_fmt = orig_formatter._fmt
                orig_style = _LOGGING_STYLECLASS_CODE[orig_formatter._style.__class__]
                hr_entered.setlogfmt(orig_fmt, orig_formatter.datefmt, orig_style)
        yield hr_entered
        for logr in loggers:
            logr.handlers = logger_states.pop(logr.name)


