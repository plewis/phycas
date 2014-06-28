import sys, os
from cStringIO import StringIO
import textwrap
from phycas import getDefaultOutFilter, OutFilter, help_double_space, current_double_space, current_follows_help
from phycas.pdfgen import PDFGenerator
import phycas.readnexus as readnexus
from phycas.utilities.CommonFunctions import getDefaultOutputStream
from phycas.utilities.io import FileFormats, DataSource, TreeCollection
from phycas import OutputFilter, name_of
import copy
try:
    _s = set()
except :
    from sets import Set as set

_fixed_terminal_width = None
_use_instance_names = True

def setFixedTerminalWidth(w):
    global _fixed_terminal_width
    _fixed_terminal_width = w

###############################################################################
def ttysize():
    global _fixed_terminal_width
    if _fixed_terminal_width:
        return (1,_fixed_terminal_width)
    if os.name == 'nt':
        # From http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/440694
        from ctypes import windll, create_string_buffer
        # stdin handle is -10
        # stdout handle is -11
        # stderr handle is -12

        h = windll.kernel32.GetStdHandle(-12)
        csbi = create_string_buffer(22)
        res = windll.kernel32.GetConsoleScreenBufferInfo(h, csbi)

        if res:
            import struct
            (bufx, bufy, curx, cury, wattr,
             left, top, right, bottom, maxx, maxy) = struct.unpack("hhhhHhhhhhh", csbi.raw)
            sizex = right - left + 1
            sizey = bottom - top + 1
            return sizey, sizex
        else:
            return None
    else:
        # From http://mail.python.org/pipermail/python-list/2006-February/365594.html
        try:
            import fcntl
            import struct
            import termios
            buf = 'abcdefgh'
            buf = fcntl.ioctl(0, termios.TIOCGWINSZ, buf)
            row, col, rpx, cpx = struct.unpack('hhhh', buf)
            return row, col
        except:
            return None
###############################################################################

def str_value_for_user(value):
    try:
        return value._brief_str()
    except:
        if isinstance(value, str):
            return repr(value)
        return str(value)

def _escape_for_latex(s):
    return str(s).replace('_', '\_').replace('<','$<$').replace('>','$>$')

#def public(mod=None):
#    "Prints a summary of the names of the globally available objects"
#    l = eval('globals().keys()')
#    l = [i for i in l if not i.startswith("_")]
#    print two_column_str(l)

def two_column_str(p):
    p.sort()
    if len(p) % 2 == 0:
        half = len(p)/2     # p has an even number of elements
    else:
        half = len(p)/2 + 1 # p has an odd number of elements
    ait = iter(p[:half])
    bit = iter(p[half:])
    o = []
    done = False
    while not done:
        a, b = ("", "")
        try:
            a = ait.next()
        except StopIteration:
            done = True
        try:
            b = bit.next()
        except StopIteration:
            done = True
        r = PhycasTablePrinter.format_help(a, b)
        o.append(r)
    return "\n".join(o)

class PhycasHelp(object):
    _phycas_cmd_classes = set()
    def __str__(self):
        PhycasTablePrinter._reset_term_width()
        n = [i.__name__.lower() for i in PhycasHelp._phycas_cmd_classes if not i.hidden()]
        t = two_column_str(n)
        return """Phycas Help

For Python Help use "python_help()"

Commands are invoked by following the name by () and then
hitting the RETURN key. Thus, to invoke the sumt command use:

sumt()

Commands (and almost everything else in python) are case-sensitive -- so
"Sumt" is _not_ the same thing as "sumt" In general, you should use the
lower case versions of the phycas command names.

The currently implemented Phycas commands are:

%s

Use <command_name>.help to see the detailed help for each command. So,

sumt.help

will display the help information for the sumt command object.
""" % t

    def _print_help(self, a, level=0):
        if isinstance(a, PhycasCommand):
            a.help()
        elif isinstance(a, PhycasCommandOutputOptions):
            print (str(a))
        elif a is None:
            print("None is the Python object used to represent the concept 'undefined' or 'not applicable'")
        elif a is True:
            print("True is the Python's representation of the boolean condition 'true'.")
        elif a is False:
            print("False is the Python's representation of the boolean condition 'false'.")
        elif isinstance(a, int) or isinstance(a, long):
            print("The integer", str(a))
        elif isinstance(a, float):
            print("The real number", str(a))
        elif isinstance(a, str):
            print("The string %s\nIf you would like to see the methods available for a string use the command 'help(str)'" % repr(a))
        elif isinstance(a, list):
            print("A python list. Use [number] to access elements in the list. For instance, x[0] is the first element of list x.\nUse 'help(list)' to see a list of available methods")
        elif isinstance(a, list):
            print("A python tuple. Use [number] to access elements in the list. For instance, x[0] is the first element of tuple x.\nUse 'help(tuple)' to see a list of available methods")
        elif isinstance(a, dict):
            print("A python dict. Use [key] to access elements in the list. For instance, x[2] returns the value associated with the key 2 (if there is such an element in the dictionary).\nUse 'help(dict)' to see a list of available methods")
        else:
            import pydoc
            d = pydoc.getdoc(a)
            if isinstance(a, type):
                m = pydoc.allmethods(a)
                pub_methods = [k for k in m.iterkeys() if not k.startswith("_")]
                if pub_methods or d:
                    if level == 0:
                        print("\nThis is a class or python type")
                    else:
                        print("Information about this type:")
                    if d:
                        print(d)
                    if pub_methods:
                        print("\nThe following public methods are available:\n%s" % '\n'.join(pub_methods))
                else:
                    print("\nThis is an undocumented class or python type without public methods")
            else:
                print("\n%s\n\nAn instance of type %s.\n" % (repr(a), a.__class__.__name__))
                if callable(a):
                    if d:
                        print(d)
                else:
                    self._print_help(a.__class__, level+1)

    def __call__(self, *args):
        """Returns the result of __str__, so that users can omit the () on the
        and simply type "cmd.help" to see the help message."""
        if len(args) == 0:
            print(str(self))
        for a in args:
            print ("Help entry for %s:" % repr(a))
            self._print_help(a)
    def __repr__(self):
        return str(self)

phycas_help = PhycasHelp()

class ExistingFileBehavior(object):
    pass

class ReplaceExistingFileBehavior(ExistingFileBehavior):
    def __str__(self):
        return "REPLACE"

class AppendExistingFileBehavior(ExistingFileBehavior):
    def __str__(self):
        return "APPEND"

class AddNumberExistingFileBehavior(ExistingFileBehavior):
    def __str__(self):
        return "ADD_NUMBER"

class SkipExistingFileBehavior(ExistingFileBehavior):
    def __str__(self):
        return "SKIP"

REPLACE = ReplaceExistingFileBehavior()
APPEND = AppendExistingFileBehavior()
ADD_NUMBER = AddNumberExistingFileBehavior()
SKIP = SkipExistingFileBehavior()

_opt_name_help_len = 30
_opt_val_help_len = 19

class PhycasTablePrinter:
    _help_str_wrapper = textwrap.TextWrapper()
    _text_wrapper = textwrap.TextWrapper()
    _help_fmt_str = "%%-%ds %%s" % (_opt_name_help_len)

    def _reset_term_width():
        tty_sz = ttysize()
        if tty_sz:
            PhycasTablePrinter._set_terminal_width(tty_sz[1])
    _reset_term_width = staticmethod(_reset_term_width)

    def _set_terminal_width(w):
        global _opt_val_help_len, _opt_name_help_len
        x =  _opt_name_help_len + 1
        new_width = max(x + _opt_name_help_len + 1, w)
        PhycasTablePrinter._help_str_wrapper.width = new_width
        PhycasTablePrinter._help_str_wrapper.subsequent_indent = " "*x
        PhycasTablePrinter._help_str_wrapper.initial_indent = " "*x
        PhycasTablePrinter._text_wrapper.width = new_width
        PhycasTablePrinter._text_wrapper.subsequent_indent = "    "
        PhycasTablePrinter._text_wrapper.initial_indent = ""
    _set_terminal_width = staticmethod(_set_terminal_width)

    def format_help(name, helpStr):
        global _opt_val_help_len, _opt_name_help_len
        hs = "\n".join(PhycasTablePrinter._help_str_wrapper.wrap(str(helpStr)))
        w = _opt_name_help_len + 1
        return PhycasTablePrinter._help_fmt_str % (name, hs[w:])
    format_help = staticmethod(format_help)

    def get_help_divider():
        global _opt_val_help_len, _opt_name_help_len
        l = PhycasTablePrinter._help_str_wrapper.width - (_opt_name_help_len + 1)
        return " ".join(["="*_opt_name_help_len, "="*l])
    get_help_divider = staticmethod(get_help_divider)
PhycasTablePrinter._set_terminal_width(80)
class PhycasOutput(object):
    def __init__(self):
        self._is_active = True
    def _silence(self):
        self._is_active = False
    def _activate(self):
        self._is_active = True
    def __bool__(self):
        return self._is_active

class FileOutputSpec(PhycasOutput):
    def __init__(self, prefix="", help_str="", filename=None):
        PhycasOutput.__init__(self)
        self.mode = ADD_NUMBER
        self.filename = filename
        self.prefix = prefix
        if filename and prefix:
            self.prefix = None
        self._help_str = help_str
        self._valid_formats = []
        self._options = None
        self._opened_filename = None
        self._opened_file = None
        self._prexisting = False
        self._opened_file_is_log_in = None # reference of outputstream that is mirroring output to the file

    def __deepcopy__(self, memo):
        c = memo.get(self)
        if not c is None:
            return c
        o = self.__class__(copy.deepcopy(self.prefix, memo),
                           copy.deepcopy(self._help_str, memo),
                           copy.deepcopy(self.filename, memo))
        opts = self._options
        c = {}
        if opts:
            c = opts._current
            o._options = PhycasCmdOpts(o, opts._optionsInOrder)
            oo = o._options
            for k, v in c.iteritems():
                oo._set_opt(k, v)
        for k, v in self.__dict__.iteritems():
            if (k not in c) and (k not in ["_options", "_opened_filename", "_opened_file", "_opened_file_is_log_in"]):
                o.__dict__[k] = copy.deepcopy(v)
        memo[self] = o
        return o

    def __bool__(self):
        return bool(self._is_active and (self.filename or self.prefix))

    def _help_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = ""
        opts_help = [PhycasTablePrinter.format_help("%s" %pref, self._help_str)]
        if self._getFilename():
            prefix_str = PhycasTablePrinter.format_help("%sprefix" % dpref, "file prefix (appropriate suffix will be added)")
            fn_str = PhycasTablePrinter.format_help("%sfilename" % dpref, "The full file name.")
            mode_str = PhycasTablePrinter.format_help("%smode" % dpref, 'Controls the behavior when the file is present. Valid settings are %s. SKIP results in no file being created, regardless of whether the file name already exists (use this to avoid creating a specific type of output file). ADD_NUMBER indicates that a number will be added to the end of the file name (or prefix) to make the name unique' % self._getValidModeNames())
            if help_double_space:
                prefix_str = '\n%s' % prefix_str
                fn_str = "\n%s" % fn_str
                mode_str = "\n%s" % mode_str
            opts_help.append(prefix_str)
            opts_help.append(fn_str)
            opts_help.append(mode_str)
            if len(self._valid_formats) > 1:
                opts_help.append(PhycasTablePrinter.format_help("%s.format" % dpref, 'Format of the file valid settings are %s' % self._getValidFormatNames()))
        if self._options:
            opts_help.extend(self._options._help_str_list(pref))
        return opts_help

    def _write_latex_item(self, out, pref="", cmdname=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called
        more recently than the last terminal width change)."""
        if pref:
            dpref = _escape_for_latex(pref + ".")
        else:
            dpref = ""

        out.write("\index{%s!%s}\item[\\bftt  %s] %s\n" % (cmdname,  _escape_for_latex(pref), _escape_for_latex(pref), _escape_for_latex(self._help_str)))
        fmt_def = " (default :%s)"
        d = self.prefix and (fmt_def % _escape_for_latex(repr(self.filename))) or ""
        out.write("\index{%s!%sprefix}\item[\\bftt  %sprefix] %s%s\n" % (cmdname, dpref, dpref, "file prefix (appropriate suffix will be added)", d))
        d = self.filename and (fmt_def % _escape_for_latex(repr(self.filename))) or ""
        out.write("\index{%s!%sfilename}\item[\\bftt  %sfilename] %s%s\n" % (cmdname, dpref, dpref, "The full file name. Specifying this field preempts `prefix` setting.", d))
        d = _escape_for_latex(self.mode)
        out.write("\index{%s!%smode}\item[\\bftt  %smode] %s%s\n" % (cmdname, dpref, dpref,
                                                   _escape_for_latex('Controls the behavior when the file is present. Valid settings are %s. SKIP results in no file being created, regardless of whether the file name already exists (use this to avoid creating a specific type of output file). ADD_NUMBER indicates that a number will be added to the end of the file name (or prefix) to make the name unique' % self._getValidModeNames()),
                                                   d))
        if len(self._valid_formats) > 1:
            d = _escape_for_latex(FileFormats.to_str(self.format))
            h = _escape_for_latex('Format of the file valid settings are %s' % self._getValidFormatNames())
            out.write("\index{%s!%smode}\item[\\bftt  %smode] %s%s\n" % (cmdname, dpref, dpref, h, d))
        if self._options:
            self._options._write_latex_item(out, pref, cmdname)

    def _current_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying current value.
        Assumes that PhycasTablePrinter._reset_term_width has been called
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = ""
        v = None
        fn = self._getFilename()
        if fn:
            if self._is_active:
                v = str_value_for_user(fn)
            else:
                v = "<silenced>"
        opts_help = [PhycasTablePrinter.format_help("%s" %pref, v)]
        if fn:
            if self.prefix:
                opts_help.append(PhycasTablePrinter.format_help("%sprefix" % dpref, str_value_for_user(self.prefix)))
            if self.filename:
                opts_help.append(PhycasTablePrinter.format_help("%sfilename" % dpref, str_value_for_user(self.filename)))
            opts_help.append(PhycasTablePrinter.format_help("%smode" % dpref, str(self.mode)))
            if len(self._valid_formats) > 1:
                opts_help.append(PhycasTablePrinter.format_help("%s.format" % dpref, FileFormats.to_str(self.format)))
        if self._options:
            opts_help.extend(self._options._current_str_list(pref))
        return opts_help

    def __setattr__(self, name, v):
        if name.startswith("_"):
            self.__dict__[name] = v
        elif name == "filename":
            if v:
                self.__dict__["prefix"] = None
            self.__dict__["filename"] = v
        elif name == "prefix":
            if v:
                self.__dict__["filename"] = None
            self.__dict__["prefix"] = v
        elif name == "mode":
            if not isinstance(v, ExistingFileBehavior):
                raise ValueError('Valid settings of the "mode" attribute are %s'  % self._getValidModeNames())
            self.__dict__["mode"] = v
        elif name == "format":
            if not (v  in self._valid_formats):
                raise ValueError('Valid settings of the "format" attribute are %s' % self._getValidFormatNames())
            self.__dict__["format"] = v
        elif self._options and name in self._options:
            self._options._set_opt(name.lower(), v)
        else:
            raise AttributeError('Output specifier has no "%s" attribute"' % name)

    def _getValidFormatNames(self):
        v = self._valid_formats
        n = [FileFormats.to_str(i) for i in v]
        lv = len(v)
        if lv == 0:
            return None
        if lv == 1:
            return n[0]
        if lv == 2:
            return "%s or %s" % (n[0], n[1])
        return "%s, or %s" % (", ".join(n[0:-1]), n[-1])

    def set(self, v):
        s = str(v)
        if s == v:
            self.filename = s
        else:
            raise ValueError("An output setting can (currently) only be set to a string that is taken to be the filename")

    def _getValidModeNames(self):
        return 'REPLACE, APPEND, ADD_NUMBER, or SKIP'

    def __str__(self):
        return "Filestream to %s" % self._getFilename()
        #PhycasTablePrinter._reset_term_width()
        #return "\n".join(self._help_str_list())

    def _getFilename(self):
        if self._opened_filename:
            return self._opened_filename
        return self._calcFilename()

    def _calcFilename(self):
        pref, suffix = None, None
        if self.filename:
            fn = self.filename
            pref, suffix = os.path.splitext(fn)
        elif self.prefix:
            pref, suffix = self.prefix, self._getSuffix()
        if pref:
            fn = "%s%s" % (pref, suffix)
            if isinstance(self.mode, AddNumberExistingFileBehavior):
                curr_index = 1
                while os.path.exists(fn):
                    fn = "%s%d%s" % (pref, curr_index ,suffix)
                    curr_index += 1
            return fn
        return None

    def _getUnchangedFilname(self):
        if isinstance(self.mode, AddNumberExistingFileBehavior):
            m = self.mode
            try:
                self.mode = REPLACE
                fn = self._calcFilename()
            finally:
                self.mode = m
            return fn
        else:
            return self._getFilename()

    def open(self, out, as_log_file=False):
        # Do not open the file if user has deactivated
        if not self.__bool__():
            return None

        # Do not open the file if user has chosen mode = SKIP
        if isinstance(self.mode, SkipExistingFileBehavior):
            return None

        self._opened_filename = self._getFilename()
        self._prexisting = os.path.exists(self._opened_filename)
        m =  isinstance(self.mode, AppendExistingFileBehavior) and "a" or "w"
        self._opened_file = open(self._opened_filename, m)
        if out is not None:
            if self._prexisting:
                if m == "a":
                    out.info("Appending to %s" %self._opened_filename)
                else:
                    out.info("Replacing %s" %self._opened_filename)
            else:
                uf = self._getUnchangedFilname()
                if uf != self._opened_filename:
                    out.info("The file %s already exists, using %s instead" % (uf, self._opened_filename))
        self._opened_file_is_log = as_log_file
        if as_log_file:
            self._opened_file_is_log_in = out
            out.add_mirror(self._opened_file)
        else:
            self._opened_file_is_log_in = None
        return self._opened_file

    def openAsLog(self, out):
        return self.open(out, True)

    def close(self):
        if self._opened_file:
            if self._opened_file_is_log_in is not None:
                self._opened_file_is_log_in.remove_mirror(self._opened_file)
            self._opened_file.close()
        self._opened_file = None
        self._opened_filename = None
        self._opened_file_is_log_in = None

    def replaceMode(self):
        return isinstance(self.mode, ReplaceExistingFileBehavior)

    def appendMode(self):
        return isinstance(self.mode, AppendExistingFileBehavior)

    def addNumberMode(self):
        return isinstance(self.mode, AddNumberExistingFileBehavior)

    def skipMode(self):
        return isinstance(self.mode, SkipExistingFileBehavior)

    def sameBaseName(self, other):
        if self.filename:
            return self.filename == other.filename
        return self.prefix and (self.prefix == other.prefix)

class TextOutputSpec(FileOutputSpec):
    def __init__(self, prefix="", suffix=".txt", help_str="", filename=None):
        FileOutputSpec.__init__(self, prefix, help_str, filename)
        self._valid_formats = [FileFormats.RAW_TEXT]
        self.__dict__["suffix"] = suffix
    def _getSuffix(self):
        return self.suffix

class ROutputSpec(FileOutputSpec):
    def __init__(self, prefix="", suffix=".R", help_str="", filename=None):
        FileOutputSpec.__init__(self, prefix, help_str, filename)
        self._valid_formats = [FileFormats.RAW_TEXT]
        self.__dict__["suffix"] = suffix
    def _getSuffix(self):
        return self.suffix

class BinaryOutputSpec(FileOutputSpec):
    def __init__(self, prefix="", help_str="", filename=None):
        FileOutputSpec.__init__(self, prefix, help_str, filename)
        self._valid_formats = [FileFormats.RAW_TEXT]
    def _getValidModeNames(self):
        return 'REPLACE, ADD_NUMBER, or SKIP'

class DevNullWriter(object):
    """Class that fulfills the TreeWriter and MatrixWriter interface, but does
    not write any data."""
    def writeTree(self, tree, name="", rooted=None):
        pass
    def writeCharacters(self, data_matrix):
        pass

    def finish(self):
        pass

class MultiWriter(object):
    """Class that fulfills the TreeWriter and MatrixWriter interface, by
    delegating the information to other writers."""
    def __init__(self, writers=()):
        """Should pass in a list of opened writers."""
        self.writers = writers
    def writeTree(self, tree, name="", rooted=None):
        for w in self.writers:
            w.writeTree(tree, name, rooted)
    def writeCharacters(self, data_matrix):
        for w in self.writers:
            w.writeCharacters(tdata_matrix)
    def finish(self):
        pass

class TreeOutputSpec(TextOutputSpec):
    def __init__(self, prefix="", help_str="", filename=None):
        FileOutputSpec.__init__(self, prefix, help_str, filename)
        self._valid_formats = [FileFormats.NEXUS]
        self.format = FileFormats.NEXUS
        self.__dict__["collection"] = None
        self._writer = None

    def _getSuffix(self):
        return ".tre"

    def open(self, taxa_labels, out):
        #should return a writer that implements the _writeTree and finish methods
        self._writer = None
        if self._getFilename():
            o = FileOutputSpec.open(self, out)
        if FileFormats.NEXUS == self.format:
            from phycas.readnexus._NexusReader import NexusWriter
            self._writer = NexusWriter(o, self._prexisting and self.appendMode(), taxa_labels)
        else:
            assert(False)
        if self.collection:
            writers = []
            if self._writer is not None:
                writers.append(self._writer)
            if isinstance(self.collection, list):
                writers.extend(self.collection)
            else:
                writers.append(self.collection)
            self._writer = MultiWriter(writers)
        return self._writer

    def close(self):
        if self._writer:
            self._writer.finish()
            self._writer = None
        self._opened_pdfgen = None
        FileOutputSpec.close(self)


class PDFOutputSpec(BinaryOutputSpec):
    def _getSuffix(self):
        return ".pdf"
    def __init__(self, prefix="", help_str="", filename=None, mode=ADD_NUMBER):
        BinaryOutputSpec.__init__(self, prefix, help_str, filename)
        self.mode = mode
        self._valid_formats = [FileFormats.PDF]
        self.format = FileFormats.PDF
        self._opened_pdfgen = None

    def open(self, w, h, out):
        FileOutputSpec.open(self, out)
        pdf = PDFGenerator(w, h)
        self._opened_pdfgen = pdf
        return pdf

    def close(self):
        if self._opened_pdfgen:
            assert(self._opened_file)
            self._opened_pdfgen.saveDocument(stream=self._opened_file)
        self._opened_pdfgen = None
        FileOutputSpec.close(self)

class PhycasCommandOutputOptions(object):
    """The PhycasCommandOutputOptions is a very simple class written to
    provide a standardized place for command output settings to reside so
    that users can clearly see what the output options are (without cluttering
    command settings with lots of settings such as tree_out_replace,
    tree_out_append...
    """
    def __init__(self, verbosity_level=None):
        self.__dict__["level"] = 0
        self.__dict__["_cached_level"] = 0
        if verbosity_level is None:
            self.level = getDefaultOutFilter()
        else:
            self.level = verbosity_level
        self.__dict__["_help_order"] = []
        self.__dict__["_outputter"] = None
        self.__dict__["_stream"] = None

    def __deepcopy__(self, memo):
        c = memo.get(self)
        if not c is None:
            return c
        o = PhycasCommandOutputOptions(self.level)
        for k, v in self.__dict__.iteritems():
            if v is None:
                o.__dict__[k] = None
            else:
                o.__dict__[k] = copy.deepcopy(v, memo)
        memo[self] = o
        return o

    def __setattr__(self, name, value):
        if name == "_outputter" or name == "_stream":
            self.__dict__[name] = value
            return
        o = self.__dict__.get(name)
        if o is not None:
            if name == "level":
                if value < OutFilter.DEBUGGING or value > OutFilter.SILENT:
                    raise ValueError("The verbosity level must be set to one of the following:\n    OutFilter.%s" % "\n    OutFilter.".join(["DEBUGGING", "VERBOSE", "NORMAL", "WARNINGS", "ERRORS", "SILENT"]))
                self.__dict__["level"] = value
                self.__dict__["_cached_level"] = value
            else:
                isintarg = isinstance(value, int) or isinstance(value, long)
                turning_off = ((value is None) or (value is False) or (isintarg and value == 0))
                turning_on = ((value is True) or (isintarg and value != 0))
                if turning_off:     # POL deleted "and (outp is not None)"
                    o._silence()
                elif turning_on:    # POL deleted  "and (outp is not None)"
                    o._activate()
                else:
                    o.set(value)
        else:
            if name in OutFilter._names:
                i = OutFilter._names.index(name)
                if value:
                    if self.level > value:
                        self.__dict__["level"] = value
                        self.__dict__["_cached_level"] = value
                elif self.level < value + 1:
                    self.__dict__["level"] = value
                    self.__dict__["_cached_level"] = value
            else:
                raise AttributeError("%s has no attribute %s" % (self.__class__.__name__, name))

    def _silence(self):
        """Silence the output from the associated function.

        Somewhat counterintuitively, this function moves the level down
        to OutFilter.ERRORS rather than OutFilter.SILENT.
        The OutFilter.SILENT level is intended only for developers
        It can be set by (for an instance called `out`):
            out._silence()
            out.level = OutFilter.SILENT
        """
        self.__dict__["level"] = OutFilter.ERRORS
        for n in self.__dict__["_help_order"]:
            self.__dict__[n]._silence()

    def _activate(self):
        self.__dict__["level"] = self._cached_level
        for n in self.__dict__["_help_order"]:
            self.__dict__[n]._silence()

    def _help_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = pref
        opts_help = [PhycasTablePrinter.format_help("%slevel" %dpref, "Controls the amount of output (verbosity) of the command")]
        for n in self.__dict__["_help_order"]:
            a = self.__dict__[n]
            opts_help.append("")
            opts_help.extend(a._help_str_list("%s%s" % (dpref, n)))
        return opts_help

    def _latex(self, pref=""):
        """Generates LaTeX description list for inclusion in the Phycas manual"""
        latex = StringIO()
        self._write_latex(latex, pref)
        return latex.getvalue()

    def _write_latex(self, out, pref="", cmdname=""):
        out.write('\\begin{description}\n')
        out.write('\index{%s!%s.level}\item[\\bftt  %s.level] Controls the amount of output (verbosity) of the command (default: OutFilter.%s)\n' % (cmdname, pref, pref, OutFilter.to_str(self.level)))
        for n in self.__dict__["_help_order"]:
            a = self.__dict__[n]
            which = '%s.%s' % (pref,n)
            a._write_latex_item(out, which, cmdname)
        out.write('\\end{description}\n')

    def _current_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = pref
        s = "OutFilter." + OutFilter.to_str(self.level)
        opts_help = [PhycasTablePrinter.format_help("%slevel" %dpref, s)]
        for n in self.__dict__["_help_order"]:
            a = self.__dict__[n]
            opts_help.append("")
            opts_help.extend(a._current_str_list("%s%s" % (dpref, n)))
        return opts_help

    def __str__(self):
       PhycasTablePrinter._reset_term_width()
       return "\n".join(self._help_str_list())

    def getStdOutputter(self):
        """Returns self._outputter, or a filter around self._stream, or (if
        both of those attributes are None) a filter around the default
        output stream.
        """
        if self._outputter:
            return self._outputter
        s = self._stream
        if s is None:
            s = getDefaultOutputStream()
        return OutputFilter(self.level, s)

class PhycasCmdOpts(object):

    def __init__(self, command=None, args=None):
        self.__dict__["_current"] = {}
        self.__dict__["_transformer"] = {}
        self.__dict__["_help_info"] = {}
        self.__dict__["_optionsInOrder"] = None
        self.__dict__["_command"] = None
        self.__dict__["_unchecked"] = set()
        if command and args:
            self._initialize(command, args)

    def _initialize(self, command=None, args=None):
        self._current = {}
        self._transformer = {}
        self._help_info = {}
        self._optionsInOrder = args
        self._command = command
        if args is not None:
            for opt in args:
                if len(opt) == 3:
                    name, default, help_str = opt
                    transf = None
                else:
                    name, default, help_str, transf = opt
                name = name.lower()
                self._help_info[name] =  [default, help_str]
                self._current[name] = default
                if transf is not None:
                    self._transformer[name] = transf
                self._command.__dict__[name] = self._current[name]

    def _latex(self):
        """Generates LaTeX description list for inclusion in the Phycas manual"""
        latex = StringIO()
        self._write_latex(latex)
        return latex.getvalue()

    def _write_latex(self, out, pref="", cmdname=""):
        out.write('\\begin{description}\n')
        self._write_latex_item(out, pref, cmdname)
        out.write('\\end{description}\n')

    def _write_latex_item(self, out, pref="", cmdname=""):
        "writes the options as items in LaTeX"
        if pref:
            dpref = _escape_for_latex(pref + ".")
        else:
            dpref = ""
        for o in self._optionsInOrder:
            name = o[0]
            o_default_value = o[1]
            default_value = None
            c = o_default_value.__class__.__name__
            if c == 'str':
                default_value = "{\\bftt '%s'}" % o_default_value
            elif c == 'Model':
                default_value = "predefined model object"
            else:
                default_value = '{\\bftt  %s}' % str(o_default_value)
            descrip = o[2]
            out.write('\index{%s!%s}\item[\\bftt  %s%s] %s' % (cmdname, _escape_for_latex(name), dpref, _escape_for_latex(name), _escape_for_latex(descrip)))
            if default_value:
                out.write(' (default: %s)\n' % _escape_for_latex(default_value))
            else:
                out.write('\n')

    def _help_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = ""
        opts_help = []
        if self._optionsInOrder is not None:
            for i in self._optionsInOrder:
                oc_name = i[0]
                name = oc_name.lower()
                n = pref and "%s%s" % (dpref, oc_name) or oc_name
                s = PhycasTablePrinter.format_help(n, i[2])
                if help_double_space and i != self._optionsInOrder[0]:
                    s = '\n%s' % s
                opts_help.append(s)
        return opts_help

    def _current_str_list(self, pref=""):
        """Generates a list of strings formatted for displaying help
        Assumes that PhycasTablePrinter._reset_term_width has been called
        more recently than the last terminal width change)."""
        if pref:
            dpref = pref + "."
        else:
            dpref = ""
        opts_help = []
        if self._optionsInOrder is not None:
            for i in self._optionsInOrder:
                oc_name = i[0]
                name = oc_name.lower()
                n = pref and "%s%s" % (dpref, oc_name) or oc_name
                s = PhycasTablePrinter.format_help(n, str_value_for_user(self._current[name]))
                if current_double_space and i != self._optionsInOrder[0]:
                    s = '\n%s' % s
                opts_help.append(s)
        return opts_help

    def __str__(self):
        PhycasTablePrinter._reset_term_width()
        return "\n".join(self._help_str_list())

    def _set_opt(self, name, value, force_valid=True):
        if name in self._current:
            transf = self._transformer.get(name)
            if transf is None:
                self._current[name] = value
            else:
                try:
                    self._current[name] = transf(self, value)
                except:
                    if force_valid:
                        try:
                            c = ""
                            c = transf.get_description()
                        except:
                            pass
                        if value.__class__.__name__ == 'str':
                            raise ValueError("'%s' is not a valid value for %s. %s" % (value, name, c))
                        else:
                            raise ValueError("%s is not a valid value for %s. %s" % (value, name, c))
                    else:
                        self._current[name] = value
            self._command.__dict__[name] = self._current[name]
        else:
            raise AttributeError("%s does not contain an attribute %s" % (self._command.__class__.__name__, name))

    def __contains__(self, name):
        return name in self._current

    def __getattr__(self, name):
        if name in self._current:
            return self._current[name]
        raise AttributeError("%s does not contain an attribute %s" % (self.__class__.__name__, name))

    def __setattr__(self, name, value):
        if name in self.__dict__:
            self.__dict__[name] = value
        else:
            self._set_opt(name, value)

    def check_unchecked(self):
        for name in self._unchecked:
            value = self._current[name]
            try:
                transf = self._transformer.get(name)
                if not transf is None:
                    self._current[name] = transf(self, value)
            except:
                pass # it is concievable that some errors here are not errors once all of the other opts have been transformed, so we set options twice
        for name in self._unchecked:
            value = self._current[name]
            self._set_opt(name, value)
        self._unchecked.clear()

    def set_unchecked(self, name, value):
        "Sets the attribute `key` to the `value`"
        if name in self._current:
            self._current[name] = value
            self._command.__dict__[name] = value
            self._unchecked.add(name)
        else:
            raise AttributeError("%s does not contain an attribute %s" % (self._command.__class__.__name__, name))

def TreeSourceValidate(opts, v):
    if v is None:
        return None
    if isinstance(v, TreeCollection):
        return v
    return TreeCollection(filename=v)

def DataSourceValidate(opts, v):
    if v is None:
        return None
    if isinstance(v, DataSource):
        return v
    return DataSource(filename=v)

def BoolArgValidate(opts, v):
    return bool(v)

def FileExistsValidate(opts, v):
    import os.path
    if not os.path.exists(v):
        raise ValueError("Expecting '%s' to be an existing file" % v)
    return v

def FileExistsOrNoneValidate(opts, v):
    if v is None:
        return None
    import os.path
    if not os.path.exists(v):
        raise ValueError("Expecting '%s' to either be None or an existing file" % v)
    return v

def FileDoesNotExistOrNoneValidate(opts, v):
    if v is None:
        return None
    import os.path
    if os.path.exists(v):
        raise ValueError("Expecting '%s' to either be None or the name of a file that does not already exist" % v)
    return v

class EnumArgValidate(object):
    def __init__(self, valid_values):
        l = []
        for i in valid_values:
            if isinstance(i, str):
                l.append(i.lower())
            else:
                l.append(i)
        self.valid_tuple = tuple(l)
        self.none_is_valid_value = False
        for_display = []
        for v in self.valid_tuple:
            if v is None:
                for_display.append('None')
                self.none_is_valid_value = True
            elif v.__class__.__name__ == 'str':
                for_display.append("'%s'" % v)
            else:
                for_display.append("%s" % str(v))
        self.display_str = ','.join(for_display)
    def __call__(self, opts, v):
        if v is None:
            if not self.none_is_valid_value:
                raise ValueError('Value cannot be None')
        else:
            if not isinstance(v, str):
                raise ValueError('Value must be a Python string')
            v = v.lower()
            if v not in self.valid_tuple:
                raise ValueError("Value must be one of the following: %s" % self.display_str)
        return v
    def get_description(self):
        return "Must be one of the following: %s" % self.display_str

class IntArgValidate(object):
    def __init__(self, min=None, max=None):
        self.min = min
        self.max = max
    def __call__(self, opts, v):
        iv = int(v)
        if (self.min is not None) and (self.min > iv):
            raise ValueError("Value must be >= %s" % str(self.min))
        if (self.max is not None) and (self.max < iv):
            raise ValueError("Value must be <= %s" % str(self.min))
        return iv
    def get_description(self):
        if self.min is None:
            if self.max is None:
                return "Must be an integer."
            else:
                return "Must be an integer no greater than %d" % self.max
        else:
            if self.max is None:
                return "Must be an integer no less than %d" % self.min
            else:
                return "Must be an integer in the range from %d to %d " % (self.min, self.max)

class FloatArgValidate(object):
    def __init__(self, min=None, max=None, lessthan=None, greaterthan=None):
        self.min = min
        self.max = max
        self.lessthan = lessthan
        self.greaterthan = greaterthan
    def __call__(self, opts, v):
        fv = float(v)
        if (self.greaterthan is not None) and (self.greaterthan >= fv):
            raise ValueError("Value must be > %s" % str(self.min))
        elif (self.min is not None) and (self.min > fv):
            raise ValueError("Value must be >= %s" % str(self.min))
        if (self.lessthan is not None) and (self.lessthan <= fv):
            raise ValueError("Value must be < %s" % str(self.min))
        elif (self.max is not None) and (self.max < fv):
            raise ValueError("Value must be <= %s" % str(self.min))
        return fv
    def get_description(self):
        s = 'Must be a real number'
        if (self.min is None) and (self.max is None) and (self.lessthan is None) and (self.greaterthan is None):
            s += '.'
        else:
            s += ' '
            if self.greaterthan is not None:
                s += 'greater than %f' % self.greaterthan
            elif self.min is not None:
                s += 'greater than or equal to %f' % self.min
            if self.lessthan is not None:
                s += ' and less than %f' % self.lessthan
            elif self.max is not None:
                s += ' and less than or equal to %f' % self.max
        return s

class ProbArgValidate(FloatArgValidate):
    def __init__(self):
        FloatArgValidate.__init__(self, min=0.0, max=1.0)

def LotArgValidate(opts, v):
    from phycas.probdist._Lot import Lot
    if isinstance(v, Lot):
        return v
    raise ValueError("Mest be a probdist.Lot instance")

class PhycasCommandHelp(object):
    def __init__(self, command, cmd_name, cmd_descrip):
        self.command = command
        self.cmd_name = cmd_name
        self.cmd_descrip = cmd_descrip

    def __str__(self):
        if current_follows_help:
            return "\n".join([self.explain(), "\n", self.current()])
        else:
            return "\n" + self.explain()

    def explain(self):
        PhycasTablePrinter._reset_term_width()
        command = self.command
        opts = command.__dict__["_options"]
        out = command.__dict__["out"]
        d = PhycasTablePrinter.get_help_divider()
        a = [self.cmd_name]
        a.extend(PhycasTablePrinter._text_wrapper.wrap(self.cmd_descrip))
        o = opts._help_str_list()
        o.sort()
        h = PhycasTablePrinter.format_help("Attribute", "Explanation")
        if len(o) > 0:
            a.extend(["\nAvailable input options:", h, d])
            a.extend(o)
            a.extend([d, ""])
        if out is not None:
            a.extend(["Available output settings:", h, d])
            a.extend(out._help_str_list("out"))
        return "\n".join(a)

    def current(self):
        PhycasTablePrinter._reset_term_width()
        command = self.command
        a = []
        opts = command.__dict__["_options"]
        out = command.__dict__["out"]
        d = PhycasTablePrinter.get_help_divider()
        h = PhycasTablePrinter.format_help("Attribute", "Current Value")
        o = opts._current_str_list()
        o.sort()
        if len(o) > 0:
            a.extend(["\nCurrent %s input settings:" % self.cmd_name , h, d])
            a.extend(o)
            a.extend([d, ""])
        if out is not None:
            a.extend(["Current %s  output settings:" % self.cmd_name, h, d])
            a.extend(out._current_str_list("out"))
        return "\n".join(a)

    def __call__(self):
        """Returns the result of __str__, so that users can omit the () on the
        and simply type "cmd.help" to see the help message."""
        print(str(self))

    def __repr__(self):
        return str(self)

class PhycasCurrentValuesHelper:
    def __init__(self, helper):
        self.helper = helper
    def __str__(self):
        return self.helper.current()
    def __call__(self):
        """Returns the result of __str__, so that users can omit the () on the
        and simply type "cmd.help" to see the help message."""
        print(str(self))

    def __repr__(self):
        return str(self)

class PhycasManualGenerator(object):
    def __init__(self, command, cmd_name, cmd_descrip):
        self.command = command
        self.cmd_name = cmd_name
        self.cmd_descrip = cmd_descrip

    def manual(self):
        f = open(self.command.help.cmd_name+'.tex', 'w')
        opts = self.command.__dict__["_options"]
        if opts:
            f.write("\paragraph[Input options]{Input options:}")
            opts._write_latex(f, '', self.command.help.cmd_name)
        out = self.command.__dict__['out']
        if out:
            f.write("\paragraph[Output options]{Output options:}")
            out._write_latex(f, 'out', self.command.help.cmd_name)
        f.close()

    def __str__(self):
        """Allows users to avoid using function call syntax"""
        self.manual()

    def __call__(self):
        """Creates a LaTeX file whose name is <command name>.tex that can be included in the Phycas manual"""
        self.manual()

class PhycasCommand(object):

    def __init__(self, option_defs, cmd_name, cmd_descrip, output_options=None):
        # the roundabout way of initializing PhycasCmdOpts is needed because
        #   _options must be in the __dict__ before the setattr is called
        #   on a PhycasCommand
        object.__init__(self)
        c = self.__class__
        if c is not PhycasCommand:
            PhycasHelp._phycas_cmd_classes.add(c)
        o = PhycasCmdOpts()
        self.__dict__["_options"] = o
        o._initialize(self, option_defs)
        self.__dict__["out"] = output_options
        h = PhycasCommandHelp(self, cmd_name, cmd_descrip)
        self.__dict__["help"] = h
        self.__dict__["current"] = PhycasCurrentValuesHelper(h)
        m = PhycasManualGenerator(self, cmd_name, cmd_descrip)
        self.__dict__["curr"] = self.__dict__["current"] # allow abbreviation for current
        self.__dict__["manual"] = m

    def hidden():
        """
        Override this function to return True if you want to keep a class's name from being
        displayed in the main phycas help display.
        """
        return False
    hidden = staticmethod(hidden)

    def __deepcopy__(self, memo):
        c = memo.get(self)
        if not c is None:
            return c
        opts = self._options
        opts_copy = opts._optionsInOrder
        h = self.help
        # this creates a new command instance, but some of the settings may
        #   reflect the initial unset state (those in which the values have
        #   been rebound in the `self` instance, thus causing a discrepancy
        #   between the initial value in opts._optionsInOrder and the current
        #   value in opts._current.  We fix this below
        out = self.out
        if out:
            out = copy.deepcopy(out, memo)
        c = PhycasCommand(opts_copy, h.cmd_name, h.cmd_descrip, out)
        # update the settings to a copy of the current settings
        opts_copy = c._options
        managed_dict = opts_copy._current
        for k, v in opts._current.iteritems():
            opts_copy.set_unchecked(k, v)
        # update the ad-hoc unmanaged attributes
        for k, v in self.__dict__.iteritems():
            if not k in managed_dict:
                if not k  in ["_options", "help", "current", "manual", "out"]:
                    c.__dict__[k] = copy.deepcopy(v, memo)
        memo[self] = c
        return c

    def __setattr__(self, name, value):
        nl = name.lower()
        o = self.__dict__["_options"]
        if nl in o:
            o._set_opt(nl, value)
        elif nl == "out" and "out" in self.__dict__:
            outp = self.__dict__["out"]
            isintarg = isinstance(value, int) or isinstance(value, long)
            turning_off = ((value is None) or (value is False) or (isintarg and value == 0))
            turning_on = ((value is True) or (isintarg and value != 0))
            if nl == "out":
                if turning_off and (outp is not None):
                    outp._silence()
                elif turning_on and (outp is not None):
                    outp._activate()
                else:
                    self.__dict__["out"] = value
        elif name in self.__dict__:
            self.__dict__[name] = value
        else:
            raise AttributeError("%s has no attribute %s" % (self.__class__.__name__, name))

    def __getattr__(self, name):
        nl = name.lower()
        try:
            return self.__dict__[nl]
        except:
            raise AttributeError("%s has no attribute %s" % (self.__class__.__name__, name))

    def set(self, **kwargs):
        "Sets the attributes of a command using keyword arguments"
        old = {}
        o = self._options
        try:
            for key, value in kwargs.iteritems():
                old[key] = getattr(self, key)
                o.set_unchecked(key, value)
            o.check_unchecked()
        except:
            for key, value in old.iteritems():
                o.set_unchecked(key, value)
            raise

    def _brief_str(self):
        global _use_instance_names
        v = None
        if _use_instance_names:
            v = name_of(self)
        if v:
            return v
        v = "with id=%d" % id(self)
        return "The %s command instance %s" % (self.help.cmd_name, v)

    def _getRNGOptions():
        return [("rng", None, "A pseudo-random number generator instance", LotArgValidate),
                ("random_seed", 0, "Determines the random number seed used; specify 0 to generate seed automatically from system clock", IntArgValidate(min=0))]
    _getRNGOptions = staticmethod(_getRNGOptions)
